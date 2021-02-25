"""
AMC pipeline to process MiSeq data for diagnostics.
Martin Haagmans https://github.com/zaag
February 2017
"""
import os
import csv
import sys
import glob
import json
import shutil
import sqlite3
import subprocess

import vcf
import matplotlib
import xlsxwriter
import numpy as np
import pandas as pd

from ngsscriptlibrary.parsing import *

from ngsscriptlibrary import metrics2db
from ngsscriptlibrary import summary2db
from ngsscriptlibrary import sangers2db
from ngsscriptlibrary import snpcheck2db
from ngsscriptlibrary import mean_std_2db
from ngsscriptlibrary import sample_in_db
from ngsscriptlibrary import get_patient_info
from ngsscriptlibrary import get_rs_gpos_dict
from ngsscriptlibrary import compare_snpchecks
from ngsscriptlibrary import get_file_locations
from ngsscriptlibrary import samplesheetinfo2db
from ngsscriptlibrary import base_perc_reads2db
from ngsscriptlibrary import get_snpcheck_serie
from ngsscriptlibrary import perc_target_covered2db
from ngsscriptlibrary import riskscore_and_genotypes_2db
from ngsscriptlibrary import parse_sangers_for_seriereport

matplotlib.use('Agg')

# Define variables from config file
GATK = config["GATK"]
REF = config["REF"]
DBSNP = config["DBSNP"]
PICARD = config["PICARD"]
TARGETREPO = config["TARGET"]
DBREPO = config['DB']
JAVA = config['JAVA']

# Define files
TARGETDB = os.path.join(TARGETREPO, 'varia', 'captures.sqlite')
SNPCHECK = os.path.join(TARGETREPO, 'varia', 'NGS-SNPcheck.vcf')
METRICSDB = os.path.join(DBREPO, 'metrics.sqlite')
SAMPLEDB = os.path.join(DBREPO, 'samplesheets.sqlite')
PATIENTDB = os.path.join(DBREPO, 'patientinfo.sqlite')

SAMPLESHEET = 'SampleSheet.csv'
STANDAARDFRAGMENTEN = 'standaardfragmenten.xlsx'

for fn in [TARGETDB, SNPCHECK, METRICSDB, SAMPLEDB, PATIENTDB, SAMPLESHEET]:
    assert os.path.isfile(fn), '{} bestaat niet'.format(fn)

def add_patient_info(todo, serie, db):
    for sample in todo.keys():
        sex, ff, dob = get_patient_info(sample, serie, db)
        if sex == 'M':
            sex = 'Man'
        elif sex == 'V':
            sex = 'Vrouw'
        todo[sample]['geslacht'] = sex
        todo[sample]['FFnummer'] = ff
        todo[sample]['geboortedatum'] = dob
    return todo

def create_recalinput_file(samples):
    """Create a file with a list of BAM-files (1 per sample) to use as input
    for the recalibration step in the pipeline.
    """
    with open('recalinput.list', 'w') as f:
        for sample in samples:
            f.write('output/{}.DM.bam\n'.format(sample))


# Get todo list for samples
input_dict = parse_samplesheet_for_pipeline(SAMPLESHEET, TARGETDB, ['Amplicon', 'Amplicons'])
input_dict = get_file_locations(input_dict, TARGETREPO)

captures = get_captures(input_dict)
pakketten = get_pakketten(input_dict)
samples = [s for s in input_dict.keys() if not input_dict[s]['amplicon']]
serie = [input_dict[_]['serie'] for _ in input_dict.keys()][0]
input_dict = add_patient_info(input_dict, serie, PATIENTDB)
standard_frags_no_initials = list()

# Rules for snakemake
onsuccess:
    shell("if [ -f recalinput.list ]; then rm -f recalinput.list ; fi ;")
    shell("if [ -z '$(ls -A tempfiles)' ]; then rmdir tempfiles ; fi ;")
    shell("if [ -z '$(ls -A alignments)'  ]; then rmdir alignments ; fi ;")

rule all:
    input:
        expand("output/{sample}.filtered.vcf", sample=samples),
        expand("output/{sample}.xlsx", sample=samples),
        expand("output/{sample}.DM.bam", sample=samples),
        "output/input.json",
        "output/MS{}_report.xlsx".format(serie),
        "output/MS{}_QC.pdf".format(serie)


rule prepare:
    input:
        {SAMPLESHEET}
    output:
        expand("reads/{sample}.R1.fastq.gz", sample=samples),
        expand("reads/{sample}.R2.fastq.gz", sample=samples)
    run:
        for sample in samples:
            samplesheetinfo2db(input_dict[sample], sample, SAMPLEDB)
            R1 = glob.glob('{s}_*R1*.gz'.format(s=sample))
            R2 = glob.glob('{s}_*R2*.gz'.format(s=sample))
            if R1 and R2:
                os.rename(R1[0], 'reads/{}.R1.fastq.gz'.format(sample))
                os.rename(R2[0], 'reads/{}.R2.fastq.gz'.format(sample))


rule dump_input_dict:
    input:
        {SAMPLESHEET}
    output:
        "output/input.json"
    run:
        with open("{}".format(output), 'w') as fp:
            json.dump(input_dict, fp, indent=4)


rule countbases:
    input:
        fq = "reads/{sample}.R1.fastq.gz"
    output:
        "tempfiles/{sample}.basepercentages.done"
    message:
        "Calculating per base percentage"
    run:
        if sample_in_db(wildcards.sample, serie, METRICSDB, 'basepercentages'):
            shell("touch {output}")
        elif not sample_in_db(wildcards.sample, serie, METRICSDB, 'basepercentages'):
            counts = get_basecounts(input.fq)
            base_perc_reads2db(wildcards.sample, serie, counts, METRICSDB)
            shell("touch {output}")


rule mapreads:
    input:
        "reads/{sample}.R1.fastq.gz",
        "reads/{sample}.R2.fastq.gz"
    output:
        temp("alignments/{sample}.sorted.bam")
    threads:
        2
    log:
        "logfiles/{sample}.BWAlignment.log"
    message:
        "Aligning reads with bwa mem"
    params:
        rg = "@RG\\tID:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tPU:{sample}\\tSM:{sample}"
    shell:
        '''(bwa mem  -R '{params.rg}' -t 1 -M {REF} {input} |\
        samtools view -Shu - |\
        samtools sort -T {wildcards.sample}.tmp -O bam - > {output}) > {log}  2>&1
        '''


rule markduplicates:
    input:
        rules.mapreads.output
    output:
        bam = "output/{sample}.DM.bam",
        metrics = temp("tempfiles/{sample}.dupmark.txt")
    log:
        "logfiles/{sample}.MarkDuplicates.log"
    message:
        "Marking duplicates with picard"
    shell:
        '''{PICARD} MarkDuplicates  I={input} O={output.bam}  \
        M={output.metrics} REMOVE_DUPLICATES=FALSE CREATE_INDEX=true > {log}  2>&1
        '''


rule hsmetrics:
    input:
        rules.markduplicates.output.bam
    output:
        temp("tempfiles/{sample}.HSMetrics.txt")
    log:
        "logfiles/{sample}.HSmetrics.log"
    message:
        "Calculating HS-metrics with picard"
    run:
        target = input_dict[wildcards.sample]['picard']
        serie = input_dict[wildcards.sample]['serie']
        if sample_in_db(wildcards.sample, serie, METRICSDB, 'hsmetrics'):
            shell("touch {output}")

        elif not sample_in_db(wildcards.sample, serie, METRICSDB, 'hsmetrics'):
            shell('''{PICARD} CalculateHsMetrics R={REF} I={input} \
            O={output[0]} TI={target} BI={target} > {log} 2>&1
            ''')
            metrics = readmetrics(wildcards.sample, serie, output[0])
            metrics2db(metrics, METRICSDB, 'hsmetrics')


rule insertsizemetrics:
    input:
        rules.markduplicates.output.bam
    output:
        tabel = temp("tempfiles/{sample}.InsertSize.txt"),
        pdf = temp("tempfiles/{sample}.InsertSize.pdf")
    log:
        "logfiles/{sample}.InsertSize.log"
    message:
        "Calculating insertsize metrics"
    run:
        serie = input_dict[wildcards.sample]['serie']
        if sample_in_db(wildcards.sample, serie, METRICSDB, 'insertsize'):
            shell("touch {output.pdf}")
            shell("touch {output.tabel}")
        elif not sample_in_db(wildcards.sample, serie, METRICSDB, 'insertsize'):
            shell('''{PICARD} CollectInsertSizeMetrics R={REF} I={input} \
            O={output.tabel} HISTOGRAM_FILE={output.pdf} > {log} 2>&1
            ''')
            metrics = readmetrics(wildcards.sample, serie, output.tabel)
            metrics2db(metrics, METRICSDB, 'insertsize')


rule alignmetrics:
    input:
        rules.markduplicates.output.bam
    output:
        temp("tempfiles/{sample}.AlignmentMetrics.txt")
    log:
        "logfiles/{sample}.AlignmentMetrics.log"
    message:
        "Calculating alignment metrics"
    run:
        serie = input_dict[wildcards.sample]['serie']

        if sample_in_db(wildcards.sample, serie, METRICSDB, 'alignment'):
            shell("touch {output}")
        elif not sample_in_db(wildcards.sample, serie, METRICSDB, 'alignment'):
            shell('''{PICARD} CollectAlignmentSummaryMetrics R={REF} I={input} \
            O={output[0]}  MAX_INSERT_SIZE=500 > {log} 2>&1
            ''')
            serie = input_dict[wildcards.sample]['serie']
            metrics = readmetrics(wildcards.sample, serie, output[0])
            metrics2db(metrics, METRICSDB, 'alignment')


rule recalibrate:
    input:
        expand("output/{sample}.DM.bam", sample=samples)
    output:
        "reads/bqsr.grp"
    log:
        "logfiles/Recal.log"
    message:
        "Recalibrating quality scores with GATK"
    threads:
        10
    run:
        create_recalinput_file(samples)
        shell('''{JAVA} -Xmx10g -jar {GATK} -R {REF} -T BaseRecalibrator -knownSites {DBSNP} \
        -I recalinput.list -o {output} -nct {threads} > {log} 2>&1
        ''')


rule callableloci:
    input:
        bam = rules.markduplicates.output.bam,
        table = rules.recalibrate.output
    output:
        bed = temp("tempfiles/{sample}.bed"),
        summary = temp("tempfiles/{sample}.summary")
    message:
        "Defining callable regions with GATK"
    log:
        "logfiles/{sample}.CallableLoci.log"
    run:
        if sample_in_db(wildcards.sample, serie, METRICSDB, 'sangers'):
            shell('touch {output.summary} {output.bed}')
        else:
            target = input_dict[wildcards.sample]['sanger']
            shell('''{JAVA} -jar {GATK} -R {REF} -T CallableLoci \
            -minDepth 30 -mdflmq 30 -mmq 20 -frlmq 0.6 -mlmq 10 \
            -I {input.bam} -o {output.bed} \
            -summary {output.summary}  \
            -L {target} --BQSR {input.table} > {log} 2>&1''')
            if input_dict[wildcards.sample]['panel'] is None:
                db_targetname = input_dict[wildcards.sample]['pakket']
            else:
                db_targetname = input_dict[wildcards.sample]['panel']
            summary2db(wildcards.sample, read_summary(output.summary),
                       db_targetname, serie, METRICSDB)


rule getsangers:
    input:
        rules.callableloci.output.bed
    output:
        sangers = "tempfiles/{sample}.sangers.txt",
        tempbed = temp("tempfiles/{sample}.tmp.bed")
    message:
        "Defining sangerfragments from callable regions"
    run:
        annot_bed = input_dict[wildcards.sample]['annot']
        panel = input_dict[wildcards.sample]['panel']
        if panel is None:
            target_name = input_dict[wildcards.sample]['pakket']
        elif panel is not None:
            target_name = panel
        sangers_for_db = list()
        if sample_in_db(wildcards.sample, serie, METRICSDB, 'sangers'):
            pass
        elif panel is not None and 'OVR' in panel:
            sangers_for_db = 'Geen sangers: OVR panel'
            sangers2db(sangers_for_db, serie, wildcards.sample, target_name, METRICSDB)
        else:
            df = annotate_callables(input, annot_bed, output.tempbed)
            sangercount = 0
            for target, data in df.groupby(df.index):
                if len(data[data['callable'] != 'CALLABLE']) > 0:
                    sangercount += 1
                    gene = ' '.join(data['gene'].unique())
                    start = data[data['callable'] != 'CALLABLE']['regionstart'].min()
                    end = data[data['callable'] != 'CALLABLE']['regionend'].max()
                    sangers_for_db.append((gene, target[0], str(start), str(end),
                                           str(target[1]), str(target[2])))
            if sangercount == 0:
                sangers_for_db = 'Geen sangers: alles callable'
            sangers2db(sangers_for_db, serie,  wildcards.sample, target_name, METRICSDB)
        shell("touch {output.sangers} {output.tempbed}")


rule depthofcoverage:
    input:
        bamfile = rules.markduplicates.output.bam,
        table = rules.recalibrate.output
    output:
        "tempfiles/{sample}.DoC",
        temp("tempfiles/{sample}.DoC.sample_cumulative_coverage_counts"),
        temp("tempfiles/{sample}.DoC.sample_cumulative_coverage_proportions"),
        temp("tempfiles/{sample}.DoC.sample_interval_statistics"),
        temp("tempfiles/{sample}.DoC.sample_statistics"),
        temp("tempfiles/{sample}.DoC.sample_summary"),
        "tempfiles/{sample}.DoC.sample_interval_summary"
    message:
        "Calculating DepthOfCoverage with GATK"
    log:
        "logfiles/{sample}.DeptOfCoverage.log"
    run:
        target = input_dict[wildcards.sample]['cnvtarget']
        shell('''{JAVA} -jar {GATK} -R {REF} -T DepthOfCoverage \
        -I {input.bamfile} -o {output[0]} \
        -L {target} \
        --minBaseQuality 20 \
        --minMappingQuality 20 \
        --printBaseCounts  \
        -dels -BQSR {input.table} > {log} 2>&1
        ''')

        if input_dict[wildcards.sample]['capispakket']:
            perc_covered = calc_perc_target_covered(output[0])
            pakket = input_dict[wildcards.sample]['capture']

        elif not input_dict[wildcards.sample]['capispakket']:
            pakket_target = input_dict[wildcards.sample]['pakkettarget']
            perc_covered = calc_perc_target_covered(output[0],
                                                    filter_regions=True,
                                                    target=pakket_target)
            pakket = input_dict[wildcards.sample]['pakket']

        mean, std = mean_doc(output[0])
        target = input_dict[wildcards.sample]['capture']
        mean_std_2db(mean, std, serie, wildcards.sample, target, METRICSDB)    
        perc_target_covered2db(wildcards.sample, perc_covered, pakket, serie, METRICSDB)


rule callvariants:
    input:
        bamfile = rules.markduplicates.output.bam,
        table = rules.recalibrate.output
    output:
        temp("tempfiles/{sample}.raw.vcf"),
        temp("tempfiles/{sample}.raw.vcf.idx")
    message:
        "Calling variants with GATK's HaplotypeCaller"
    log:
        "logfiles/{sample}.HaplotypeCaller.log"
    run:
        target = input_dict[wildcards.sample]['varcal']
        shell('''{JAVA} -jar {GATK} -R {REF} -T HaplotypeCaller \
        -I {input.bamfile} -o {output[0]} -L {target} -ip 500 \
        -pairHMM VECTOR_LOGLESS_CACHING -BQSR {input.table}  > {log} 2>&1
        ''')


rule filtervariants:
    input:
        vcf = rules.callvariants.output[0],
        index = rules.callvariants.output[1]
    output:
        "output/{sample}.filtered.vcf"
    log:
        "logfiles/{sample}.VariantFiltration.log"
    message:
        "Annotating variants with GATK's QC annotations"
    shell:
        '''{JAVA} -jar {GATK} -R {REF} -T VariantFiltration \
        -V {input.vcf} -o {output} \
        --clusterWindowSize 20 --clusterSize 6 \
        --filterExpression "DP < 30 " --filterName "LowCoverage" \
        --filterExpression "QUAL < 50.0 " --filterName "LowQual" \
        --filterExpression "QD < 4.0 " --filterName "QD" \
        --filterExpression "SOR > 10.0 " --filterName "SOR" > {log} 2>&1
        '''


rule callsnpcheck:
    input:
        rules.markduplicates.output.bam
    output:
        vcf = "tempfiles/{sample}.snpcheck.vcf",
        index = "tempfiles/{sample}.snpcheck.vcf.idx"
    message:
        "Calling SNP checks sample"
    log:
        "logfiles/{sample}.CallSNPcheck.log"
    run:
        if sample_in_db(wildcards.sample, serie, METRICSDB, 'snpcheck'):
            shell('touch {output.vcf} {output.index}')  
        elif not sample_in_db(wildcards.sample, serie, METRICSDB, 'snpcheck'):
            shell('''{JAVA} -jar {GATK} -R {REF} -T HaplotypeCaller \
                -I {input} -o {output.vcf} -L {SNPCHECK} -alleles {SNPCHECK} \
                -gt_mode GENOTYPE_GIVEN_ALLELES -GQB 0 > {log} 2>&1
                ''')
            


rule comparesnpcheckssample:
    input:
        vcf = rules.callsnpcheck.output.vcf,
        alt = "SNPcheck/{sample}.qpcrsnpcheck"
    output:
        temp("SNPcheck/{sample}.snpcheck.done.txt")
    message:
        "Comparing SNP checks sample"
    log:
        "logfiles/{sample}.CompareSNPcheck.log"
    run:
        if sample_in_db(wildcards.sample, serie, METRICSDB, 'snpcheck'):
            shell('touch {output}')
        else:
            store_in_db = dict()
            store_in_db['NGS'] = parse_ngssnpcheck(input.vcf)
            store_in_db['ALT'] = parse_taqman(input.alt)
            store_in_db['COMP'] = compare_snpchecks(store_in_db['ALT'], store_in_db['NGS'])
            store_in_db_json = json.dumps(store_in_db)
            snpcheck2db(wildcards.sample, serie, store_in_db_json, METRICSDB)
            shell('touch {output}')


rule addcnvdata:
    input:
        rules.depthofcoverage.output[-1]
    output:
        temp("tempfiles/{sample}.cnv_data_added.txt")
    message:
        "Add CNV data"
    log:
        "logfiles/{sample}.CNV_add_data.log"
    run:
        if input_dict[wildcards.sample]['cnvscreening']:
            target = input_dict[wildcards.sample]['capture']
            shell("CNV -c {target} -s {serie} -n {input} --sample {wildcards.sample} --addonly > {log} 2>&1")
            shell("touch {output}")
        else:
            shell("touch {output}")


rule cnvdetection:
    input:
        expand(rules.addcnvdata.output, sample=samples)
    output:
        "tempfiles/cnvdetection.txt"
    message:
        "CNV detection"
    log:
        "logfiles/CNV_detector.log"
    run:
        for capture in captures:
            screening = False
            for s in samples:
                workdir = os.path.join(os.getcwd(), 'output', 'CNV_{}'.format(capture))
                if input_dict[s]['capture'] == capture:
                    screening = input_dict[s]['cnvscreening']
            if screening:
                os.makedirs(workdir, exist_ok=True)
                shell("CNV -c {capture} -s {serie} -o {workdir} > {log} 2>&1")
        shell("touch {output}")


rule callriskscore:
    input:
        rules.markduplicates.output.bam
    output:
        vcf = "tempfiles/{sample}.riskscore.vcf",
        index = "tempfiles/{sample}.riskscore.vcf.idx"
    message:
        "Calling risk score loci"
    log:
        "logfiles/{sample}.CallRiskScore.log"
    run:
        in_database = sample_in_db(wildcards.sample, serie, METRICSDB, 'riskscore')
        if in_database or not input_dict[wildcards.sample]['riskscore']:
            shell('touch {output.vcf} {output.index}')
        elif not in_database and input_dict[wildcards.sample]['riskscore']:
            riskcore_vf = input_dict[wildcards.sample]['riskscorevcf']
            shell('''{JAVA} -jar {GATK} -R {REF} -T HaplotypeCaller \
                -I {input} -o {output.vcf} -L {riskcore_vf} -alleles {riskcore_vf} \
                -gt_mode GENOTYPE_GIVEN_ALLELES -GQB 0 > {log} 2>&1
                ''')

rule calculateriskscore:
    input:
        vcf = rules.callriskscore.output.vcf
    output:
        "tempfiles/{sample}.riskscore.caculated.done"
    message:
        "Calculating risk score"
    run:
        in_database = sample_in_db(wildcards.sample, serie, METRICSDB, 'riskscore')
        if in_database or not input_dict[wildcards.sample]['riskscore']:
            shell('touch {output}') 
        elif not in_database and input_dict[wildcards.sample]['riskscore']:
            genesis = input_dict[wildcards.sample]['genesis']
            target = input_dict[wildcards.sample]['pakket']
            R = RiskScore(input.vcf, genesis)
            score = R.get_score()
            genotypes = json.dumps(R.get_genotypes())
            riskscore_and_genotypes_2db(score, genotypes, serie, wildcards.sample, target, METRICSDB)
            with open(output[0], 'w') as f:
                f.write(f'{score}\n{genotypes}\n')

            

include:
    'callmosaic.rules'

include:
    'samplereport.rules'

include:
    'seriereport.rules'

include:
    'qcplots.rules'

      

    

