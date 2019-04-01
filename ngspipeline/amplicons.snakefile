import os
import csv
import glob
import pysam
import xlsxwriter

from decimal import Decimal
from collections import namedtuple
from collections import defaultdict

from ngsscriptlibrary import mosaic
from ngsscriptlibrary import samplesheetinfo2db
from ngsscriptlibrary.parsing import parse_samplesheet_for_pipeline

# Define variables from config file
GATK = config["GATK"]
REF = config["REF"]
TARGETREPO = config["TARGET"]
DBREPO = config['DB']
JAVA = config['JAVA']

# Define files
METRICSDB = os.path.join(DBREPO, 'amplicons.sqlite')
SAMPLEDB = os.path.join(DBREPO, 'samplesheets.sqlite')
TARGETDB = os.path.join(TARGETREPO, 'varia', 'captures.sqlite')

SAMPLESHEET = 'SampleSheet.csv'


def get_file_locations(todo, targetrepo):
    """Read dict with todo and add targetfile location for amplicon analysis.
    Return dict.
    """
    for s in todo.keys():
        if not todo[s]['amplicon']:
            continue
        genesis = todo[s]['genesis']

        if genesis.startswith('chr'):
            target = genesis
        else:
            target = os.path.join(targetrepo, 'amplicons', f'{genesis}_target.bed')
        todo[s]['target'] = target
    return todo


def parse_doc(fn, ref, loci):
    data = dict()
    with open(fn) as f:
        f_reader = csv.reader(f, delimiter='\t')
        _header = next(f_reader)
        for line in f_reader:
            nonref = list()
            locus, _TD, _ADS, DP, basecounts = line

            if locus not in loci:
                continue
            refbase = ref[locus]
            for bases in basecounts.split(' '):
                base, cov = bases.split(':')
                cov = int(cov)
                try:
                    basep = int(cov) / int(DP)
                except ZeroDivisionError:
                    basep = 0
                basep = Decimal(basep)
                if refbase == base:
                    refp = basep
                elif refbase != base:
                    nonref.append((base, basep))

            locus_data = dict()
            locus_data['DP'] = int(DP)
            locus_data['refP'] = refp
            locus_data['refbase'] = refbase
            for base, basep in nonref:
                locus_data[base] = basep
            data[locus] = locus_data
    return data


def get_loci_from_docfile(docfile):
    loci = list()
    with open(docfile) as f:
        header = next(f)
        for line in f:
            locus, *_ = line.split()
            loci.append(locus)
    return loci


def get_ref_dict(loci, ref):
    refd = dict()
    for locus in loci:
        base = pysam.faidx(ref,
                           '{}:{}-{}'.format(locus.split(':')[0],
                                             locus.split(':')[1],
                                             locus.split(':')[1])
                           ).split('\n')[1]
        refd[locus] = base
    return refd


def create_output(outfile, data, loci):

    with open(outfile, 'w') as f_out:
        f_out.write('Sample')
        for locus in loci:
            f_out.write(f'\t{locus} DP\tA\tC\tT\tG\tD\tI\trefP')
        f_out.write('\n')

        for sampleID, sample_data in data.items():
            f_out.write(f'{sampleID}\t')
            insertions = sample_data['insertions']
            for locus in loci:
                dp = sample_data[locus]['DP']
                refp = sample_data[locus]['refP']
                refbase = sample_data[locus]['refbase']
                f_out.write(f'{dp}\t')

                for base in 'A C T G D'.split():
                    if base == refbase:
                        f_out.write('.\t')
                    else:
                        perc = sample_data[locus][base]
                        f_out.write(f'{perc:.4}\t')
                if locus in insertions:
                    insertion = Decimal(insertions[locus] / dp)
                else:
                    insertion = Decimal(0)
                f_out.write(f'{insertion:.4}\t')
                refp = refp - insertion
                f_out.write(f'{refp:.4}\t')
            f_out.write('\n')

def create_output_excel(outfile, data, loci):
    wb = xlsxwriter.Workbook(outfile)
    wb.set_properties({
        'title':    'Amplicons',
        'subject':  'MiSEQUENCING',
        'author':   'Scipio Africanus',
        'comments': 'Created with Python and XlsxWriter'})

    lowcov = wb.add_format()
    lowcov.set_bg_color('orange')
    positive = wb.add_format()
    positive.set_bg_color('green')
    headerformat = wb.add_format()
    headerformat.set_font_size(16)
    underlined = wb.add_format()
    underlined.set_bottom()    

    ws = wb.add_worksheet('Percentages')
    ws_analyse = wb.add_worksheet('Analyse info')

    ws.set_column('A:B', 18)
    ws_analyse.set_column('A:B', 20)

    analyse_info = [('Referentie', 'lifescope.hg19.fa'), 
                    ('Aligner', 'bwa-mem 0.7.12-r1039'), 
                    ('GATK', 'GATK 3.8')]

    row = 0 
    col = 0

    ws_analyse.write(row, col, 'ANALYSE INFO', headerformat)
    row += 2

    for name, version in analyse_info:
        ws_analyse.write(row, col, name)
        col += 1
        ws_analyse.write(row, col, version)
        row += 1
        col = 0

    row = 0 
    col = 0


    ws.write(row, col, 'PERCENTAGE PER LOCUS', headerformat)
    row += 2

    ws.write(row, col, 'Sample', underlined)
    col += 1

    for locus in loci:
        for _ in [f'{locus} DP', 'A', 'C', 'T', 'G', 'D', 'I', 'reP']:
            ws.write(row, col, _, underlined)
            col += 1
            
    row += 1
    col = 0

    for sampleID, sample_data in data.items():
        ws.write(row, col, f'{sampleID}')
        col += 1 
        
        insertions = sample_data['insertions']
        
        for locus in loci:
            dp = sample_data[locus]['DP']
            refp = sample_data[locus]['refP']
            refbase = sample_data[locus]['refbase']
            if int(dp) < 500:
                ws.write(row, col, f'{dp}', lowcov)
            else:
                ws.write(row, col, f'{dp}')
            col += 1

            for base in 'A C T G D'.split():
                if base == refbase:
                    ws.write(row, col, '.')
                    col += 1
                else:
                    perc = sample_data[locus][base]
                    if perc >= Decimal(0.01):
                        ws.write(row, col, f'{perc:.4}', positive)
                    else:
                        ws.write(row, col, f'{perc:.4}')
                    col += 1
            if locus in insertions:
                insertion = Decimal(insertions[locus] / dp)
            else:
                insertion = Decimal(0)
            if insertion >= Decimal(0.01):
                ws.write(row, col, f'{insertion:.4}', positive)
            else:
                ws.write(row, col, f'{insertion:.4}')
            col += 1
            refp = refp - insertion
            ws.write(row, col, f'{refp:.4}')
            col += 1
        row += 1
        col = 0
    
    
    wb.close()

input_dict = parse_samplesheet_for_pipeline(SAMPLESHEET, TARGETDB)
input_dict = get_file_locations(input_dict, TARGETREPO)
samples = [s for s in input_dict.keys() if input_dict[s]['amplicon']]
serie = [input_dict[_]['serie'] for _ in input_dict.keys()][0]

rule all:
    input:
        f"output/MS{serie}_amplicon.xlsx"


rule prepare:
    input:
        {SAMPLESHEET}
    output:
        expand("reads/{sample}.R1.fastq.gz", sample=samples),
        expand("reads/{sample}.R2.fastq.gz", sample=samples)
    run:
        for sample in samples:
            samplesheetinfo2db(input_dict[sample], sample, SAMPLEDB)
            R1 = glob.glob(f'{sample}_*R1*.gz')
            R2 = glob.glob(f'{sample}_*R2*.gz')
            if R1 and R2:
                os.rename(R1[0], f'reads/{sample}.R1.fastq.gz')
                os.rename(R2[0], f'reads/{sample}.R2.fastq.gz')


rule mapreads:
    input:
        "reads/{sample}.R1.fastq.gz",
        "reads/{sample}.R2.fastq.gz"
    output:
        bam = "output/{sample}.sorted.bam",
        index = "output/{sample}.sorted.bam.bai",
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
        samtools sort -T {wildcards.sample}.tmp -O bam - > {output.bam} \
        && samtools index {output.bam}) > {log}  2>&1
        '''


rule depthofcoverage:
    input:
        bam = rules.mapreads.output.bam
    output:
        doc = temp("tempfiles/{sample}.DoC"),
        cumcovcounts = temp("tempfiles/{sample}.DoC.sample_cumulative_coverage_counts"),
        cumcovprop = temp("tempfiles/{sample}.DoC.sample_cumulative_coverage_proportions"),
        samplestats = temp("tempfiles/{sample}.DoC.sample_statistics"),
        samplesum = temp("tempfiles/{sample}.DoC.sample_summary")
    message:
        "Calculating DepthOfCoverage with GATK"
    log:
        "logfiles/{sample}.DeptOfCoverage.log"
    run:
        target = input_dict[wildcards.sample]['target']
        shell('''{JAVA} -jar {GATK} -R {REF} -T DepthOfCoverage \
        -I {input.bam} -o {output[0]} \
        -L {target} \
        --minBaseQuality 20 \
        --minMappingQuality 20 \
        --omitIntervalStatistics \
        --printBaseCounts \
        --countType COUNT_FRAGMENTS \
        -dels  > {log} 2>&1
        ''')


rule calculatepercentages:
    input:
        expand(rules.depthofcoverage.output.doc, sample=samples)
    output:
        excel = f"output/MS{serie}_amplicon.xlsx"
    message:
        "Calculating ref and var percentages"
    log:
        "logfiles/DoC2Percentage.log"        
    run:
        data = dict()
        loci = get_loci_from_docfile(input[0])
        refd = get_ref_dict(loci, REF)
        for sample in samples:
            target = input_dict[sample]['target']
            if target.startswith('chr'):
                insertions = mosaic.get_indel_dict_for_locus(f'output/{sample}.sorted.bam', target)
            else:
                insertions = mosaic.get_indel_dicts(f'output/{sample}.sorted.bam', target)
            data[sample]  = parse_doc(f'tempfiles/{sample}.DoC', refd, loci)
            data[sample]['insertions'] = insertions
        create_output_excel(output.excel, data, loci)
