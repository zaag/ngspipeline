import os
import csv
import pandas as pd
import subprocess
import sqlite3
import vcf
import sys

sys.path.append('/home/manager/Documents/ngsscriptlibrary')
sys.path.append('/data/dnadiag/ngsscriptlibrary')

from ngsscriptlibrary import TargetDatabase

def get_rs_gpos_dict():
    rs_gpos_dict = dict()
    basedir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(basedir, 'docs', 'snpcheck.csv')) as f:
        reader = csv.reader(f)
        for line in reader:
            rs, locus, genotype = line
            rs_gpos_dict[rs] = locus
            rs_gpos_dict[locus] = rs
    return rs_gpos_dict


def get_header(samplesheet):
    "Read samplesheet and find line with Sample_ID. Return integer."
    with open(samplesheet, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('Sample_ID'):
                return i


def parse_samplesheet(samplesheet, db):
    """Read samplesheet part with sample info. Use ngsscriptlibrary to get
    targets to use and analyses to perform. Return dict.
    """
    samples = list()
    with open(samplesheet) as f:
        reader = csv.reader(f)
        for i, line in enumerate(reader):
            if i > get_header(samplesheet):
                samples.append(line)
    TD = TargetDatabase(db)
    samples_todo = dict()
    for line in samples:
        if 'R' in line[0]:
            continue
        genesis = line[-1]
        genesis = genesis.replace('.NGS', '')
        samples_todo[line[0]] = TD.get_todo(genesis)
        samples_todo[line[0]]['serie'] = line[-2]
        vcapture = samples_todo[line[0]]['capture']
        oid = TD.get_oid_for_vcapture(vcapture)
        samples_todo[line[0]]['oid'] = oid
    return samples_todo


def get_file_locations(todo, targetrepo):
    """Read dict with targets and analyses and add correct file locations.
    Return dict
    """
    for s in todo.keys():
        picard = '{}_target.interval_list'.format(todo[s]['capture'])
        picard = os.path.join(targetrepo, 'captures', picard)

        cnvtarget = '{}_target.bed'.format(todo[s]['capture'])
        cnvtarget = os.path.join(targetrepo, 'captures', cnvtarget)

        annot = '{}_target.annotated'.format(todo[s]['capture'])
        annot = os.path.join(targetrepo, 'captures', annot)

        cap_is_pakket = todo[s]['capispakket']

        if cap_is_pakket:
            varcal = '{}_generegions.bed'.format(todo[s]['capture'])
            varcal = os.path.join(targetrepo, 'captures', varcal)
            sanger = '{}_target.bed'.format(todo[s]['capture'])
            sanger = os.path.join(targetrepo, 'captures', sanger)
            pakket = '{}_target.bed'.format(todo[s]['capture'])
            pakket = os.path.join(targetrepo, 'captures', pakket)
        elif not cap_is_pakket:
            varcal = '{}_generegions.bed'.format(todo[s]['pakket'])
            varcal = os.path.join(targetrepo, 'pakketten', varcal)
            sanger = '{}_target.bed'.format(todo[s]['pakket'])
            sanger = os.path.join(targetrepo, 'pakketten', sanger)
            pakket = '{}_target.bed'.format(todo[s]['pakket'])
            pakket = os.path.join(targetrepo, 'pakketten', pakket)

        if todo[s]['panel'] is not None:
            sanger = '{}_target.bed'.format(todo[s]['panel'])
            sanger = os.path.join(targetrepo, 'panels', sanger)

        todo[s]['annot'] = annot
        todo[s]['picard'] = picard
        todo[s]['cnvtarget'] = cnvtarget
        todo[s]['pakkettarget'] = pakket
        todo[s]['varcal'] = varcal
        todo[s]['sanger'] = sanger

    return todo


def create_recalinput_file(samples):
    """Create a file with a list of BAM-files (1 per sample) to use as input
    for the recalibration step in the pipeline.
    """
    with open('recalinput.list', 'w') as f:
        for sample in samples:
            f.write('output/{}.DM.bam\n'.format(sample))


def sample_in_db(sample, serie, db, table):
    "Check if sample-serie combi is in database. Return boolean"
    sql = """SELECT * FROM {}
    WHERE(SAMPLE='{}' AND SERIE={})
    """.format(table, sample, serie)

    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute(sql)

    if c.fetchone():
        conn.close()
        return True
    else:
        conn.close()
        return False


def get_pakketten(todo):
    pakketten = list()
    for s in todo.keys():
        pakket = todo[s]['pakket']
        if pakket not in pakketten:
            pakketten.append(pakket)
    return pakketten


def get_captures(todo):
    captures = list()
    for s in todo.keys():
        capture = todo[s]['capture']
        if capture not in captures:
            captures.append(capture)
    return captures


def parse_sangers(sangerfile):
    """Read file into dataframe. Create a dict with loci as keys and WT, HET
    or HOM as value. Return a dict.
    """
    sanger = dict()
    sangerout = pd.read_csv(sangerfile, sep='\t')
    sangerout.set_index('mut Effect', inplace=True)

    for locus in sangerout.index:
        if str(locus) == 'nan':
            next
        elif 'het' in sangerout.loc[locus]['Nuc Change']:
            sanger[locus] = 'HET'
        elif 'homo' in sangerout.loc[locus]['Nuc Change']:
            sanger[locus] = 'HOM'
        elif ('het' not in sangerout.loc[locus]['Nuc Change'] and
              'homo' not in sangerout.loc[locus]['Nuc Change']):
            sanger[locus] = 'WT'
    return sanger


def parse_taqman(taqmanfile):
    """Parse CSV file with rsID, genotype info. Create a dict with rsID as keys
    and WT, HET, HOM or NoTaqman (no call) as value. Return a dict.
    """
    rsid_to_locus = get_rs_gpos_dict()
    taqman = dict()
    with open (taqmanfile, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            rs, call = line
            if rs.endswith(')'):
                rs = rs.split(' ')[0]
            else:
                rs = rs.rsplit()[0]
            if call != 'WT' and call != 'HET' and call != 'HOM':
                call = 'NoTaqMan'
            locus = rsid_to_locus[rs]
            taqman[locus] = call
    return taqman


def parse_ngssnpcheck(ngssnpcheck):
    """Parse vcf with vcfreader. Create a dict with loci as key and WT, HET,
    HOM or NC (no coverage) as value. Return a dict.
    """
    d = dict()
    d['0/0'] = 'WT'
    d['0/1'] = 'HET'
    d['1/1'] = 'HOM'
    d['./.'] = 'NC'

    ngsdict = dict()
    vcfreader = vcf.Reader(open(ngssnpcheck, 'r'))

    for record in vcfreader:
        for call in record:
            if d[call['GT']] == 'NC' and record.INFO['DP'] > 29:
                ngsdict['{}:{}'.format(record.CHROM, record.POS)] = 'WT'
            else:
                ngsdict['{}:{}'.format(record.CHROM,
                                       record.POS)] = d[call['GT']]
    return ngsdict


def compare_snpchecks(sangerdict, ngsdict):
    """Compare values from 2 dicts with overlapping keys. NGS-dict contains
    all keys from sangerdict. Create dict with loci as keys and ok or ERROR as
    values. Ok when both values are the same, ERROR if not. Return a dict.
    """
    out = dict()
    for k, v in sangerdict.items():
        try:
            if ngsdict[k] == v:
                out[k] = 'ok'
            elif ngsdict[k] != v:
                out[k] = 'ERROR'
        except KeyError as e:
            out[k] = 'NoNGS'
        if v == 'NoTaqMan':
            out[k] = 'NoTaqMan'
    return out


def print_extra_ngscalls(compfile, extrangsfile):
    with open(extrangsfile) as fin, open(compfile, 'a') as fout:
        for line in fin:
            line.strip()
            fout.write(line)


def read_summary(summary):
    "Parse GATK CallableLoci summary. Return dict."
    d = dict()
    with open(summary, 'r') as f:
        next(f)
        for line in f:
            x, y = line.lstrip().split()
            d[x] = y
    return d


def summary2db(sampleID, data, target, serie, db):
    "Insert dict with CallableLoci summary into database."
    data['TARGET'] = target
    data['SAMPLE'] = sampleID
    data['SERIE'] = serie
    conn = sqlite3.connect(db)
    c = conn.cursor()
    sql = """INSERT INTO callable ({headers})
    VALUES ('{values}');
    """.format(headers=', '.join(data.keys()),
               values="""', '""".join(data.values()))
    try:
        c.execute(sql)
    except sqlite3.IntegrityError as e:
        print(e)
    else:
        conn.commit()
    conn.close()


def annotate_callables(bed, annottarget, tempbed):
    "Intersect annotated target with callable region. Return dataframe."
    bash = 'bedtools intersect -wa  -a {} -wb -b {} > {}'.format(
        annottarget, bed, tempbed)
    subprocess.call(bash, shell=True)
    df = pd.DataFrame.from_csv(tempbed, sep='\t', header=None, index_col=None)
    df = df.drop([4], axis=1)
    df.columns = ['chromosome', 'targetstart', 'targetend',
                  'gene', 'regionstart', 'regionend', 'callable']

    integerlist = ['targetstart', 'targetend',
                   'regionstart', 'regionend']

    df[integerlist] = df[integerlist].apply(pd.to_numeric)
    df['regionsize'] = df['regionend'] - df['regionstart']
    df['targetsize'] = df['targetend'] - df['targetstart']
    df.set_index(['chromosome', 'targetstart', 'targetend'], inplace=True)
    return df


def get_basecounts(fastq):
    """Use subprocess/bash to get all bases from a fastq in a string.
    Count A, C, G, T and divide by total length of string. Return tuple.
    """
    if fastq.endswith('.gz'):
        bash = "zcat {} | paste - - - - | cut -f 2 | tr -d '\n'".format(fastq)
    else:
        bash = "cat {} | paste - - - - | cut -f 2 | tr -d '\n'".format(fastq)

    proc = subprocess.Popen(bash,
                            stdout=subprocess.PIPE,
                            shell=True)

    out, err = proc.communicate()
    out = str(out)

    countA = float(out.count('A') / len(out))
    countT = float(out.count('T') / len(out))
    countC = float(out.count('C') / len(out))
    countG = float(out.count('G') / len(out))

    return(countA, countT, countC, countG)


def get_base_count_filelist(outputfile, filelist):
    "Feed all fastqs in filelist and combine output."

    if len(filelist) == 0:
        return 'No fastq R1 files found'

    else:
        counts = [(pd.DataFrame({f: get_basecounts(f)},
                                index=['A', 'T', 'C', 'G']
                                )
                   ) for f in filelist if 'R2' not in f]

        df = pd.concat(counts, axis=1).transpose()
        df.index = [path.split('/')[-1] for path in df.index]
        df.index = [fn.split('_')[0] for fn in df.index]
        df.to_csv(outputfile, sep='\t')

def get_loci(target):
    "Read BED file and create chr:pos for each targetbase. Return list."
    loci = list()
    with open(target) as f:
        for line in f:
            chromosome, start, end, *_ = line.split()
            start = int(start)
            end = int(end)
            while start <= end:
                locus = '{}:{}'.format(chromosome, start)
                loci.append(locus)
                start += 1
    return loci


def calc_perc_target_covered(docfile, filter_regions=False, target=None):
    """Calculate targetpercentage covered with >30 reads from GATK's
    DepthOfCoverage output. Filter regions if requested. Return float.
    """
    df = pd.read_csv(docfile, usecols=['Total_Depth', 'Locus'],
                     sep='\t', index_col='Locus')
    if filter_regions:
        loci = get_loci(target)
        df = df.loc[loci]
    return float((1 - df[df['Total_Depth'] < 30].count() / df.count()) * 100)


def perc_target_covered2db(sampleID, percentage, pakket, serie, db):
    "Create database connection and write perc. target covered >30x to table."
    conn = sqlite3.connect(db)
    c = conn.cursor()
    sql = '''INSERT INTO procenttargetcovered
    (SAMPLE, TARGET, SERIE, PERCENTAGE)
    VALUES ('{}', '{}', '{}', '{}')
    '''.format(sampleID, pakket, serie, percentage)
    try:
        c.execute(sql)
    except sqlite3.IntegrityError as e:
        pass
    else:
        conn.commit()
    conn.close()


def get_line(picardfile, phrase):
    """Read file. Report line number that starts with phrase and the first
    blank line after that line. Return tuple of integers
    """
    with open(picardfile, 'r') as f:
        start = 1000
        end = 1000
        for i, line in enumerate(f):
            if line.startswith(phrase):
                start = i
            if line.startswith('\n'):
                end = i
            if start < end:
                return (start, end)


def readmetrics(sampleID, serie, metrics):
    "Read block (line x through y) from file with pandas. Return dataframe."
    start, end = get_line(metrics, '## METRICS CLASS')
    df = pd.read_csv(metrics, sep='\t', header=start, nrows=end - start - 2)
    df['SAMPLE'] = sampleID
    df['SERIE'] = serie
    df.set_index(['SAMPLE', 'SERIE'], inplace=True)
    return df


def metrics2db(df, db, table):
    "Create database connection and write picard results to table."
    conn = sqlite3.connect(db)
    try:
        df.to_sql(table, conn, if_exists='append')
    except sqlite3.OperationalError as e:
        print(e)
    conn.close()
