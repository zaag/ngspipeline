import sqlite3

def connect_and_commit_sql(db, sql):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute(sql)
    conn.commit()
    conn.close()


def create_table_for_alignmentmetrics(db):
    sql = """CREATE TABLE IF NOT EXISTS alignment (
    SAMPLE TEXT NOT NULL,
    SERIE TEXT NOT NULL,
    CATEGORY TEXT NOT NULL,
    TOTAL_READS INTEGER NOT NULL,
    PF_READS INTEGER NOT NULL,
    PCT_PF_READS INTEGER NOT NULL,
    PF_NOISE_READS INTEGER NOT NULL,
    PF_READS_ALIGNED INTEGER NOT NULL,
    PCT_PF_READS_ALIGNED REAL NOT NULL,
    PF_ALIGNED_BASES INTEGER NOT NULL,
    PF_HQ_ALIGNED_READS INTEGER NOT NULL,
    PF_HQ_ALIGNED_BASES INTEGER NOT NULL,
    PF_HQ_ALIGNED_Q20_BASES INTEGER NOT NULL,
    PF_HQ_MEDIAN_MISMATCHES INTEGER NOT NULL,
    PF_MISMATCH_RATE REAL NOT NULL,
    PF_HQ_ERROR_RATE REAL NOT NULL,
    PF_INDEL_RATE REAL NOT NULL,
    MEAN_READ_LENGTH INTEGER NOT NULL,
    READS_ALIGNED_IN_PAIRS INTEGER NOT NULL,
    PCT_READS_ALIGNED_IN_PAIRS REAL NOT NULL,
    BAD_CYCLES INTEGER NOT NULL,
    STRAND_BALANCE REAL NOT NULL,
    PCT_CHIMERAS REAL NOT NULL,
    PCT_ADAPTER INTEGER NOT NULL,
    LIBRARY TEXT,
    READ_GROUP TEXT,
    PRIMARY KEY(SAMPLE,SERIE,CATEGORY))
    """
    connect_and_commit_sql(db, sql)

def create_table_for_hsmetrics(db):
    sql = """CREATE TABLE hsmetrics (
    SAMPLE TEXT NOT NULL,
    SERIE TEXT NOT NULL,
    BAIT_SET TEXT NOT NULL,
    GENOME_SIZE INTEGER NOT NULL,
    BAIT_TERRITORY INTEGER NOT NULL,
    TARGET_TERRITORY INTEGER NOT NULL,
    BAIT_DESIGN_EFFICIENCY REAL NOT NULL,
    TOTAL_READS INTEGER NOT NULL,
    PF_READS INTEGER NOT NULL,
    PF_UNIQUE_READS INTEGER NOT NULL,
    PCT_PF_READS INTEGER NOT NULL,
    PCT_PF_UQ_READS INTEGER NOT NULL,
    PF_UQ_READS_ALIGNED INTEGER NOT NULL,
    PCT_PF_UQ_READS_ALIGNED REAL NOT NULL,
    PF_BASES_ALIGNED INTEGER NOT NULL,
    PF_UQ_BASES_ALIGNED INTEGER NOT NULL,
    ON_BAIT_BASES INTEGER NOT NULL,
    NEAR_BAIT_BASES INTEGER NOT NULL,
    OFF_BAIT_BASES INTEGER NOT NULL,
    ON_TARGET_BASES INTEGER NOT NULL,
    PCT_SELECTED_BASES REAL NOT NULL,
    PCT_OFF_BAIT REAL NOT NULL,
    ON_BAIT_VS_SELECTED REAL NOT NULL,
    MEAN_BAIT_COVERAGE REAL NOT NULL,
    MEAN_TARGET_COVERAGE REAL NOT NULL,
    MEDIAN_TARGET_COVERAGE INTEGER NOT NULL,
    PCT_USABLE_BASES_ON_BAIT REAL NOT NULL,
    PCT_USABLE_BASES_ON_TARGET REAL NOT NULL,
    FOLD_ENRICHMENT REAL NOT NULL,
    ZERO_CVG_TARGETS_PCT REAL NOT NULL,
    PCT_EXC_DUPE INTEGER NOT NULL,
    PCT_EXC_MAPQ REAL NOT NULL,
    PCT_EXC_BASEQ INTEGER NOT NULL,
    PCT_EXC_OVERLAP INTEGER NOT NULL,
    PCT_EXC_OFF_TARGET REAL NOT NULL,
    FOLD_80_BASE_PENALTY TEXT NOT NULL,
    PCT_TARGET_BASES_1X REAL NOT NULL,
    PCT_TARGET_BASES_2X INTEGER NOT NULL,
    PCT_TARGET_BASES_10X INTEGER NOT NULL,
    PCT_TARGET_BASES_20X INTEGER NOT NULL,
    PCT_TARGET_BASES_30X INTEGER NOT NULL,
    PCT_TARGET_BASES_40X INTEGER NOT NULL,
    PCT_TARGET_BASES_50X INTEGER NOT NULL,
    PCT_TARGET_BASES_100X INTEGER NOT NULL,
    HS_LIBRARY_SIZE REAL,
    HS_PENALTY_10X INTEGER NOT NULL,
    HS_PENALTY_20X INTEGER NOT NULL,
    HS_PENALTY_30X INTEGER NOT NULL,
    HS_PENALTY_40X INTEGER NOT NULL,
    HS_PENALTY_50X INTEGER NOT NULL,
    HS_PENALTY_100X INTEGER NOT NULL,
    AT_DROPOUT INTEGER NOT NULL,
    GC_DROPOUT REAL NOT NULL,
    HET_SNP_SENSITIVITY REAL NOT NULL,
    HET_SNP_Q INTEGER NOT NULL,
    LIBRARY TEXT,
    READ_GROUP TEXT,
    PRIMARY KEY(SAMPLE,SERIE))
    """
    connect_and_commit_sql(db, sql)

def create_table_for_insertsizemetrics(db):
    sql = """CREATE TABLE insertsize (
    SAMPLE TEXT NOT NULL,
    SERIE TEXT NOT NULL,
    MEDIAN_INSERT_SIZE INTEGER NOT NULL,
    MEDIAN_ABSOLUTE_DEVIATION INTEGER NOT NULL,
    MIN_INSERT_SIZE INTEGER NOT NULL,
    MAX_INSERT_SIZE INTEGER NOT NULL,
    MEAN_INSERT_SIZE REAL NOT NULL,
    STANDARD_DEVIATION REAL NOT NULL,
    READ_PAIRS INTEGER NOT NULL,
    PAIR_ORIENTATION TEXT NOT NULL,
    WIDTH_OF_10_PERCENT INTEGER NOT NULL,
    WIDTH_OF_20_PERCENT INTEGER NOT NULL,
    WIDTH_OF_30_PERCENT INTEGER NOT NULL,
    WIDTH_OF_40_PERCENT INTEGER NOT NULL,
    WIDTH_OF_50_PERCENT INTEGER NOT NULL,
    WIDTH_OF_60_PERCENT INTEGER NOT NULL,
    WIDTH_OF_70_PERCENT INTEGER NOT NULL,
    WIDTH_OF_80_PERCENT INTEGER NOT NULL,
    WIDTH_OF_90_PERCENT INTEGER NOT NULL,
    WIDTH_OF_99_PERCENT INTEGER NOT NULL,
    LIBRARY TEXT,
    READ_GROUP TEXT,
    PRIMARY KEY(SAMPLE,SERIE))
    """
    connect_and_commit_sql(db, sql)

def create_table_for_callable_summary(db):
    sql = """CREATE TABLE IF NOT EXISTS callable
    (REF_N INTEGER NOT NULL,
    CALLABLE INTEGER NOT NULL,
    NO_COVERAGE INTEGER NOT NULL,
    LOW_COVERAGE INTEGER NOT NULL,
    EXCESSIVE_COVERAGE INTEGER NOT NULL,
    POOR_MAPPING_QUALITY INTEGER NOT NULL,
    TARGET TEXT NOT NULL,
    SAMPLE TEXT NOT NULL,
    SERIE TEXT NOT NULL,
    PRIMARY KEY(SAMPLE, TARGET, SERIE))
    """
    connect_and_commit_sql(db, sql)

def create_table_for_perc_covered(db):
    sql = """CREATE TABLE IF NOT EXISTS procenttargetcovered
    (SAMPLE TEXT NOT NULL,
    TARGET TEXT NOT NULL,
    SERIE TEXT NOT NULL,
    PERCENTAGE REAL NOT NULL,
    PRIMARY KEY(SAMPLE, TARGET, SERIE))
    """
    connect_and_commit_sql(db, sql)



class TargetDBCreator:

    def __init__(self, db):
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()

    def create_genesis_table(self):
        genesis = """CREATE TABLE genesis (
        genesis TEXT NOT NULL,
        capture TEXT NOT NULL,
        pakket TEXT NOT NULL,
        panel TEXT,
        cnvscreening BOOLEAN NOT NULL DEFAULT 0,
        cnvdiagnostiek BOOLEAN NOT NULL DEFAULT 0,
        mozaiekdiagnostiek BOOLEAN NOT NULL DEFAULT 0,
        PRIMARY KEY (genesis))
        """
        self.c.execute(genesis)
        self.conn.commit()

    def create_aandoeningen_table(self):
        aandoeningen = """CREATE TABLE aandoeningen (
        genesis TEXT NOT NULL,
        aandoening TEXT NOT NULL,
        PRIMARY KEY (genesis))
        """
        self.c.execute(aandoeningen)
        self.conn.commit()

    def create_captures_table(self):
        captures = """CREATE TABLE captures (
        capture TEXT NOT NULL,
        versie INTEGER NOT NULL,
        oid INTEGER NOT NULL,
        verdund BOOLEAN NOT NULL,
        grootte INTEGER NOT NULL,
        genen TEXT,
        PRIMARY KEY (capture, versie))
        """
        self.c.execute(captures)
        self.conn.commit()

    def create_pakketten_table(self):
        pakketten = """CREATE TABLE pakketten (
        pakket TEXT NOT NULL,
        versie INTEGER NOT NULL,
        grootte INTEGER,
        genen TEXT,
        PRIMARY KEY (pakket, versie))
        """
        self.c.execute(pakketten)
        self.conn.commit()

    def create_panels_table(self):
        panels = """CREATE TABLE panels (
        panel TEXT NOT NULL,
        versie INTEGER NOT NULL,
        grootte INTEGER,
        genen TEXT,
        PRIMARY KEY (panel, versie))
        """
        self.c.execute(panels)
        self.conn.commit()

    def create_all(self):
        try:
            self.create_genesis_table()
        except sqlite3.OperationalError as e:
            print(e)
        try:
            self.create_aandoeningen_table()
        except sqlite3.OperationalError as e:
            print(e)
        try:
            self.create_captures_table()
        except sqlite3.OperationalError as e:
            print(e)
        try:
            self.create_pakketten_table()
        except sqlite3.OperationalError as e:
            print(e)
        try:
            self.create_panels_table()
        except sqlite3.OperationalError as e:
            print(e)


class MetricsDBReader:

    def __init__(self, db, sampleID, serie, target):
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()
        self.sampleID = sampleID
        self.serie = serie
        self.target = target

    def get_alignmetrics(self):
        read1 = '''SELECT DISTINCT
        TOTAL_READS,
        round(PCT_CHIMERAS*100),
        round(MEAN_READ_LENGTH),
        round(PCT_READS_ALIGNED_IN_PAIRS*100)
        FROM {table}
        WHERE (CATEGORY = 'FIRST_OF_PAIR'
        AND sample = '{dnr}'
        AND serie = '{serie}')
        '''.format(table='alignment', dnr=self.sampleID, serie=self.serie)

        read2 = '''SELECT DISTINCT
        round(MEAN_READ_LENGTH),
        round(PCT_READS_ALIGNED_IN_PAIRS*100)
        FROM {table}
        WHERE (CATEGORY = 'SECOND_OF_PAIR'
        AND sample = '{dnr}'
        AND serie = '{serie}')
        '''.format(table='alignment', dnr=self.sampleID, serie=self.serie)

        self.c.execute(read1)
        read1out = [val for tup in self.c.fetchall() for val in tup]

        self.c.execute(read2)
        read2out = [val for tup in self.c.fetchall() for val in tup]
        return read1out + read2out

    def get_hsmetrics(self):
        sql = '''SELECT DISTINCT
        TARGET_TERRITORY,
        BAIT_TERRITORY,
        round(PCT_PF_UQ_READS*100),
        round(PCT_SELECTED_BASES*100)
        FROM {table}
        WHERE (sample = '{dnr}'
        AND serie = '{serie}')
        '''.format(table='hsmetrics', dnr=self.sampleID, serie=self.serie)
        self.c.execute(sql)
        out = [val for tup in self.c.fetchall() for val in tup]
        return out

    def get_callable(self):
        total_sql = """SELECT REF_N, CALLABLE,
        NO_COVERAGE, LOW_COVERAGE, EXCESSIVE_COVERAGE
        FROM callable
        WHERE (SAMPLE='{s}' AND TARGET='{t}')
        """.format(s=self.sampleID, t=self.target)

        callable_sql = """SELECT CALLABLE
        FROM callable
        WHERE (SAMPLE='{s}' AND TARGET='{t}')
        """.format(s=self.sampleID, t=self.target)

        self.c.execute(total_sql)
        all_out = [val for tup in self.c.fetchall() for val in tup]

        self.c.execute(callable_sql)
        callable_out = [val for tup in self.c.fetchall() for val in tup]

        return (sum(callable_out) / sum(all_out)) * 100

    def get_perc_target_covered(self):
        percentage = """SELECT PERCENTAGE
        FROM procenttargetcovered
        WHERE (SAMPLE='{s}' AND TARGET='{t}')
        """.format(s=self.sampleID, t=self.target)
        self.c.execute(percentage)
        return self.c.fetchone()[0]

    def get_all(self):
        print('Alignment: {}'.format(self.get_alignmetrics()),
              'HSmetrics: {}'.format(self.get_hsmetrics()),
              'Callable: {}'.format(self.get_callable()),
              'Perctargetcovered: {}'.format(self.get_perc_target_covered())
              )
