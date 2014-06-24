#!/usr/bin/env python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


"""Contains general utility code for the metapaths project"""

try:
    import os, re

    from shutil import rmtree
    from optparse import make_option
    from os import path

    from libs.python_modules.utils.metapathways_utils import *
    from libs.python_modules.utils.sysutil import pathDelim
    from libs.python_modules.utils.utils import *
except:
    print "Cannot load some modules"
    sys.exit(0)
   
PATHDELIM = pathDelim()
class SampleData():
    """Contains the sample related data """

    runlogger = None
    stepslogger = None
    errorlogger = None
    runstatslogger = None


    input_fp = None
    output_dir = None
    sample_name = None

    preprocessed_dir = None
    orf_prediction_dir = None
    genbank_dir = None
    output_run_statistics_dir = None
    blast_results_dir = None
    output_mltreemap_calculations_dir = None
    output_results = None
    output_results_annotation_table_dir  = None
    output_results_megan_dir  = None
    output_results_sequin_dir  = None
    output_results_rpkm_dir  = None
    output_results_mltreemap_dir  = None
    mltreemap_image_output = None
    output_fasta_pf_dir= None
    output_results_pgdb_dir  = None
    output_results_rRNA_dir  = None
    output_results_tRNA_dir  = None


    def __init__(self):
        pass

    def setInputOutput(self, inputFile = None, sample_output_dir = None):
        if inputFile == None and sample_output_dir == None:
            return False

        self.input_fp = inputFile
        self.sample_name = re.sub(r'[.][a-zA-Z]*$','',self.input_fp)
        self.sample_name = path.basename(self.sample_name)
        self.sample_name = re.sub('[.]','_',self.sample_name)

        self.output_dir = sample_output_dir

        self.preprocessed_dir = self.output_dir + PATHDELIM + "preprocessed" + PATHDELIM
        self.orf_prediction_dir =  self.output_dir + PATHDELIM + "orf_prediction"  + PATHDELIM
        self.genbank_dir =  self.output_dir + PATHDELIM + "genbank"  + PATHDELIM
        self.output_run_statistics_dir = self.output_dir + PATHDELIM + "run_statistics"  +PATHDELIM
        self.blast_results_dir =  self.output_dir +  PATHDELIM + "blast_results"  + PATHDELIM
        self.output_mltreemap_calculations_dir = self.output_dir +  PATHDELIM + "mltreemap_calculations"  + PATHDELIM
        self.output_results = self.output_dir + PATHDELIM + "results" + PATHDELIM 
        self.output_results_annotation_table_dir  = self.output_results +  PATHDELIM + "annotation_table"  + PATHDELIM
        self.output_results_megan_dir  = self.output_results + PATHDELIM + "megan"  + PATHDELIM
        self.output_results_sequin_dir  = self.output_results + PATHDELIM + "sequin"  + PATHDELIM
        self.output_results_rpkm_dir  = self.output_results + PATHDELIM + "rpkm"  + PATHDELIM
        self.output_results_mltreemap_dir  = self.output_results +  PATHDELIM + "mltreemap"  + PATHDELIM
        self.mltreemap_image_output = self.output_results_mltreemap_dir  + PATHDELIM + "tables_and_figures" + PATHDELIM 
        self.output_fasta_pf_dir=  self.output_dir + PATHDELIM + "ptools" + PATHDELIM
        self.output_results_pgdb_dir  = self.output_results + PATHDELIM + "pgdb"  + PATHDELIM
        self.output_results_rRNA_dir  = self.output_results +  PATHDELIM + "rRNA"  + PATHDELIM
        self.output_results_tRNA_dir  = self.output_results +  PATHDELIM + "tRNA"   + PATHDELIM
        self.run_stats_file = self.output_run_statistics_dir + PATHDELIM + self.sample_name + ".run.stats.txt"


    def setParameter(self,  parameter, value):
        
        setattr(self, parameter, value)


    def prepareToRun(self):
        self._createFolders()
        self._createLogFiles()


    def _createLogFiles(self):
        self.runlogger = WorkflowLogger(generate_log_fp(self.output_dir, basefile_name='metapathways_run_log'), open_mode='a')
        self.stepslogger = WorkflowLogger(generate_log_fp(self.output_dir, basefile_name='metapathways_steps_log'),open_mode='a')
        self.errorlogger = WorkflowLogger(generate_log_fp(self.output_dir, basefile_name='errors_warnings_log'),open_mode='a')

        print self.output_run_statistics_dir

        self.runstatslogger = WorkflowLogger(generate_log_fp(self.output_run_statistics_dir,  basefile_name = self.sample_name + '.run.stats'),open_mode='a')

    def _createFolders(self):
        checkOrCreateFolder(self.preprocessed_dir)
        checkOrCreateFolder(self.orf_prediction_dir)
        checkOrCreateFolder(self.genbank_dir)
        checkOrCreateFolder(self.output_run_statistics_dir)
        checkOrCreateFolder(self.blast_results_dir)
        checkOrCreateFolder(self.output_mltreemap_calculations_dir)
        checkOrCreateFolder(self.output_results)
        checkOrCreateFolder(self.output_results_annotation_table_dir)
        checkOrCreateFolder(self.output_results_megan_dir)
        checkOrCreateFolder(self.output_results_sequin_dir)
        checkOrCreateFolder(self.output_results_rpkm_dir)
        checkOrCreateFolder(self.output_results_mltreemap_dir)
        #checkOrCreateFolder(mltreemap_image_output) 
        checkOrCreateFolder(self.output_fasta_pf_dir)
        checkOrCreateFolder(self.output_results_pgdb_dir)
        checkOrCreateFolder(self.output_results_rRNA_dir)
        checkOrCreateFolder(self.output_results_tRNA_dir)


