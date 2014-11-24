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
    from os import path, _exit

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


    input_file = None
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
    ncbi_params_file = None
    ncbi_sequin_sbt = None
    bwa_folder = None

    stages = []
    stages_context = {}

    def __init__(self):
        pass

    def getContexts(self):
        contexts = []
        for name in self.stages:
           if name in self.stages_context:
              contexts.append( self.stages_context[name] )
        return contexts

    def getContextBlocks(self):
        contextBlocks  =[]

        for block in self.stages:
          contexts = []
          for name in block:
             if name in self.stages_context:
                 contexts.append( self.stages_context[name] )
          contextBlocks.append(contexts)

        return contextBlocks


    def numJobs(self):
        return len(self.stages)

    def clearJobs(self):
       self.stages = []
       self.stages_context = {}



    def setInputOutput(self, inputFile = None, sample_output_dir = None):
        if inputFile == None and sample_output_dir == None:
            return False

        self.input_file = inputFile
        self.sample_name = re.sub(r'[.][a-zA-Z]*$','',self.input_file)
        self.sample_name = path.basename(self.sample_name)
        self.sample_name = re.sub('[.]','_',self.sample_name)

        self.rpkm_input_dir  = path.dirname(inputFile) + PATHDELIM  + "reads" 

        self.output_dir = sample_output_dir

        self.preprocessed_dir = self.output_dir + PATHDELIM + "preprocessed" + PATHDELIM
        self.orf_prediction_dir =  self.output_dir + PATHDELIM + "orf_prediction"  + PATHDELIM
        self.genbank_dir =  self.output_dir + PATHDELIM + "genbank"  + PATHDELIM
        self.output_run_statistics_dir = self.output_dir + PATHDELIM + "run_statistics"  +PATHDELIM
        self.blast_results_dir =  self.output_dir +  PATHDELIM + "blast_results"  + PATHDELIM
        self.output_mltreemap_calculations_dir = self.output_dir +  PATHDELIM + "mltreemap_calculations"  + PATHDELIM
        self.bwa_folder = self.output_dir + PATHDELIM + "bwa" + PATHDELIM 
        self.output_results = self.output_dir + PATHDELIM + "results" + PATHDELIM 
        self.output_results_annotation_table_dir  = self.output_results +  PATHDELIM + "annotation_table"  + PATHDELIM
        self.output_results_megan_dir  = self.output_results + PATHDELIM + "megan"  + PATHDELIM
        self.output_results_sequin_dir  = self.output_results + PATHDELIM + "sequin"  + PATHDELIM
        self.output_results_rpkm_dir  = self.output_results+ PATHDELIM  + "rpkm" 
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


    def getType(self):
        if not hasattr(self, 'SEQ_TYPE'):
           return None

        if not hasattr(self, 'FILE_TYPE'):
           return None

        if self.FILE_TYPE=='FASTA':
            if self.SEQ_TYPE =='NUCL':
               return self.SEQ_TYPE + '-' + self.FILE_TYPE

            if self.SEQ_TYPE =='AMINO':
               return self.SEQ_TYPE + '-' + self.FILE_TYPE

        if self.FILE_TYPE=='GENBANK':
            if self.SEQ_TYPE =='NOT-USED':
               return 'AMINO' '-' + self.FILE_TYPE + '-' + 'UNANNOT'

        return 'UNKNOWN'


    def _createLogFiles(self):
        self.runlogger = WorkflowLogger(generate_log_fp(self.output_dir, basefile_name='metapathways_run_log'), open_mode='a')
        self.stepslogger = WorkflowLogger(generate_log_fp(self.output_dir, basefile_name='metapathways_steps_log'),open_mode='a')
        self.errorlogger = WorkflowLogger(generate_log_fp(self.output_dir, basefile_name='errors_warnings_log'),open_mode='a')
        self.runstatslogger = WorkflowLogger(generate_log_fp(self.output_run_statistics_dir,  basefile_name = self.sample_name + '.run.stats'),open_mode='a')

    def _createFolders(self):
        checkOrCreateFolder(self.preprocessed_dir)
        checkOrCreateFolder(self.orf_prediction_dir)
        checkOrCreateFolder(self.genbank_dir)
        checkOrCreateFolder(self.output_run_statistics_dir)
        checkOrCreateFolder(self.blast_results_dir)
        checkOrCreateFolder(self.output_mltreemap_calculations_dir)
        checkOrCreateFolder(self.bwa_folder)
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

    def addPipeLineStage(self, stepName, inputs= [], outputs = [], status = 'yes'):

        stagecontext = Context()
        stagecontext.inputs  = inputs
        stagecontext.outputs = outputs
        stagecontext.name= stepName
        stagecontext.status = context

        stages_context[stepName] = stagecontext

    def addContextBlock(self, contextBlock):
          self.contextBlock.append(contextBlock)

    def addContexts(self, contextBlock):
        stages = []

        for contexts in contextBlock:
           for context in contexts:
              stages.append(context.name)
              self.stages_context[context.name] = context

        self.stages.append( stages)


    def hasPToolsInput(self):
        """ checks if the ptools folder has the right inputs"""
        files = [ '0.pf', '0.fasta', 'genetic-elements.dat', 'organism-params.dat']

        for file in files:
            if not doesFileExist( self.output_fasta_pf_dir +PATHDELIM + file):
               return False

        return True
               

    def hasGenbankFile(self):
        """ checks if genbank file has already been created"""
        output_annot_gbk= self.genbank_dir + PATHDELIM + self.sample_name +  '.gbk'
        return doesFileExist(output_annot_gbk)

    def hasSequinFile(self):
        """ checks if sequin file has already been created"""
        output_annot_sequin= self.output_results_sequin_dir + PATHDELIM + self.sample_name +  '.tbl'
        return doesFileExist(output_annot_sequin)
        


