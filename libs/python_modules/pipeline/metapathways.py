#!/usr/bin/python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


try:
   from optparse import make_option
   from os import makedirs,  path, listdir, remove, rename, _exit
   import os, sys, errno, shutil, re
   from glob import glob
   from datetime import date
   #from metapaths_utils  import pars[s._command_line_parameters

   from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
   #from libs.python_modules.utils.utils import *, hasInput, createFolderIfNotFound
   from libs.python_modules.utils.utils import *
   from libs.python_modules.parsers.parse  import parse_metapaths_parameters
   from libs.python_modules.pipeline.metapathways_pipeline import print_commands,  execute_tasks
   from libs.python_modules.pipeline.MetaPathways_gather_run_stats import MetaPathways_gather_run_stats
   from libs.python_modules.utils.metapathways_utils import fprintf, printf, eprintf, remove_existing_pgdb, exit_process, WorkflowLogger, generate_log_fp

   from libs.python_modules.pipeline.sampledata import *
   from libs.python_modules.pipeline.jobscreator import *
   from libs.python_modules.pipeline.commands import *
   import libs.python_scripts


except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)



PATHDELIM = pathDelim()


def copyFile(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

def dry_run_status( commands ):
    for command in commands:
        printf("%s", command[0])
        if command[4] == True:
           printf("%s", " Required")
        else:
           printf("%s", " Not Required")
    printf("\n")


def get_refdb_name( dbstring ):
    dbstring = dbstring.rstrip()
    dbstring = dbstring.lstrip()
    dbstring = dbstring.lower() 
    return dbstring



def format_db(formatdb_executable, seqType, raw_sequence_file, formatted_db,  algorithm):
     _temp_formatted_db  =  formatted_db+ "__temp__"

     """ format with 4GB file size """
     if algorithm=='BLAST':
         cmd='%s -dbtype %s --max_file_sz 4294967296  -in %s -out %s' %(formatdb_executable, seqType, raw_sequence_file, _temp_formatted_db)

     if algorithm=='LAST':
         # dirname = os.path.dirname(raw_sequence_file)    
         cmd='%s -s 4G -p -c %s  %s' %(formatdb_executable, _temp_formatted_db, raw_sequence_file)

     result= getstatusoutput(cmd)
     temp_fileList = glob(_temp_formatted_db + '*') 
     try:
        for tempFile in temp_fileList:
           file = re.sub('__temp__','', tempFile)
           rename( tempFile, file);

     except:
        return False

     if result[0]==0:
        return True
     else:
        return False


# convert an input gbk file to fna faa and gff file
def  convert_gbk_to_fna_faa_gff(input_gbk, output_fna, output_faa, output_gff, config_settings):
    cmd = "%s  -g %s --output-fna %s --output-faa %s --output-gff %s" %((config_settings['METAPATHWAYS_PATH'] \
                 + config_settings['GBK_TO_FNA_FAA_GFF']), input_gbk, output_fna, output_faa, output_gff) 
    return cmd

# convert an input gff file to fna faa and gff file
def  convert_gff_to_fna_faa_gff(inputs, outputs,  config_settings):
    cmd = "%s " %(config_settings['METAPATHWAYS_PATH']+ config_settings['GFF_TO_FNA_FAA_GFF'])
    for  source, target in zip(inputs, outputs):
       cmd += ' --source ' + source + ' --target ' + target
    return cmd



def  make_sure_map_file_exists(config_settings, dbname, globallogger = None):
    dbmapFile = config_settings['REFDBS'] + PATHDELIM + 'functional' + PATHDELIM + 'formatted' + PATHDELIM + dbname + "-names.txt"
    seqFilePath = config_settings['REFDBS'] + PATHDELIM + 'functional' + PATHDELIM + dbname
    if not doFilesExist( [dbmapFile ] ):
         eprintf("WARNING: Trying to create database map file for %s\n", dbname)
         if globallogger!= None:
            globallogger.write("WARNING: Trying to create database map file for %s\n" %( dbname) )

         if not doFilesExist( [seqFilePath] ):
            eprintf("ERROR : You do not even have the raw sequence for Database  %s to format!\n", dbname)
            eprintf("      : Make sure you have the file %s\n", seqFilePath)

            if globallogger!= None:
               globallogger.write("ERROR \t You do not even have the raw sequence for Database  %s to format!\n" %( dbname))
               globallogger.write("Make sure you have the file %s\n" %( seqFilePath))

            exit_process()

         mapfile = open(dbmapFile,'w')
         seqFile = open(seqFilePath,'r')
         for line in seqFile:
             if re.match(r'>', line):
                 fprintf(mapfile, "%s\n",line.strip())
         seqFile.close()
         mapfile.close()

    return dbmapFile


# create the command to make the MLTreeMap Images
def create_MLTreeMap_Imagemaker(mltreemap_image_output, mltreemap_final_outputs, config_settings):
    executable_path = config_settings['MLTREEMAP_IMAGEMAKER']
    if not path.isfile( executable_path):
        executable_path = config_settings['METAPATHWAYS_PATH'] + executable_path
    cmd= "%s -i %s -o %s -m a"  %(executable_path, mltreemap_final_outputs, mltreemap_image_output) 
    return cmd


# gather mltreemap calculations 
def create_MLTreeMap_Hits(mltreemap_output_dir, output_folder, config_settings):
    cmd= "%s -i %s -o %s"  %(config_settings['METAPATHWAYS_PATH'] + config_settings['MLTREEMAP_HITS'], mltreemap_output_dir, output_folder +PATHDELIM + 'sequence_to_cog_hits.txt') 
    return cmd


#gets the parameter value from a category as.ecified in the 
# parameter file
def get_parameter(params, category, field, default = None):
    if params == None:
      return default

    if category in params:
        if field in params[category]:
            return params[category][field]
        else:
            return default
    return default


# parameter file
def get_make_parameter(params,category, field, default = False):
    if category in params:
        if field in params[category]:
            return params[category][field]
        else:
            return default
    return default

def get_pipeline_steps(steps_log_file):
    try:
       logfile = open(steps_log_file, 'r')
    except IOError:
       eprintf("Did not find %s!\n", logfile) 
       eprintf("Try running in \'complete\' run-type\n")
    else:
       lines = logfile.readlines()

    pipeline_steps = None
    return pipeline_steps


def write_run_parameters_file(fileName, parameters):
    try:
       paramFile = open(fileName, 'w')
    except IOError:
       eprintf("Cannot write run parameters to file %s!\n", fileName)
       exit_process("Cannot write run parameters to file %s" %(fileName) )

#       16s_rRNA      {'min_identity': '40', 'max_evalue': '0.000001', 'min_bitscore': '06', 'refdbs': 'silva_104_rep_set,greengenes_db_DW'}
    paramFile.write("\nRun Date : " + str(date.today()) + " \n")

    paramFile.write("\n\nNucleotide Quality Control parameters[s.n")
    paramFile.write( "  min length" + "\t" + str(parameters['quality_control']['min_length']) + "\n")

    paramFile.write("\n\nORF prediction parameters[s.n")
    paramFile.write( "  min length" + "\t" + str(parameters['orf_prediction']['min_length']) + "\n")
    paramFile.write( "  algorithm" + "\t" + str(parameters['orf_prediction']['algorithm']) + "\n")


    paramFile.write("\n\nAmino acid quality control and annotation parameters[s.n")
    paramFile.write( "  min bit score" + "\t" + str(parameters['annotation']['min_score']) + "\n")
    paramFile.write( "  min seq length" + "\t" + str(parameters['annotation']['min_length']) + "\n")
    paramFile.write( "  annotation reference dbs" + "\t" + str(parameters['annotation']['dbs']) + "\n")
    paramFile.write( "  min BSR" + "\t" + str(parameters['annotation']['min_bsr']) + "\n")
    paramFile.write( "  max evalue" + "\t" + str(parameters['annotation']['max_evalue']) + "\n")

    paramFile.write("\n\nPathway Tools parameters[s.n")
    paramFile.write( "  taxonomic pruning " + "\t" + str(parameters['ptools_settings']['taxonomic_pruning']) + "\n")

    paramFile.write("\n\nrRNA search/match parameters[s.n")
    paramFile.write( "  min identity" + "\t" + str(parameters['rRNA']['min_identity']) + "\n")
    paramFile.write( "  max evalue" + "\t" + str(parameters['rRNA']['max_evalue']) + "\n")
    paramFile.write( "  rRNA reference dbs" + "\t" + str(parameters['rRNA']['refdbs']) + "\n")

    paramFile.close()


# checks if the necessary files, directories  and executables really exis.or not
def check_config_settings(config_settings, file, globalerrorlogger = None):
   essentialItems= ['METAPATHWAYS_PATH', 'EXECUTABLES_DIR', 'RESOURCES_DIR']
   missingItems = []

   for key, value in  config_settings.items():
      # make sure  MetaPathways directory is present
      if key in ['METAPATHWAYS_PATH' ]:
         if not path.isdir( config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 1.Currently it is set to \"%s\"\n",  config_settings[key] )  

            if globalerrorlogger!=None:
               globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n"  %(key, file))  
               globalerrorlogger.write("       Currently it is set to \"%s\"\n" %(config_settings[key] )  )
            missingItems.append(key) 
         continue


      # make sure  REFDB directories are present
      if key in [ 'REFDBS' ]:
         if not path.isdir( config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 2.Currently it is set to \"%s\"\n", config_settings[key] )  
            if globalerrorlogger!=None:
                globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key,file))
                globalerrorlogger.write("Currently it is set to \"%s\"\n" %( config_settings[key]) )  
            missingItems.append(key) 
         continue

      # make sure EXECUTABLES_DIR directories are present
      if key in [ 'EXECUTABLES_DIR']:
         if not path.isdir( config_settings['METAPATHWAYS_PATH'] + PATHDELIM +  config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 3.Currently it is set to \"%s\"\n", config_settings[key] )  
            if globalerrorlogger!=None:
               globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key, file))  
               globalerrorlogger.write("Currently it is set to \"%s\"\n" %( config_settings[key] )) 
            missingItems.append(key) 
         continue

      # make sure RESOURCES_DIR directories are present
      if key in [ 'RESOURCES_DIR']:
         if not path.isdir( config_settings['METAPATHWAYS_PATH'] + PATHDELIM + config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 4.Currently it is set to \"%s\"\n",  config_settings['METAPATHWAYS_PATH'] + PATHDELIM + config_settings[key] )  
            print  config_settings['METAPATHWAYS_PATH'], config_settings[key] 
            if globalerrorlogger!=None:
               globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key, file))
               globalerrorlogger.write("Currently it is set to \"%s\"\n" %( config_settings[key]))  
            missingItems.append(key) 
         continue

      # make sure  MetaPaths directory is present
      if key in ['PYTHON_EXECUTABLE' , 'PATHOLOGIC_EXECUTABLE' ]:
         if not path.isfile( config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 5.Currently it is set to \"%s\"\n", config_settings[key] )  
            if globalerrorlogger!=None:
               globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key, file)) 
               globalerrorlogger.write("Currently it is set to \"%s\"\n" %( config_settings[key] ) )
            missingItems.append(key) 
         continue

      # ignore pgdb folder for now
      if key in ['PGDB_FOLDER' ]:
          continue
      
      # check if the desired file exists. if not, then print a message
      if not path.isfile( config_settings['METAPATHWAYS_PATH'] + PATHDELIM +  value)\
        and  not path.isfile( config_settings['METAPATHWAYS_PATH'] + PATHDELIM + config_settings['EXECUTABLES_DIR'] + PATHDELIM + value ) :
           eprintf("ERROR:Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
           eprintf("5.Currently it is set to \"%s\"\n", config_settings['METAPATHWAYS_PATH'] + value ) 
           if globalerrorlogger!=None:
              globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key, file) )
              globalerrorlogger.write("Currently it is set to \"%s\"\n" %(config_settings['METAPATHWAYS_PATH'] + value)) 
           missingItems.append(key) 
           continue
     
   stop_execution = False
   for item in missingItems:
      if item in essentialItems:
         eprintf("ERROR\t Essential field in setting %s is missing in configuration file!\n", item)
         if globalerrorlogger!=None:
            globalerrorlogger.write("ERROR\tEssential field in setting %s is missing in configuration file!\n" %(item))
         stop_execution = True

   if stop_execution ==True:
      eprintf("ERROR: Terminating execution due to missing essential  fields in configuration file!\n")
      if globalerrorlogger!=None:
         globalerrorlogger.write("ERROR\tTerminating execution due to missing essential  fields in configuration file!\n")
      exit_process()

   

# This function reads the pipeline configuration file and sets the 
# paths to differenc scripts and executables the pipeline call
def read_pipeline_configuration( file, globallogger ):
    patternKEYVALUE = re.compile(r'^([^\t\s]+)[\t\s]+\'(.*)\'')
    try:
       configfile = open(file, 'r')
    except IOError:
       eprintf("ERROR :Did not find pipeline config %s!\n", file) 
       globalerrorlogger.write("ERROR\tDid not find pipeline config %s!\n" %(file)) 
    else:
       lines = configfile.readlines()

    config_settings = {}
    for line in lines:
        if not re.match("#",line) and len(line.strip()) > 0 :
           line = line.strip()
           result = patternKEYVALUE.search(line)
           
           try:
              if len(result.groups()) == 2:
                 fields = result.groups()
              else:
                 eprintf("     The following line in your config settings files is not set up yet\n")
                 eprintf("     Please rerun the pipeline after setting up this line\n")
                 eprintf("     Error in line : %s\n", line)
                 globalerrorlogger(
                      "WARNING\t\n"+\
                      "     The following line in your config settings files is not set up yet\n"+\
                      "     Please rerun the pipeline after setting up this line\n"+\
                      "     Error in line : %s\n" %(line))

                 exit_process()
           except:
                 eprintf("     The following line in your config settings files is not set up yet\n")
                 eprintf("     Please rerun the pipeline after setting up this line\n")
                 eprintf("     Error ine line : %s\n", line)
                 globalerrorlogger(
                      "WARNING\t\n"+\
                      "     The following line in your config settings files is not set up yet\n"+\
                      "     Please rerun the pipeline after setting up this line\n"+\
                      "     Error in line : %s\n" %(line))
                 exit_process()
              
           if PATHDELIM=='\\':
              config_settings[fields[0]] = re.sub(r'/',r'\\',fields[1])   
           else:
              config_settings[fields[0]] = re.sub(r'\\','/',fields[1])   

           
    config_settings['METAPATHWAYS_PATH'] = config_settings['METAPATHWAYS_PATH'] + PATHDELIM
    config_settings['REFDBS'] = config_settings['REFDBS'] + PATHDELIM
    
    check_config_settings(config_settings, file, globallogger);
    config_settings['configuration_file'] = file

    return config_settings

#check for empty values in parameter settings 
def  checkMissingParam_values(params, choices, logger = None):
     reqdCategoryParams = { 
                            'annotation': {'dbs': False}, 
                            'orf_prediction':{}, 
                            'rRNA':{},
                            'metapaths_steps':{}
                         }

     success  = True
     for category in choices:
       for parameter in choices[category]:
         if (not params[category][parameter]) and\
            ( (category in reqdCategoryParams) and\
               (parameter in reqdCategoryParams[category]) and   reqdCategoryParams[category][parameter]) :
            print category, parameter
            print reqdCategoryParams
            print reqdCategoryParams[category]
            eprintf('ERROR: Empty parameter %s of type %s\n'  %(parameter, category))
            eprintf('Please select at least one database for %s\n'  %(category))
            if logger!=None:
               logger.write('ERROR\tEmpty parameter %s of type %s\n'  %(parameter, category))
               logger.write('Please select at least one database for %s\n'  %(category))
            success = False

     return success

# check if all of the metapaths_steps have 
# settings from the valid list [ yes, skip stop, redo]

def  checkParam_values(allcategorychoices, parameters, runlogger = None):
     for category in allcategorychoices:
        for choice in allcategorychoices[category]:
           if choice in parameters: 

             if not parameters[choice] in allcategorychoices[category][choice]:
                 logger.write('ERROR\tIncorrect setting in your parameter file')
                 logger.write('for step %s as %s' %(choice, parameters[choices]))
                 eprintf("ERROR: Incorrect setting in your parameter file" +\
                         "       for step %s as %s", choice, parameters[choices])
                 exit_process()

def checkMetapathsteps(params, runlogger = None):
     choices = { 'metapaths_steps':{}, 'annotation':{}, 'INPUT':{} }

     choices['INPUT']['format']  = ['fasta', 'gbk_unannotated', 'gbk_annotated', 'gff_unannotated', 'gff_annotated']

     choices['annotation']['algorithm'] =  ['last', 'blast'] 

     choices['metapaths_steps']['PREPROCESS_FASTA']   = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['ORF_PREDICTION']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['GFF_TO_AMINO']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['FILTERED_FASTA']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['COMPUTE_REFSCORE']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['BLAST_REFDB'] = ['yes', 'skip', 'stop', 'redo', 'grid']
     choices['metapaths_steps']['PARSE._BLAST'] = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['SCAN_rRNA']   = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['STATS_rRNA']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['ANNOTATE']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['PATHOLOGIC_INPUT']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['GENBANK_FILE']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['CREATE_SEQUIN_FILE']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['CREATE_REPORT_FILES']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['SCAN_tRNA']   = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['MLTREEMAP_CALCULATION']   = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['MLTREEMAP_IMAGEMAKER']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['PATHOLOGIC']  = ['yes', 'skip', 'stop', 'redo']


     if params['metapaths_steps']:
        checkParam_values(choices, params['metapaths_steps'], runlogger)

     checkparams = {}
     checkparams['annotation'] = []
     checkparams['annotation'].append('dbs') 

     if not checkMissingParam_values(params, checkparams, runlogger):
        exit_process("Missing parameters")


def  copy_fna_faa_gff_orf_prediction( source_files, target_files, config_settings) :

     for source, target in zip(source_files, target_files):  

         sourcefile = open(source, 'r')
         targetfile = open(target, 'w')
         sourcelines = sourcefile.readlines()
         for line in sourcelines:
            fprintf(targetfile, "%s\n", line.strip())

         sourcefile.close()
         targetfile.close()


#################################################################################
########################### BEFORE BLAST ########################################
#################################################################################
def run_metapathways(samplesData, output_dir, all_samples_output_dir, globallogger,\
                     command_line_params, params, metapaths_config, status_update_callback,\
                     config_file, run_type, config_settings = None, block_mode = False, runid = ""):

    jobcreator = JobCreator(params, config_settings)

    sorted_samplesData_keys = sorted(samplesData.keys())
    for input_file in sorted_samplesData_keys:
      s =  samplesData[input_file]
      jobcreator.addJobs(s, block_mode = block_mode)


    if block_mode:
       eprintf("==============  RUNNING STEPS IN BLOCK 0 ================\n")
       for input_file in sorted_samplesData_keys:
         s =  samplesData[input_file]
         s.stepslogger.printf("\n\n==============  BEGIN RUN " + s.sample_name + " " + runid + " BLOCK0 ================\n")
         sample_name_banner = "PROCESSING INPUT " + input_file
         eprintf('\n'+ '#'*len(sample_name_banner) + "\n")
         eprintf( '\n' + sample_name_banner +  ' [STEPS BLOCK 0] ' + '\n')
         try:
            pass
            execute_tasks(s, verbose = command_line_params['verbose'], block = 0)    
         except:
            pass

       for input_file in sorted_samplesData_keys:
         s =  samplesData[input_file]
         s.stepslogger.printf("\n\n==============  BEGIN RUN " + s.sample_name + " " + runid + " BLOCK1 ================\n")
         sample_name_banner = "PROCESSING INPUT " + input_file
         eprintf('\n' + '#'*len(sample_name_banner) + "\n")
         eprintf( '\n' + sample_name_banner +  ' [STEPS BLOCK 1] ' + '\n')
         try:
            pass
            execute_tasks(s, verbose = command_line_params['verbose'], block = 1)    
         except:
            pass

       for input_file in sorted_samplesData_keys:
         s =  samplesData[input_file]
         s.stepslogger.printf("\n\n==============  BEGIN RUN " + s.sample_name + " " + runid + " BLOCK2 ================\n")
         sample_name_banner = "PROCESSING INPUT " + input_file
         eprintf('\n' + '#'*len(sample_name_banner) + "\n")
         eprintf( '\n' + sample_name_banner +  ' [STEPS BLOCK 2] ' + '\n')
         try:
            pass
            execute_tasks(s, verbose = command_line_params['verbose'], block = 2)    
         except:
            pass


    else:
       for input_file in sorted_samplesData_keys:
         s =  samplesData[input_file]
         s.stepslogger.printf("\n\n==============  BEGIN RUN " + s.sample_name + " " + runid + "  ==================\n")
         sample_name_banner = "PROCESSING INPUT " + input_file
         eprintf('#'*len(sample_name_banner) + "\n")
         eprintf( '\n' + sample_name_banner + '\n')
         try:
            execute_tasks(s, verbose = command_line_params['verbose'], block = 0)    
         except:
            pass
         try:
            execute_tasks(s, verbose = command_line_params['verbose'], block = 1)    
         except:
            pass
         try:
            execute_tasks(s, verbose = command_line_params['verbose'], block = 2)    
         except:
            pass

    return


