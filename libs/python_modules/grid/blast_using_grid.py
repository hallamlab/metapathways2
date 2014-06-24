#!/usr/bin/python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


try:
   from optparse import make_option
   from os import makedirs, path, listdir, remove
   import os, sys, re, errno, shutil
   from glob import glob
   from datetime import date
   import traceback

   from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
   from libs.python_modules.parsers.parse  import parse_metapaths_parameters
   from libs.python_modules.utils.metapathways_utils import fprintf, printf,  WorkflowLogger, generate_log_fp
   from libs.python_modules.utils import *
   from libs.python_modules.grid.BlastService import BlastService
   from libs.python_modules.grid.GridParam import GridParam
   from libs.python_modules.grid.BlastBroker import BlastBroker
   import libs.python_scripts

except:
     print traceback.print_exc(10)
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
        #printf("%s", command[0])
        if command[4] == True:
           printf("%s", " Required")
        else:
           printf("%s", " Not Required")
    printf("\n")


def format_db(formatdb_executable, seqType, refdb_sequence_file, algorithm):
     if algorithm=='BLAST':
         cmd='%s -dbtype %s -in %s' %(formatdb_executable, seqType, refdb_sequence_file)

     if algorithm=='LAST':
         dirname = os.path.dirname(refdb_sequence_file)    
         cmd='%s -p -c %s  %s' %(formatdb_executable, refdb_sequence_file, refdb_sequence_file)
         

     result= getstatusoutput(cmd)
     if result[0]==0:
        return True
     else:
        return False

# decide if a command should be run if it is overlay,
# when results are alread computed decide not to run
def shouldRunStep(run_type, expected_output):
    if  run_type =='overlay'  and  path.exists(expected_output):
        return False
    else:
        return True 
    #print enable_flag


def formatted_db_exists(dbname):
    fileList = glob(dbname) 
    if len(fileList) > 0:
       return True
    else: 
       return False

def check_if_db_exists(dbname):
    if path.exists(dbname):
       return True
    else: 
       return False

def  make_sure_map_file_exists(dbmapfile):
    if not doFilesExist( [dbmapfile ] ):
         print 'WARNING: ' + 'Creating the database map file'
         fullRefDbName = re.sub(r'-names.txt','',dbmapfile)
         mapfile = open(dbmapfile,'w')
         fullRefDbFile = open(fullRefDbName,'r')
         for line in fullRefDbFile:
             if re.match(r'>', line):
                 fprintf(mapfile, "%s\n",line.strip())
         mapfile.close()
         fullRefDbFile.close()


#gets the parameter value from a category as specified in the 
# parameter file
def get_parameter(config_params, category, field, default = None):
    if config_params == None:
      return default

    if category in config_params:
        if field in config_params[category]:
            return config_params[category][field]
        else:
            return default
    return default


# parameter file
def get_make_parameter(config_params,category, field, default = False):
    if category in config_params:
        if field in config_params[category]:
            return config_params[category][field]
        else:
            return default
    return default

def get_pipeline_steps(steps_log_file):
    #print steps_log_file
    try:
       logfile = open(steps_log_file, 'r')
    except IOError:
       print "Did not find " + logfile + "!" 
       print "Try running in \'complete\' run-type"
    else:
       lines = logfile.readlines()
#       for line in lines:
#           print line
    pipeline_steps = None
    return pipeline_steps



# This function reads the pipeline configuration file and sets the 
# paths to differenc scripts and executables the pipeline call
def read_pipeline_configuration( file ):
    try:
       configfile = open(file, 'r')
    except IOError:
       print "Did not find pipeline config " + file + "!" 
    else:
       lines = configfile.readlines()

    config_settings = {}
    for line in lines:
        if not re.match("#",line) and len(line.strip()) > 0 :
           line = line.strip()
           line = re.sub('\t+', ' ', line)
           line = re.sub('\s+', ' ', line)
           line = re.sub('\'', '', line)
           line = re.sub('\"', '', line)
           fields = re.split('\s', line)
           if not len(fields) == 2:
              print """     The following line in your config settings files is set set up yet"""
              print """     Please rerun the pipeline after setting up this line"""
              print """ Error ine line :""" + line
              sys.exit(-1);
              
#           print fields[0] + "\t" + fields[1]
           if PATHDELIM=='\\':
              config_settings[fields[0]] = re.sub(r'/',r'\\',fields[1])   
           else:
              config_settings[fields[0]] = re.sub(r'\\','/',fields[1])   

           
    config_settings['METAPATHWAYS_PATH'] = config_settings['METAPATHWAYS_PATH'] + PATHDELIM
    config_settings['REFDBS'] = config_settings['REFDBS'] + PATHDELIM
    
    #check_config_settings(config_settings, file);

    return config_settings


# has the results to use
def hasResults(expected_output):
    if  path.exists(expected_output):
        return True
    else:
        return False


# has the results to use
def hasResults1(dir , expected_outputs):
    if  doFilesExist(expected_outputs, dir =  dir):
        return True
    else:
        return False


# if the directory is empty then there is not precomputed results
# and so you should decide to run the command
def shouldRunStepOnDirectory(run_type, dirName):
    dirName = dirName + PATHDELIM + '*'
    files = glob(dirName)
    if len(files)==0:
      return True
    else:
      return False

# if the command is "redo" then delete all the files
# in the folder and then delete the folder too
def removeDirOnRedo(command_Status, origFolderName):
    if command_Status=='redo' and path.exists(origFolderName) :
       folderName = origFolderName + PATHDELIM + '*'
       files = glob(folderName)
       for  f in files:
         remove(f)
       if path.exists(origFolderName): 
         shutil.rmtree(origFolderName)

# if the command is "redo" then delete the file
def removeFileOnRedo(command_Status, fileName):
    if command_Status=='redo' and path.exists(fileName) :
        remove(fileName)
        return True
    else:
        return False


# remove all the files in the directory on Redo
def cleanDirOnRedo(command_Status, folderName):
    if command_Status=='redo':
       cleanDirectory(folderName)


# remove all the files in the directory
def cleanDirectory( folderName):
    folderName = folderName + PATHDELIM + '*'
    files = glob(folderName)
    for  f in files:
       remove(f)

# if folder does not exist then create one
def checkOrCreateFolder( folderName ):
    if not path.exists(folderName) :
        makedirs(folderName)
        return False
    else:
        return True

#does the file Exist?
def doFilesExist( fileNames, dir="" ):
    for fileName in fileNames:
       file = fileName
       if dir!='':
         file = dir + PATHDELIM + fileName
       if not path.exists(file):
          return False
    return True


# check if all of the metapaths_steps have 
# settings from the valid list [ yes, skip stop, redo]

def  checkParam_values(allcategorychoices, parameters):
     for category in allcategorychoices:
        for choice in allcategorychoices[category]:
           if choice in parameters: 
#             print choice + " " + parameters[choice]  
#             print allcategorychoices[category][choice] 

             if not parameters[choice] in allcategorychoices[category][choice]:
                 print "ERROR: Incorrect setting in your parameter file"
                 print "       for step " + choice + " as  " + parameters[choices]
                 sys.exit(0)

def checkMetapaths_Steps(config_params):
     choices = { 'metapaths_steps':{}, 'annotation':{}, 'INPUT':{} }

     choices['INPUT']['format']  = ['fasta', 'gbk_unannotated', 'gbk_annotated', 'gff_unannotated', 'gff_annotated']

     choices['annotation']['algorithm'] =  ['last', 'blast'] 

     choices['metapaths_steps']['PREPROCESS_FASTA']   = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['ORF_PREDICTION']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['GFF_TO_AMINO']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['FILTERED_FASTA']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['COMPUTE_REFSCORE']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['BLAST_REFDB'] = ['yes', 'skip', 'stop', 'redo', 'grid']
     choices['metapaths_steps']['PARSE_BLAST'] = ['yes', 'skip', 'stop', 'redo']
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


     if config_params['metapaths_steps']:
        checkParam_values(choices, config_params['metapaths_steps'])


def isValidInput(output_dir, samples_and_input, dbs, gridSettings, config_settings,  messagelogger): 

    if not doesFolderExist(output_dir):   
        messagelogger.write("ERROR: Output folder \"%s\" does not exist!\n" %(output_dir)) 
        return False
    else:
        messagelogger.write("OK: Output folder \"%s\" exists!\n" %(output_dir)) 

    for sample, inputfile in samples_and_input.iteritems():
       if not doesFolderExist(output_dir + PATHDELIM + sample):   
          messagelogger.write("ERROR: Sample folder  \"%s\" for sample \"%s\" not found!\n" %(output_dir + PATHDELIM + sample, sample)) 
          return False
       else:
          messagelogger.write("OK: Sample folder  \"%s\" for sample \"%s\" exists!\n" %(output_dir + PATHDELIM + sample, sample)) 

       if not doesFileExist(inputfile):   
          messagelogger.write("ERROR: Input file \"%s\" does not exist!\n" %(inputfile)) 
          return False
       else:
          messagelogger.write("OK: Input file \"%s\" exists!\n" %(inputfile)) 

    for db in  dbs:
       if not doesFileExist(config_settings['REFDBS'] + PATHDELIM + 'functional' + PATHDELIM + db):   
          messagelogger.write("ERROR: Database file \"%s\" does not exist!\n" %(config_settings['REFDBS'] + PATHDELIM + db)) 
          return False
       else:
          messagelogger.write("OK: Database file \"%s\" found!\n" %(config_settings['REFDBS'] + PATHDELIM + db)) 

    for gridsetting in gridSettings:
         if 'server' in gridsetting:
            if not 'user' in gridsetting:
                 messagelogger.write("ERROR: User in grid servers \"%s\" not specified!\n" %(gridsetting['server'])) 
                 return False
         else:
                 messagelogger.write("OK: Specification for Grid\"%s\" with user \"%s\" found!\n" %(gridsetting['server'], gridsetting['user'])) 
    return True

#################################################################################
###########################  BLAST ##############################################
#################################################################################
def blast_in_grid(input_files, output_dir, config_params, metapaths_config, config_file, run_type):

    algorithm = get_parameter(config_params, 'annotation', 'algorithm', default='BLAST').upper()
    messagelogger = WorkflowLogger(generate_log_fp(output_dir, basefile_name='metapathways_messages', suffix='txt'),\
                    open_mode='w')

    command_Status=  get_parameter(config_params,'metapaths_steps','BLAST_REFDB')

    config_settings = read_pipeline_configuration( config_file )

#   preprocessed_dir = output_dir + PATHDELIM + "preprocessed" + PATHDELIM
    orf_prediction_dir =   "orf_prediction"  
#   genbank_dir =  output_dir + PATHDELIM + "genbank"  + PATHDELIM
    output_run_statistics_dir = output_dir + PATHDELIM + "run_statistics"  +PATHDELIM
    blast_results_dir =  output_dir +  PATHDELIM + "blast_results"  + PATHDELIM
    output_results = output_dir + PATHDELIM + "results" + PATHDELIM 
    #---

    # create the sample and input pairs 
    samples_and_input = {}
    for input_file in input_files:
       sample_name = re.sub(r'[.][a-zA-Z]*$','',input_file)
       sample_name = path.basename(sample_name)
       sample_name = re.sub('[.]','_',sample_name)
       samples_and_input[sample_name] =  output_dir + PATHDELIM + sample_name + PATHDELIM + orf_prediction_dir + PATHDELIM +  sample_name + ".qced.faa"   
    
    

    # BLAST THE ORFs AGAINST THE REFERENCE DATABASES  FOR FUNCTIONAL ANNOTATION
    dbstring = get_parameter(config_params, 'annotation', 'dbs', default=None)
    dbs= dbstring.split(",")

    #parse the grid settings from the param file
    gridEnginePATTERN = re.compile(r'(grid_engine\d+)')
    trueOrYesPATTERN = re.compile(r'^[yYTt]')

    gridSettings = []
    for key in config_params:
       match = gridEnginePATTERN.match(key)
       if match ==None:
           continue
       if 'active' in config_params[key]:
           trueOrYes =  trueOrYesPATTERN.match(config_params[key]['active'])
           if trueOrYes:  # this grid is inactive
               # proceed with adding the grid
               match = gridEnginePATTERN.match(key)
               if match:
                  gridSettings.append(config_params[key])

    
    if not isValidInput(output_dir, samples_and_input, dbs, gridSettings, config_settings = config_settings,\
         messagelogger = messagelogger): 
       sys.exit(0)
       
    blastbroker = BlastBroker(messagelogger) # setup the broker with a message logger
    blastbroker.setBaseOutputFolder(output_dir)  #set up the output folder 
    blastbroker.addSamples(samples_and_input)   # add the samples and the input files
    
    # add databases against the samples
    for sample in samples_and_input:
       for db in dbs:
          blastbroker.addDatabase(sample, db)
       blastbroker.addAlgorithm(sample, algorithm)   # add the algorithms
       
    # setup services and add them to the Broker 
    for gridsetting in gridSettings:
        gridsetting['messagelogger']=messagelogger
        gridsetting['MetaPathwaysDir']=config_settings['METAPATHWAYS_PATH']
        gridsetting['base_output_folder']=blastbroker.base_output_folder
        gridsetting['blast_db_folder']=config_settings['REFDBS'] + PATHDELIM + 'functional'

        try:
          blastservice = BlastService(gridsetting)
        except:
          print traceback.format_exc(10)

        blastbroker.addService(blastservice)

    # create the work space folders
    if  blastbroker.are_working_folders_available():
       messagelogger.write("STATUS: Local working folders for Grid found!\n")
    elif blastbroker.create_working_folders():
       messagelogger.write("OK: Successfully created the grid related local working folders!\n")
    else:
       messagelogger.write("ERROR: Cannot create the grid working folders!\n")
       messagelogger.write("ERROR: Exiting blast in grid mode!\n")
       return

    
    # check if the input files are already split
    messagelogger.write("STATUS: Checking if input files are already split!\n")
#    for s in blastbroker.getSamples():
#       if not blastbroker.doesValidSplitExist(s):
#          messagelogger.write("STATUS: Did not find any previously split files for sample \"%s\"!\n" %(s))
#          if not blastbroker.splitInput(s): #if not then split
#             messagelogger.write("ERROR: Cannot split the files for some or all of the samples!\n")
#             sys.exit(0)
#          else:
#             messagelogger.write("SUCCESS: Successfully split the files for some or all of the samples!\n")
#       else:
#          messagelogger.write("OK: Found previously split files for sample \"%s\"!\n" %(s))
#           
    messagelogger.write("STATUS: Competed checks for file splits!\n")

    batch_size = int(get_parameter(config_params, 'grid_submission', 'batch_size', default=1000))
    blastbroker.setBatchSize(batch_size)
    
    
    # check if the input files are already split
    for s in blastbroker.getSamples():
       if not blastbroker.doesValidSplitExist(s):
          messagelogger.write("STATUS: Did not find any previously split files for sample \"%s\"!\n" %(s))
          if not blastbroker.splitInput(s): #if not then split
             print ("ERROR: Cannot split the files for some or all of the samples!\n")
             messagelogger.write("ERROR: Cannot split the files for some or all of the samples!\n")
             sys.exit(0)
          else:
             messagelogger.write("SUCCESS: Successfully split the files for some or all of the samples!\n")
       else:
          messagelogger.write("OK: Found previously split files for sample \"%s\"!\n" %(s))
           
    # load the list of splits
    blastbroker.load_list_splits()
    messagelogger.write("SUCCESS: Successfully loaded the list of file splits!\n")
    
    # create the databse and split combinations as jobs for each sample
    blastbroker.createJobs(redo=False)
    messagelogger.write("SUCCESS: Successfully created the (split, database) pairs!\n")
    
    # make sure you loaded the latest job lists on file
    blastbroker.load_job_lists()
    messagelogger.write("SUCCESS: Successfully recovered the old/existing job list!\n")

    # for each sample load the submitted and completed lists
    # and compute the loadper Server
    blastbroker.load_job_status_lists()
    messagelogger.write("SUCCESS: Successfully loaded the status of the jobs!\n")

    blastbroker.compute_performance()

    try:
       blastbroker.compute_server_loads()
    except:
       print traceback.format_exc(10)

    #print blastbroker.list_jobs_submitted
    #print blastbroker.list_jobs_completed
    #blastbroker.launch_AWS_grid()
    blastbroker.setupStatsVariables()

    messagelogger.write("STATUS: Getting ready to submit jobs to the servers!\n")
    blastbroker.Do_Work()
    #blastbroker.stop_AWS_grid()

    blastbroker.Delete_Remote_Directories()
    

    #print output_dir    
    #print samples_and_input
    #print dbs
    #print gridSettings
    
    message = "\n6. Blasting using Grid ORFs against reference database - "


    #################################
#   Now call the commands

