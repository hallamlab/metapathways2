#!/usr/bin/env python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


"""Contains general utility code for the metapaths project"""

try:
    from shutil import rmtree
    from optparse import make_option
    from parameters import *
    from configuration import *
    from tools import *
    from libs.python_modules.utils.sysutil import pathDelim, getstatusoutput
    from libs.python_modules.utils.utils import *
    from os import path, _exit, rename
except:
    print "Cannot load some modules"
    sys.exit(0)
   
PATHDELIM = pathDelim()



def staticDiagnose(configs, params, logger = None ):
    """
    Diagnozes the pipeline basedon the configs and params for  
    binaries, scripts and resources               
    """

    """ makes sure that the choices in  parameter file are valid """
    errors = checkParams(params, log=logger)
    if errors:
       print 'errors'
       return 

    """ Get the configurations for the executables/scripts and databasess """
    _configuration = Configuration() 
    configuration = _configuration.getConfiguration()
    

    """ the place holders for the tools required to make the run """
    _tools = Tools()
    tools = _tools.getTools()

    """ load the actual executables """
    executables = matchToolsFromConfigs(configs, tools, logger )

    """ make sure all the executables exist """
    executablesExist(executables, logger)

    parameters = Parameters()

    #print  parameters.getRunSteps()
    """ check if the required set of executables exists """
    missingList = checkRequiredExecutables(parameters.getRunSteps(), _tools, params, logger)

    """ check if the required standard databases exists """

#    print  parameters.getRunSteps( activeOnly = True)
    checkForRequiredDatabases(tools, params, configs, 'functional',  logger = logger)

    checkForRequiredDatabases(tools, params, configs, 'taxonomic',  logger = logger)

def checkForRequiredDatabases(tools, params, configs, dbType, logger =None):
    """ checks the 
          -- database folder structure
          -- checks for raw sequences
          -- checks for formatted sequences
          -- formats if necessary 
    """
    
    if dbType=='functional':
       dbstring = get_parameter(params, 'annotation', 'dbs', default=None)
       _algorithm = get_parameter(params, 'annotation', 'algorithm', default=None)


    if dbType=='taxonomic':
       dbstring = get_parameter(params, 'rRNA', 'refdbs', default=None)

    dbs= [x.strip() for x in dbstring.split(",") ]
    refdbspath  = configs['REFDBS']

    """ checks refdb path """
    if not check_if_refDB_path_valid(refdbspath, logger = None):
        return False
        
    """ checks raw sequences for dbtype functional/taxonimic """
    if isRefDBNecessary(params, dbType):
        if not check_for_raw_sequences(dbs, refdbspath, dbType,  logger = logger):
          return False

        for db in dbs:
           algorithm = ""
           print dbType, "dbtype"
           if dbType=='taxonomic':
               algorithm = 'BLAST'
               seqType = 'nucl'
           elif dbType=='functional':
               algorithm= _algorithm
               seqType = 'prot'
           else:
               algorithm = None

           """ is db formatted ? """
           if not isDBformatted(db, refdbspath, dbType, seqType,   algorithm):
              """ if note formatted then format it """
              print 'Formatting DB ' , db
              if not formatDB(tools, db, refdbspath, seqType, dbType, algorithm, logger = logger):
                 return False


def formatDB(tools, db, refdbspath, seqType, dbType, algorithm, logger = None):
     """ Formats the sequences for the specified algorithm """

     formatdb_executable = tools['BLAST_REFDB']['exec']['BLAST']['FORMATDB_EXECUTABLE']
     if seqType=='nucl':
            formatdb_executable = tools['BLAST_REFDB']['exec']['BLAST']['FORMATDB_EXECUTABLE']
     if seqType=='prot':
        if algorithm=='LAST':
            formatdb_executable = tools['BLAST_REFDB']['exec']['LAST']['LASTDB_EXECUTABLE']
        if algorithm=='BLAST':
            formatdb_executable = tools['BLAST_REFDB']['exec']['BLAST']['FORMATDB_EXECUTABLE']


     formatted_db = refdbspath + PATHDELIM + dbType + PATHDELIM + 'formatted'  + PATHDELIM + db
     raw_sequence_file = refdbspath + PATHDELIM + dbType + PATHDELIM + db

     _temp_formatted_db  =  formatted_db+ "__temp__"

     """ format with 4GB file size """
     if algorithm=='BLAST':
         cmd='%s -dbtype %s -max_file_sz 4294967296  -in %s -out %s' %(formatdb_executable, seqType, raw_sequence_file, _temp_formatted_db)

     if algorithm=='LAST':
         # dirname = os.path.dirname(raw_sequence_file)    
         cmd='%s -s 4G -p -c %s  %s' %(formatdb_executable, _temp_formatted_db, raw_sequence_file)


     print cmd
     result= getstatusoutput(cmd)
     temp_fileList = glob(_temp_formatted_db + '*') 

     try: 
        for tempFile in temp_fileList:
           file = re.sub('__temp__','', tempFile)
           rename(tempFile, file);
     except:
        return False


     if result[0]==0:
        return True 
     else:
        return False


def isRefDBNecessary(params, dbType ):
    """ decide yes or no based on the params settings yes or redo """
    if dbType=="functional": 
        status = get_parameter(params, 'metapaths_steps', 'BLAST_REFDB', default=None)
        if status in [ 'yes', 'redo' ]:
           return True

    if dbType=="taxonomic": 
        status = get_parameter(params, 'metapaths_steps', 'SCAN_rRNA', default=None)
        if status in [ 'yes', 'redo' ]:
           return True

    return False


def isDBformatted(db, refdbspath, dbType, seqType,  algorithm):
    """ check if the DB is formatted """
    """Checks if the formatted database for the specified algorithm exits """
    dbPath = refdbspath + PATHDELIM + dbType + PATHDELIM + 'formatted'  
    dbname = dbPath + PATHDELIM + db 
    print seqType, algorithm
    suffixes = getSuffixes(algorithm, seqType) 

    print algorithm, suffixes
    if not suffixes :
       return False

    for suffix in suffixes:
       allfileList = glob(dbname + '*.' + suffix)

       fileList = []
       tempFilePattern = re.compile(r''+ dbname + '\d*.' + suffix +'$');


       for aFile in allfileList:
           searchResult =  tempFilePattern.search(aFile)
           if searchResult:
             fileList.append(aFile)

       if len(fileList)==0 :
          eprintf("ERROR :  if formatted correctely then expected the files with pattern %s\n", dbname + suffix)
          return False

    return True

def check_if_refDB_path_valid(refdbspath, logger = None):
    """ it checks for the validity of the refdbs path structure 
       refdbpath  /functional
                      /formatted
                  /tanxonomic
                      /formatted
    """

    status = True
    if not doesFolderExist(refdbspath):
        print "does not exist "
        return False

    dbTypes = [ 'functional', 'taxonomic' ]  

    if not doesFolderExist(refdbspath):
        print "does not exist "
        return False

    """ now check if respective dbtype folders are available """
    status = True
    for dbType in dbTypes:
       if not doesFolderExist(refdbspath + PATHDELIM + dbType):
          print dbType + " folder " + sQuote(refdbspath + PATHDELIM + dbType) +  " does not exist "
          status = False

    if status == False: 
       return status;

    """ now check if path to drop the formatted dbs are available """
    for dbType in dbTypes:
       if not doesFolderExist(refdbspath + PATHDELIM + dbType + PATHDELIM + 'formatted'):
          print dbType + " folder " + sQuote(refdbspath + PATHDELIM + dbType + PATHDELIM + 'formatted') +  " does not exist "
          status = False

    return status


def check_for_raw_sequences(dbs, refdbspath, dbType,  logger = None):
    """ check for the raw sequence file """
    status = True
    for db in dbs:
       fullPath =  refdbspath + PATHDELIM + dbType + PATHDELIM +  db 
       if not doesFolderExist(fullPath):
            print dbType + " raw seq does not exist in " + fullPath
            status = False

    return status 
    

def get_parameter(params, category, field, default = None):
    """  gets the parameter value from a category 
       as specified in the  parameter file """

    if params == None:
      return default

    if category in params:
        if field in params[category]:
            return params[category][field]
        else:    
            return default
    return default



def checkRequiredExecutables(steps, tools, params, logger = None):
    """  check the required executables in the steps """
    missingList = []
    for step in steps:
       missingList +=  executablesExist(tools.getExecutables(step, params), logger)

    return missingList


def executablesExist( executables, logger = None ):
    missingList = []
    for name, script in executables.iteritems():
      if path.exists(script):
           pass
      else:
           print "ERROR ", name, script
           missingList.append(script)

    return missingList

def matchToolsFromConfigs(configs, tools, logger = None ):
    """ iterate througs each of the configs item and fill it to 
        to the actual value in the tools 
    """
    executables = {}
    """ iterate through each config key """
    for config_key, config_value in configs.iteritems():
       for param_step, placeHolder in tools.iteritems(): 
         if param_step in tools:
           for script in tools[param_step]['exec']: 

               """ does not have alternatives """
               if not type(tools[param_step]['exec'][script]) is dict: 
                   if config_key == script: 
                      tools[param_step]['exec'][script]=config_value
                      executables[script]= config_value
               else: 
                   """ go a level deeper """
                   for sub_script in tools[param_step]['exec'][script]: 
                     if config_key == sub_script: 
                        tools[param_step]['exec'][script][sub_script]=config_value
                        executables[sub_script]= config_value

    return executables


def getRequiredTools(params, configs,  tools, configuration):
    if not 'metapaths_steps' in params:
       return None

    for key, value in params['metapaths_steps'].iteritems():
        if value in [  'skip', 'redo', 'yes' ]:

           if not key in tools:
              #print "ERROR : " + key + " is missing in class Tools file!"
              continue

           if type(tools[key]['exec'][key]) is dict: 
              for execname in tools[key]['exec']: 
                 print execname
          
           #if not tools[key]['exec'] in configs:
              #print "ERROR : Exec in " + tools[key]['exec']  + " is missing in class Configuration file!"
           #   tools[key]['exec'] = None
           #   continue

           #print configs[tools[key]['exec']]
           #tools[key]['exec'] = [ configs[tools[key]['exec']] ]

    print 'tools' 
    print tools

def _checkParams(key, params, paramsAccept, log= None, errors= None):

    """  make sure that every parameter in the params is valid recursively 

     This is initialed by the checkParams() function 

     store the erros in the erros dictionary 
     """


    """ if not level to go deeper  then the leaves of the dict are reached"""
    if not type(params) is dict: 
        if  type(paramsAccept) is dict:
          try:
             if not params in paramsAccept:
                if errors!=None:
                    errors[key] = False
          except:
                pass
        return

    """  make sure that every parameter in the params is valid recursively """
    for key, value in params.iteritems(): 
        if type(paramsAccept) is dict:
           if len(key) and key in paramsAccept:
               _checkParams(key, params[key], paramsAccept[key], log= log, errors = errors)

def checkParams(params, log= None):
    """ makes sure that all the params provides are valid or acceptable """
    """ when the choices are not any of the acceptable 
    values then it is considered erroneous"""

    _paramsAccept = Parameters()
    paramsAccept = _paramsAccept.getAcceptableParameters() 
    errors = {}

    for key, value in params.iteritems(): 
       _checkParams(key, params, paramsAccept, log= log, errors = errors)

    return errors


def getSuffixes(algorithm, seqType) :
    """ Get the suffixes for the right algorithm with the right 
        sequence type 
    """

    suffixes = {}
    suffixes['LAST'] = {}
    suffixes['BLAST'] = {}
    suffixes['BLAST']['nucl'] = ['nhr', 'nsq', 'nin']
    suffixes['BLAST']['prot'] = ['phr', 'psq', 'pin']

    suffixes['LAST']['prot'] = [ 'des', 'sds', 'suf', 'bck', 'prj', 'ssp', 'tis' ]
    suffixes['LAST']['nucl'] = []

    if not algorithm in suffixes:
        return None


    if not seqType in suffixes[algorithm]:
        return None

    return suffixes[algorithm][seqType]

