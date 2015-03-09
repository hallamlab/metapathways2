from __future__ import division

__author__ = "Kishori M Konwar Niels W Hanson"
__copyright__ = "Copyright 2014, MetaPathways"
__credits__ = [""]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar Niels W Hanson"
__status__ = "Release"

#from  libs.starcluster.test import  teststarcluster as sctest
#import sys

try:
     import sys, traceback, re, inspect, signal, shutil 
     from os import makedirs, sys, listdir, environ, path, _exit
     #from commands import getstatusoutput
     from optparse import OptionParser
     
     from libs.python_modules.utils import metapathways_utils
     from libs.python_modules.utils.utils import *
     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, eprintf, halt_process, exit_process, WorkflowLogger, generate_log_fp
     from libs.python_modules.parsers.parse  import parse_metapaths_parameters, parse_parameter_file
     from libs.python_modules.pipeline.metapathways_pipeline import print_commands, print_to_stdout, no_status_updates
     from libs.python_modules.utils.sysutil import pathDelim
     from libs.python_modules.pipeline.metapathways import run_metapathways, get_parameter, read_pipeline_configuration
     from libs.python_modules.annotate import *
     from libs.python_modules.grid.blast_using_grid import blast_in_grid

     from libs.python_modules.diagnostics.parameters import *
     from libs.python_modules.diagnostics.diagnoze import *
     from libs.python_modules.pipeline.sampledata import *
except:
   print """ Could not load some user defined  module functions"""
   print """ Make sure your typed \"source MetaPathwaysrc\""""
   print """ """
   #print traceback.print_exc(10)
   sys.exit(3)


cmd_folder = path.abspath(path.split(inspect.getfile( inspect.currentframe() ))[0])

PATHDELIM =  str(pathDelim())

#print cmd_folder
#if not sys.platform.startswith('win'):
#    res =getstatusoutput('source  '+ cmd_folder +'/'+'.metapathsrc')
#    if( int(res[0])==0 ): 
#       print 'Ran ' + cmd_folder +'/'+'.metapathsrc ' + ' file successfully!'
#    else:
#       print 'Error : ' + res[1] 
#       print 'while running  ' + cmd_folder +'/'+'.metapathsrc ' + ' file!'

#sys.path.insert(0,cmd_folder + "/libs/python_modules/")
#sys.path.insert(1, cmd_folder + "/libs/")
#print sys.path

#config = load_config()
metapaths_config = """config/template_config.txt""";
metapaths_param = """config/template_param.txt""";

script_info={}
script_info['brief_description'] = """A workflow script for making PGDBs from metagenomic sequences"""
script_info['script_description'] = \
    """ This script starts a MetaPathways pipeline run. It requires an input directory of fasta or genbank files
    containing sequences to process, an output directory for results to be placed. It also requires the
    configuration files, template_config.txt and template_param.txt in the config/ directory, to be updated with the
    location of resources on your system.
    """
script_info['script_usage'] = []


usage=  sys.argv[0] + """ -i input_dir -o output_dir -p parameters.txt
For more options:  ./MetaPathways.py -h"""

parser = None
def createParser():
    global parser
    parser = OptionParser(usage)
    parser.add_option("-i", "--input_file", dest="input_fp",
                      help='the input fasta file/input dir [REQUIRED]')
    parser.add_option("-o", "--output_dir", dest="output_dir",
                      help='the input fasta file/input dir [REQUIRED]')
    parser.add_option('-p','--parameter_fp', dest="parameter_fp",
                       help='path to the parameter file [REQUIRED]')
    parser.add_option("-c", "--config_filer", dest="config_file",
                      help='pipeline_configuratin file [OPTIONAL,  default : \"MetaPathways/template_config.txt\"]')
    parser.add_option('-r','--run-type', dest="run_type", default='safe',
                       choices=['safe', 'overlay', 'overwrite','dry-run'], 
                       help= '\n(a) \'overwrite\' -- wipes out the previous runs with the same name\n'+
                             '\n(b)\'overlay\' -- recomputes the steps that are not present \n' +
                             '\n(d)\'safe\' -- safe mode does not run on an existing run folder\n')

    #ith out of order completion \ time-stamps in the \'workflow_log.txt\' 
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="print lots of information on the stdout [default]")
    
    parser.add_option("-b", "--block-mode",
                      action="store_true", dest="block_mode", default=True,
                      help="processes the samples by blocking the stages before and after functional search [default off]")

    parser.add_option("-d", "--delay", dest="delay", type='int',  default=0,
                      help="number of seconds to sleep once the run is done")
    
    
    parser.add_option("-P", "--print-only",
                      action="store_true", dest="print_only", default=False,
                      help="print only  the commands [default False]")
    
    parser.add_option("-n", "--ncbi-header", dest="ncbi_header", 
                      help="NCBI sequin submission parameters file" )
    
    parser.add_option("-s", "--subset", dest="sample_subset", action="append", default=[],
                      help="Processes only samples in the list  subset specified [ -s sample1 -s sample2 ]" )
    
    parser.add_option("--runid", dest="runid",  default="",
                      help="Any string to represent the runid [ default Empty string ]" )
    #parser.add_option("-s", "--ncbi-sbt-file", dest="ncbi_sbt", 
    #                  help="the NCBI sbt location created by the \"Create Submission Template\" form: http://www.ncbi.nlm.nih.gov/WebSub/template.cgi" )



def valid_arguments(opts, args):
    """ checks if the supplied arguments are adequate """
    if (opts.input_fp == None and opts.output_dir ==None )  or\
     opts.output_dir == None:
       return True
    else:
       return False

def derive_sample_name(filename):
    basename = path.basename(filename) 
    
    shortname = re.sub('[.]gbk$','',basename, re.IGNORECASE) 
    shortname = re.sub('[.](fasta|fas|fna|faa|fa)$','',shortname, re.IGNORECASE) 
    return shortname
    


def remove_unspecified_samples(input_output_list, sample_subset,  globalerrorlogger = None):
   """ keep only the samples that are specified  before processing  """

   shortened_names = {}
   input_sample_list = input_output_list.keys()
   for sample_name in input_sample_list:
      if not derive_sample_name(sample_name) in sample_subset and  sample_subset:
         del input_output_list[sample_name]



def check_for_error_in_input_file_name(shortname, globalerrorlogger=None):

    """  creates a list of  input output pairs if input is  an input dir """
    clean = True
    if not re.search(r'^[a-zA-Z]',shortname):
         eprintf("ERROR\tSample name %s must begin with an alphabet!\n",shortname)
         if globalerrorlogger:
            globalerrorlogger.printf("ERROR\tSample name %s must begin with an alphabet!\tConsider prefixing an alphabet to the front\n",shortname)
         clean = False

    if re.search(r'[.]',shortname):
         eprintf("ERROR\tSample name %s contains a '.' in its name!\n",shortname)
         if globalerrorlogger:
            globalerrorlogger.printf("ERROR\tSample name %s contains a '.' in its name!\n",shortname)
         clean = False

    if len(shortname)<2:
         eprintf("ERROR\tSample name %s is too short!\n",shortname)
         if globalerrorlogger:
             globalerrorlogger.printf("ERROR\tSample name %s is too short1\n",shortname)
         clean = False

    if clean:
         return clean

    errmessage = """Sample names before the  suffixes .fasta, .fas, .fna, .faa or .gbk, must  consist only of alphabets, digits and _; and should consist of at least two characters """
    eprintf("ERROR\t%s\n",errmessage)
    if globalerrorlogger:
        globalerrorlogger.printf("ERROR\t%s\n",errmessage)
    #    exit_process(errmessage + "Exiting!" + "\n", logger=globalerrorlogger)
    return False


def create_an_input_output_pair(input_file, output_dir,  globalerrorlogger=None):
    """ creates an input output pair if input is just an input file """
       
    input_output = {}

    if not re.search(r'.(fasta|fas|fna|faa|gbk|gff|fa)$',input_file, re.IGNORECASE):
       return input_output

    shortname = None 
    shortname = re.sub('[.]gbk$','',input_file, re.IGNORECASE) 
    shortname = re.sub('[.](fasta|fas|fna|faa|fa)$','',input_file, re.IGNORECASE) 
    #    shortname = re.sub('[.]gff$','',input_file, re.IGNORECASE) 

    shortname = re.sub(r'.*' + PATHDELIM ,'',shortname) 

    if  check_for_error_in_input_file_name(shortname, globalerrorlogger=globalerrorlogger):
       input_output[input_file] = path.abspath(output_dir) + PATHDELIM + shortname

    return input_output


def create_input_output_pairs(input_dir, output_dir,  globalerrorlogger=None):
    """  creates a list of  input output pairs if input is  an input dir """
    fileslist =  listdir(input_dir)

    gbkPatt = re.compile('[.]gbk$',re.IGNORECASE) 
    fastaPatt = re.compile('[.](fasta|fas|fna|faa|fa)$',re.IGNORECASE) 
    gffPatt = re.compile('[.]gff$',re.IGNORECASE) 

    input_files = {}
    for input_file in fileslist:
       
       shortname = None 
       result = None

       result =  gbkPatt.search(input_file)
       if result:
         shortname = re.sub('[.]gbk$','',input_file, re.IGNORECASE) 

       if result==None:
          result =  fastaPatt.search(input_file)
          if result:
             shortname = re.sub('[.](fasta|fas|fna|faa|fa)$','',input_file, re.IGNORECASE) 

       if shortname == None:
          continue

       if re.search('.(fasta|fas|fna|faa|gff|gbk|fa)$',input_file, re.IGNORECASE):
          if check_for_error_in_input_file_name(shortname, globalerrorlogger=globalerrorlogger):
             input_files[input_file] = shortname

    paired_input = {} 
    for key, value in input_files.iteritems():
       paired_input[input_dir + PATHDELIM + key] = path.abspath(output_dir) + PATHDELIM + value

    return paired_input

def removeSuffix(sample_subset_in):
    sample_subset_out = []
    for sample_name in sample_subset_in:
       mod_name = re.sub('.(fasta|fas|fna|faa|gff|gbk|fa)$','',sample_name)
       sample_subset_out.append(mod_name)

    return sample_subset_out


def openGrades():
    pass

def openRank():
    pass

def halt_on_invalid_input(input_output_list, filetypes, sample_subset):

    for samplePath in input_output_list.keys():
       sampleName =  path.basename(input_output_list[samplePath]) 

       ''' in the selected list'''
       if not sampleName in sample_subset:
          continue

       if filetypes[samplePath][0]=='UNKNOWN':
          eprintf("ERROR\tIncorrect input sample %s. Check for bad characters or format\n!", samplePath)
          return False

    return True
          



def report_missing_filenames(input_output_list, sample_subset, logger=None):
    foundFiles = {}
    for samplePath in input_output_list.keys():
       sampleName =  path.basename(input_output_list[samplePath]) 
       foundFiles[sampleName] =True

    for sample_in_subset in sample_subset:
       if not sample_in_subset in foundFiles:
          eprintf("ERROR\tCannot find input file for sample %s\n!", sample_in_subset)
          if logger:
             logger.printf("ERROR\tCannot file input for sample %s!\n", sample_in_subset)

# main function

def sigint_handler(signum, frame):
    eprintf("Received TERMINATION signal\n")
    exit_process()

def main(argv):
    global parser
    (opts, args) = parser.parse_args()
    if valid_arguments(opts, args):
       print usage
       sys.exit(0)

    signal.signal(signal.SIGINT, sigint_handler)
    signal.signal(signal.SIGTERM, sigint_handler)

    eprintf("COMMAND : %s\n", sys.argv[0] + ' ' +  ' '.join(argv))
    # initialize the input directory or file
    input_fp = opts.input_fp 
    output_dir = path.abspath(opts.output_dir)
    verbose = opts.verbose
    print_only = opts.print_only

    sample_subset = removeSuffix(opts.sample_subset)

    run_type = opts.run_type.strip()


    '''no need to remove the whole directory'''
#    if run_type == 'overwrite':
#       force_remove_dir=True
#    else:
#       force_remove_dir=False

    if opts.config_file:
       config_file= opts.config_file
    else:
       config_file = cmd_folder + PATHDELIM + metapaths_config
    
    if opts.ncbi_header and opts.ncbi_sbt:
       if not path.exists(opts.ncbi_header):
          print "Could not open or missing NCBI header file " + opts.ncbi_header
          print "Either disable option to CREATE_SEQUIN_FILE or provide a valid header file"
          sys.exit(0)

       if  not path.exists(opts.ncbi_sbt):
          print """You must must have a sbt file obtained from the NCBI \"Create Submission Template\" form \n 
                 http://www.ncbi.nlm.nih.gov/WebSub/template.cgi """ + opts.ncbi_sbt
          sys.exit(0)

       ncbi_sequin_params = path.abspath(opts.ncbi_header)
       ncbi_sequin_sbt = path.abspath(opts.ncbi_sbt)
    else:
       ncbi_sequin_params = None
       ncbi_sequin_sbt = None

    # try to load the parameter file    
    try:
       if opts.parameter_fp:
          parameter_fp= opts.parameter_fp
       else:
          parameter_fp = cmd_folder + PATHDELIM + metapaths_param
    except IOError:
        raise IOError, ( "Can't open parameters file (%s). Does it exist? Do you have read access?" % opts.parameter_fp )

    
    try:
       if run_type in ['overlay', 'safe'] and not path.exists(output_dir):
             makedirs(output_dir)
    except OSError:
        print ""
        print "ERROR: Cannot create output directory \"" + output_dir + "\"\n"+\
              "       Perhaps directory \"" + output_dir  + "\" already exists.\n" +\
              "       Please choose a different directory, or \n" +\
              "       run with the option \"-r  overwrite\" to force overwrite it."
        sys.exit(1)

        
    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates
    
    command_line_params={}
    command_line_params['verbose']= opts.verbose

    params=parse_metapaths_parameters(parameter_fp)
    """ load the sample inputs  it expects either a fasta 
        file or  a directory containing fasta and yaml file pairs
    """

    globalerrorlogger = WorkflowLogger(generate_log_fp(output_dir, basefile_name= 'global_errors_warnings'), open_mode='w') 

    input_output_list = {}
    if path.isfile(input_fp):   
       """ check if it is a file """
       input_output_list = create_an_input_output_pair(input_fp, output_dir,  globalerrorlogger=globalerrorlogger)
    else:
       if path.exists(input_fp):   
          """ check if dir exists """
          input_output_list = create_input_output_pairs(input_fp, output_dir, globalerrorlogger=globalerrorlogger)
       else:   
          """ must be an error """
          eprintf("ERROR\tNo valid input sample file or directory containing samples exists .!")
          eprintf("ERROR\tAs provided as arguments in the -in option.!\n")
          exit_process("ERROR\tAs provided as arguments in the -in option.!\n")
   
    """ these are the subset of sample to process if specified
        in case of an empty subset process all the sample """

    # remove all samples that are not specifed unless sample_subset is empty
    remove_unspecified_samples(input_output_list, sample_subset, globalerrorlogger = globalerrorlogger)

    # add check the config parameters 
    sorted_input_output_list = sorted(input_output_list.keys())

    filetypes = check_file_types(sorted_input_output_list) 

    #stop on in valid samples
    if not halt_on_invalid_input(input_output_list, filetypes, sample_subset):
       globalerrorlogger.printf("ERROR\tInvalid inputs found. Check for file with bad format or characters!\n")
       halt_process(opts.delay)

    # make sure the sample files are found
    report_missing_filenames(input_output_list, sample_subset, logger=globalerrorlogger)


    #check the pipeline configuration
    config_settings = read_pipeline_configuration(config_file, globalerrorlogger)

    parameter =  Parameters()
    if not staticDiagnose(config_settings, params, logger = globalerrorlogger):
        eprintf("ERROR\tFailed to pass the test for required scripts and inputs before run\n")
        globalerrorlogger.printf("ERROR\tFailed to pass the test for required scripts and inputs before run\n")
        halt_process(opts.delay)

    
    samplesData = {}
    # PART1 before the blast

    block_mode = opts.block_mode
    runid = opts.runid

    try:
         # load the sample information 
         print "RUNNING MetaPathways version 2.5.1"
         if len(input_output_list): 
              for input_file in sorted_input_output_list:
                sample_output_dir = input_output_list[input_file]
                algorithm = get_parameter(params, 'annotation', 'algorithm', default='LAST').upper()
   
                s = SampleData() 
                s.setInputOutput(inputFile = input_file, sample_output_dir = sample_output_dir)
                s.setParameter('algorithm', algorithm)
                s.setParameter('ncbi_params_file', ncbi_sequin_params)
                s.setParameter('ncbi_sequin_sbt', ncbi_sequin_sbt)
                s.setParameter('FILE_TYPE', filetypes[input_file][0])
                if params["INPUT"]['format'] in ["gbk-annotated", "gff-annotated"]:
                    s.setParameter('ANNOTATED', True)
                s.setParameter('SEQ_TYPE', filetypes[input_file][1])
                s.clearJobs()
   
                if run_type=='overwrite' and  path.exists(sample_output_dir):
                   shutil.rmtree(sample_output_dir)
                   makedirs(sample_output_dir)
                if not  path.exists(sample_output_dir):
                   makedirs(sample_output_dir)
   
                s.prepareToRun()
                samplesData[input_file] = s
   
              # load the sample information 
              run_metapathways(
                   samplesData,
                   sample_output_dir,
                   output_dir,
                   globallogger = globalerrorlogger,
                   command_line_params=command_line_params,
                   params=params,
                   metapaths_config=metapaths_config,
                   status_update_callback=status_update_callback,
                   config_file=config_file,
                   run_type = run_type, 
                   config_settings = config_settings,
                   block_mode = block_mode,
                   runid = runid
              )
         else: 
              eprintf("ERROR\tNo valid input files/Or no files specified  to process in folder %s!\n",sQuote(input_fp) )
              globalerrorlogger.printf("ERROR\tNo valid input files to process in folder %s!\n",sQuote(input_fp) )
   
        
         # blast the files
     
         blasting_system =    get_parameter(params,  'metapaths_steps', 'BLAST_REFDB', default='yes')
         if blasting_system =='grid':
            #  blasting the files files on the grids
             input_files = sorted_input_output_list
             blast_in_grid(
                   sampleData[input_file],
                   input_files, 
                   path.abspath(opts.output_dir),   #important to use opts.
                   params=params,
                   metapaths_config=metapaths_config,
                   config_file=config_file,
                   run_type = run_type,
                   runid = runid
                )
     
    except:
       exit_process(str(traceback.format_exc(10)), logger= globalerrorlogger )


    
    eprintf("            ***********                \n")
    eprintf("INFO : FINISHED PROCESSING THE SAMPLES \n")
    eprintf("             THE END                   \n")
    eprintf("            ***********                \n")
    halt_process(opts.delay)

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])    
    

