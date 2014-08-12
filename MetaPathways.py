from __future__ import division

__author__ = "Kishori M Konwar Niels W Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
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
     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, eprintf, halt_process, exit_process,WorkflowLogger, generate_log_fp
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
metapaths_config = """template_config.txt""";

script_info={}
script_info['brief_description'] = """A workflow script for making PGDBs from metagenomic sequences"""
script_info['script_description'] = """ takes a sequence file and performs all processing steps through building the OTU table.
             REQUIRED: You must have a fas and an yaml file  and  a custom parameters file:"""
script_info['script_usage'] = []


usage=  sys.argv[0] + """  -i input_file -o outdir  -p parameters.txt 
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
                      action="store_true", dest="block_mode", default=False,
                      help="processes the samples by blocking the stages before and after functional search [default off]")
    
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
     opts.output_dir == None or opts.parameter_fp == None :
       return True
    else:
       return False

def remove_unspecified_samples(input_output_list, sample_subset, format, globalerrorlogger = None):
   """ keep only the samples that are specified  before processing  """
   shortened_names = {}

   for input_file in input_output_list.keys():
      shortname = None 
      if format in ['gbk-unannotated', 'gbk-annotated']:
          shortname = re.sub('[.]gbk$','',input_file, re.IGNORECASE) 
      elif format =='fasta':
          shortname = re.sub('[.](fasta|fas|fna|faa|fa)$','',input_file, re.IGNORECASE) 

      if shortname==None:
         continue

      shortname = re.sub(r'[.]','_',shortname) 
      shortened_names[shortname] = input_file

   shortened_subset_names = [] 
   for sample_in_subset in sample_subset:
      shortname = None 
      if format in ['gbk-unannotated', 'gbk-annotated']:
          shortname = re.sub('[.]gbk$','',sample_in_subset, re.IGNORECASE) 
      elif format =='fasta':
          shortname = re.sub('[.](fasta|fas|fna|faa|fa)$','',sample_in_subset, re.IGNORECASE) 

      if shortname==None:
         continue

      if check_for_error_in_input_file_name(shortname, globalerrorlogger=globalerrorlogger):
         shortened_subset_names.append(shortname)

   samples_to_keep = {} 

   for keep_sample in shortened_subset_names:
      sampleMatchPAT = re.compile(r'' + keep_sample + '$') 
      for sample  in shortened_names:
         result = sampleMatchPAT.search(sample, re.IGNORECASE)
         if result:
            samples_to_keep[shortened_names[sample]]= True
            break
    
   input_sample_list = input_output_list.keys()
   for sample in input_sample_list:
      if not sample in samples_to_keep:
         del input_output_list[sample]


def check_for_error_in_input_file_name(shortname, globalerrorlogger=None):

    """  creates a list of  input output pairs if input is  an input dir """
    clean = True
    if not re.search(r'^[a-zA-Z]',shortname):
         eprintf("ERROR\tSample name %s must begin with an alphabet!\n",shortname)
         if globalerrorlogger:
            globalerrorlogger.printf("ERROR\tSample name %s must begin with an alphabet!\n",shortname)
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

    errmessage = """ Sample names before the  suffixes .fasta, .fas, .fna, .faa or .gbk, must  consist only of alphabets, digits and _; and should consist of at least two characters """
    eprintf("ERROR\t%s\n",errmessage)
    if globalerrorlogger:
        globalerrorlogger.printf("ERROR\t%s\n",errmessage)
        exit_process("ERROR\t" + errmessage + "Exiting!" + "\n")
    return False


def create_an_input_output_pair(input_file, output_dir, format, globalerrorlogger=None):
    """ creates an input output pair if input is just an input file """
       
    input_output = {}

    if not re.search(r'.(fasta|fas|fna|faa|gbk|gff|fa)$',input_file, re.IGNORECASE):
       return input_output

    shortname = None 
    if format in ['gbk-unannotated', 'gbk-annotated']:
        shortname = re.sub('[.]gbk$','',input_file, re.IGNORECASE) 
    elif format =='fasta':
        shortname = re.sub('[.](fasta|fas|fna|faa|fa)$','',input_file, re.IGNORECASE) 
    else:
        shortname = re.sub('[.]gff$','',input_file, re.IGNORECASE) 

    shortname = re.sub(r'.*' + PATHDELIM ,'',shortname) 

    check_for_error_in_input_file_name(shortname, globalerrorlogger=globalerrorlogger)

    input_output[input_file] = path.abspath(output_dir) + PATHDELIM + shortname

    return input_output


def create_input_output_pairs(input_dir, output_dir, format, globalerrorlogger=None):
    """  creates a list of  input output pairs if input is  an input dir """
    fileslist =  listdir(input_dir)
    gbkPatt = re.compile('[.]gbk$',re.IGNORECASE) 

    fastaPatt = re.compile('[.](fasta|fas|fna|faa|fa)$',re.IGNORECASE) 
    gffPatt = re.compile('[.]gff$',re.IGNORECASE) 

    input_files = {}
    for input_file in fileslist:
       shortname = None 
       if format in ['gbk-unannotated', 'gbk-annotated']:
          result =  gbkPatt.search(input_file)
          if result:
             shortname = re.sub('[.]gbk$','',input_file, re.IGNORECASE) 

       elif format in [ 'fasta' ]:
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

def openGrades():
    pass

def openRank():
    pass

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
    sample_subset= opts.sample_subset

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
        parameter_f = opts.parameter_fp
    except IOError:
        raise IOError,\
         "Can't open parameters file (%s). Does it exist? Do you have read access?"\
         % opts.parameter_fp

    
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

    params=parse_metapaths_parameters(parameter_f)
    format = params['INPUT']['format']

    """ load the sample inputs  it expects either a fasta 
        file or  a directory containing fasta and yaml file pairs
    """

    globalerrorlogger = WorkflowLogger(generate_log_fp(output_dir, basefile_name= 'global_errors_warnings'), open_mode='w') 

    input_output_list = {}
    if path.isfile(input_fp):   
       """ check if it is a file """
       input_output_list = create_an_input_output_pair(input_fp, output_dir, format, globalerrorlogger = globalerrorlogger)
    else:
       if path.exists(input_fp):   
          """ check if dir exists """
          input_output_list = create_input_output_pairs(input_fp, output_dir, format, globalerrorlogger=globalerrorlogger)
       else:   
          """ must be an error """
          eprintf("ERROR\tNo valid input sample file or directory containing samples exists .!")
          eprintf("ERROR\tAs provided as arguments in the -in option.!\n")
          exit_process("ERROR\tAs provided as arguments in the -in option.!\n")
   
    """ these are the subset of sample to process if specified
        in case of an empty subset process all the sample """
    if sample_subset:
       remove_unspecified_samples(input_output_list, sample_subset, format, globalerrorlogger = globalerrorlogger)


    # add check the config parameters 
    sorted_input_output_list = sorted(input_output_list.keys())

    config_settings = read_pipeline_configuration(config_file, globalerrorlogger)

    parameter =  Parameters()
    if not staticDiagnose(config_settings, params, logger = globalerrorlogger):
        eprintf("ERROR\tFailed to pass the test for required scripts and inputs before run\n")
        globalerrorlogger.printf("ERROR\tFailed to pass the test for required scripts and inputs before run\n")
        exit_process("ERROR\tFailed to pass the test for required scripts and inputs before run\n")

    
    
    samplesData = {}
    # PART1 before the blast

    block_mode = opts.block_mode
    runid = opts.runid

    try:
         # load the sample information 
         if len(input_output_list): 
              for input_file in sorted_input_output_list:
                sample_output_dir = input_output_list[input_file]
                algorithm = get_parameter(params, 'annotation', 'algorithm', default='LAST').upper()
   
                s = SampleData() 
                s.setInputOutput(inputFile = input_file, sample_output_dir = sample_output_dir)
                s.setParameter('algorithm', algorithm)
                s.setParameter('ncbi_params_file', ncbi_sequin_params)
                s.setParameter('ncbi_sequin_sbt', ncbi_sequin_sbt)
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
              eprintf("ERROR\tNo input files in the specified folder %s to process!\n",sQuote(input_fp) )
              globalerrorlogger.printf("ERROR\tNo input files in the specified folder %s to process!\n",sQuote(input_fp) )
   
        
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
       globalerrorlogger.write( "ERROR\t" + str(traceback.format_exc(10)))
       exit_process("ERROR:" + str(traceback.format_exc(10)))


    
    eprintf("            ***********                \n")
    eprintf("INFO : FINISHED PROCESSING THE SAMPLES \n")
    eprintf("             THE END                   \n")
    eprintf("            ***********                \n")
    halt_process(4)

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])    
    

