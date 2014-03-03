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
     import sys, traceback
     from os import makedirs, sys, listdir, environ, path
     import re 
     import inspect
     #from commands import getstatusoutput
     from optparse import OptionParser
     import shutil 
     
     from libs.python_modules import metapaths_utils
     from libs.python_modules.metapaths_utils  import parse_command_line_parameters
     from libs.python_modules.parse  import parse_metapaths_parameters, parse_parameter_file
     from libs.python_modules.metapaths_pipeline import print_commands, call_commands_serially, print_to_stdout, no_status_updates
     from libs.python_modules.sysutil import pathDelim
     
     from libs.python_modules.metapaths import run_metapathways_before_BLAST, run_metapathways_at_BLAST, run_metapathways_after_BLAST, get_parameter
     from libs.python_modules.annotate import *
     from libs.python_modules.blast_using_grid import blast_in_grid
except:
   print """ Could not load some user defined  module functions"""
   print """ Make sure your typed \"source MetaPathwaysrc\""""
   print """ """
   print traceback.print_exc(10)
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

usage= """./MetaPathways.py  -i input_file -o outdir  -p parameters.txt 
For more options:  ./MetaPathways.py -h"""
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
                         '\n(c)\'dry-run\' -- shows the steps that are going to be computed or not\n' +
                         '\n(d)\'safe\' -- safe mode does not run on an existing run folder\n')
#ith out of order completion \ time-stamps in the \'workflow_log.txt\' 
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="print lots of information on the stdout [default]")
parser.add_option("-P", "--print-only",
                  action="store_true", dest="print_only", default=False,
                  help="print only  the commands [default False]")

parser.add_option("-n", "--ncbi-header", dest="ncbi_header", 
                  help="NCBI sequin submission parameters file" )

parser.add_option("-s", "--ncbi-sbt-file", dest="ncbi_sbt", 
                  help="the NCBI sbt location created by the \"Create Submission Template\" form: http://www.ncbi.nlm.nih.gov/WebSub/template.cgi" )



# checks if the supplied arguments are adequate
def valid_arguments(opts, args):
    if (opts.input_fp == None and opts.output_dir ==None )  or\
     opts.output_dir == None or opts.parameter_fp == None :
       return True
    else:
       return False


#creates an input output pair if input is just an input file
def create_an_input_output_pair(input_file, output_dir):
    input_output = {}
    shortname = re.sub('[.](fasta|fas|fna|faa|gbk|gff|fa)$','',input_file, re.I) 
    shortname = re.sub(r'.*' + PATHDELIM ,'',shortname) 
    shortname = re.sub(r'[.]','_',shortname) 
    if re.search(r'.(fasta|fas|fna|faa|gbk|gff|fa)$',input_file, re.I):
       input_output[input_file] = path.abspath(output_dir) + PATHDELIM + shortname

    return input_output


#creates a list of  input output pairs if input is  an input dir
def create_input_output_pairs(input_dir, output_dir):
    fileslist =  listdir(input_dir)

    input_files = {}
    for file in fileslist:
       shortname = re.sub('.(fasta|fas|fna|faa|gbk|gff|fa)$','',file,re.I) 
       shortname = re.sub(r'[.]','_',shortname) 
       if re.search('.(fasta|fas|fna|faa|gff|gbk|fa)$',file, re.I):
         input_files[file] = shortname

    paired_input = {} 
    for key, value in input_files.iteritems():
            paired_input[input_dir + PATHDELIM + key] = path.abspath(output_dir) + PATHDELIM + value

    return paired_input

def openGrades():
    pass

def openRank():
    pass

# main function


def main(argv):

    (opts, args) = parser.parse_args()
    if valid_arguments(opts, args):
       print usage
       sys.exit(0)

    # initialize the input directory or file
    input_fp = opts.input_fp 
    output_dir = path.abspath(opts.output_dir)
    verbose = opts.verbose
    print_only = opts.print_only

    run_type = opts.run_type.strip()

    if run_type == 'overwrite':
       force_remove_dir=True
    else:
       force_remove_dir=False

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
        parameter_f = open(opts.parameter_fp)
    except IOError:
        raise IOError,\
         "Can't open parameters file (%s). Does it exist? Do you have read access?"\
         % opts.parameter_fp

    
    if force_remove_dir:
        try: 
           if path.exists(output_dir):
              shutil.rmtree(output_dir)
        except OSError:
           print "ERROR: Cannot remove directory: " + output_dir
           sys.exit(1)

    try:
       if (run_type in ['overlay', 'safe'] or force_remove_dir) and not path.exists(output_dir):
             makedirs(output_dir)
       elif run_type in ['safe'] and path.exists(output_dir):
        print ""
        print "ERROR: Cannot create output directory \"" + output_dir + "\"\n"+\
              "       Perhaps output directory already exists.\n" +\
              "       Please choose a different directory, or \n" +\
              "       run with the option \"-r  overwrite\" to force overwrite it."
        sys.exit(1)
    except OSError:
        print ""
        print "ERROR: Cannot create output directory \"" + output_dir + "\"\n"+\
              "       Perhaps directory \"" + output_dir  + "\" already exists.\n" +\
              "       Please choose a different directory, or \n" +\
              "       run with the option \"-r  overwrite\" to force overwrite it."
        sys.exit(1)
        
    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially
    
    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates
    
    command_line_params={}
    command_line_params['verbose']= opts.verbose


    # load the sample inputs  it expects either a fasta file or  a directory 
    # containing fasta and yaml file pairs

    input_output_list = {}
    if path.isfile(input_fp):   # check if it is a file
       input_output_list = create_an_input_output_pair(input_fp, output_dir)
    else:
       if path.exists(input_fp):   # check if dir exists
          input_output_list = create_input_output_pairs(input_fp, output_dir)
       else:   # must be an error
          print "No valid input sample file or directory containing samples exists .!"
          print "As provided as arguments in the -in option.!"
          sys.exit(1)
   
    #print input_output_list
    #sys.exit(1)

    config_params=parse_metapaths_parameters(parameter_f)

    # add check the config parameters 
    
    sorted_input_output_list = sorted(input_output_list.keys())
    # PART1 before the blast
    if len(input_output_list): 
      for input_file in sorted_input_output_list:
        output_dir = input_output_list[input_file]
        if run_type=='overwrite' and  path.exists(output_dir):
           shutil.rmtree(output_dir)
           makedirs(output_dir)
        if not  path.exists(output_dir):
           makedirs(output_dir)
        print ""
        print "##################################################"
        print "PROCESSING INPUT " + input_file
        run_metapathways_before_BLAST(
           input_file, 
           output_dir,
           command_handler=command_handler,
           command_line_params=command_line_params,
           config_params=config_params,
           metapaths_config=metapaths_config,
           status_update_callback=status_update_callback,
           config_file=config_file,
           ncbi_sequin_params = ncbi_sequin_params,
           ncbi_sequin_sbt = ncbi_sequin_sbt,
           run_type = run_type
        )
    else: 
        print "ERROR :No input files were found in the input folder"
        sys.exit(0)

    # blast the files

    blasting_system =    get_parameter(config_params,  'metapaths_steps', 'BLAST_REFDB', default='yes')

    if blasting_system =='grid':
       #  blasting the files files on the grids
        input_files = sorted_input_output_list
        blast_in_grid(
              input_files, 
              path.abspath(opts.output_dir),   #important to use opts.
              config_params=config_params,
              metapaths_config=metapaths_config,
              config_file=config_file,
              run_type = run_type
           )

    else:
       #  blasting  the files files locally
       for input_file in sorted_input_output_list:
           output_dir = input_output_list[input_file]
           run_metapathways_at_BLAST(
              input_file, 
              output_dir,
              command_handler=command_handler,
              command_line_params=command_line_params,
              config_params=config_params,
              metapaths_config=metapaths_config,
              status_update_callback=status_update_callback,
              config_file=config_file,
              ncbi_sequin_params = ncbi_sequin_params,
              ncbi_sequin_sbt = ncbi_sequin_sbt,
              run_type = run_type
           )

    # after blasting  the files
    for input_file in sorted_input_output_list:
        output_dir = input_output_list[input_file]
        run_metapathways_after_BLAST(
           input_file, 
           output_dir,
           command_handler=command_handler,
           command_line_params=command_line_params,
           config_params=config_params,
           metapaths_config=metapaths_config,
           status_update_callback=status_update_callback,
           config_file=config_file,
           ncbi_sequin_params = ncbi_sequin_params,
           ncbi_sequin_sbt = ncbi_sequin_sbt,
           run_type = run_type
        )

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])    

