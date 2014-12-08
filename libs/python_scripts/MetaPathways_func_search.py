#!/usr/bin/python

"""This script run the orf prediction """

try:
   import  sys, re, csv, traceback
   from os import path, _exit, rename
   import logging.handlers
   from optparse import OptionParser, OptionGroup
   from libs.python_modules.utils.sysutil import pathDelim
   from libs.python_modules.utils.metapathways_utils  import fprintf, printf, eprintf,  exit_process
   from libs.python_modules.utils.sysutil import getstatusoutput

   from libs.python_modules.utils.pathwaytoolsutils import *

except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     print traceback.print_exc(10)
     sys.exit(3)

PATHDELIM=pathDelim()



def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

def files_exist( files , errorlogger = None):
    status = True    
    for file in files:
       if not path.exists(file):
          if errorlogger:
             errorlogger.write( 'ERROR\tCould not find ptools input  file : ' +  file )
          status = False
    return not status


usage = sys.argv[0] + """ -i input -o output [algorithm dependent options]"""

parser = None
def createParser():
    global parser
    epilog = """This script is used for running a homology search algorithm such as BLAST or LAST
              on a set of query amino acid sequences against a target of  reference protein sequences.
              Currently it supports the BLASTP and LAST algorithm. Any other homology search algorithm 
              can be added by first adding the new algorithm name in upper caseusing in to the 
              choices parameter in the algorithm option of this script.
              The results are put in a tabular form in the folder blast_results, with individual files 
               for each of the databases. The files are named as "<samplename>.<dbname>.<algorithm>out"
              In the case of large number of amino acid sequences, this step of the computation can be 
               also done using multiple grids (to use batch processing system) """

    epilog = re.sub(r'\s+',' ', epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    # Input options

    parser.add_option('--algorithm', dest='algorithm', default='BLAST', choices = ['BLAST', 'LAST'],
                           help='the homology search algorithm')

    blast_group =  OptionGroup(parser, 'BLAST parameters')

    blast_group.add_option('--blast_query', dest='blast_query', default=None,
                           help='Query amino acid sequences for BLASTP')

    blast_group.add_option('--blast_db', dest='blast_db', default=None,
                           help='Target reference database sequenes for BLASTP')

    blast_group.add_option('--blast_out', dest='blast_out', default=None,
                           help='BLAST output file')

    blast_group.add_option('--blast_outfmt', dest='blast_outfmt', default='6',
                           help='BLASTP output format [default 6, tabular]')

    blast_group.add_option('--blast_evalue', dest='blast_evalue', default=None,
                           help='The e-value cutoff for the BLASTP')

    blast_group.add_option('--blast_num_threads', dest='blast_num_threads', default=None,
                           help='Number of BLAST threads')

    blast_group.add_option('--blast_max_target_seqs', dest='blast_max_target_seqs', default=None,
                           help='Maximum number of target hits per query')

    blast_group.add_option('--blast_executable', dest='blast_executable',  default=None,
                           help='The BLASTP executable')

    parser.add_option_group(blast_group)

    last_group =  OptionGroup(parser, 'LAST parameters')

    last_group.add_option('--last_query', dest='last_query', default=None,
                           help='Query amino acid sequences for LAST')

    last_group.add_option('--last_db', dest='last_db', default=None,
                           help='Target reference database sequenes for LAST')

    last_group.add_option('--last_f', dest='last_f', default='0',
                           help='LAST output format [default 0, tabular]')

    last_group.add_option('--last_o', dest='last_o', default=None,
                           help='LAST output file')

    last_group.add_option('--last_executable', dest='last_executable',  default=None,
                           help='The LAST executable')

    parser.add_option_group(last_group)


def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser

    options, args = parser.parse_args(argv)

    if options.algorithm == 'BLAST':
       _execute_BLAST(options, logger = errorlogger)
    elif options.algorithm == 'LAST':
        _execute_LAST(options, logger = errorlogger)
    else:
        eprintf("ERROR\tUnrecognized algorithm name for FUNC_SEARCH\n")
        if errorlogger:
            errorlogger.printf("ERROR\tUnrecognized algorithm name for FUNC_SEARCH\n")
        exit_process("ERROR\tUnrecognized algorithm name for FUNC_SEARCH\n")


def  _execute_LAST(options, logger = None):
    args= [ ]

    if options.last_executable :
       args.append( options.last_executable )

    if options.last_f:
       args += [ "-f", options.last_f ]
    
    if options.last_o:
       args += [ "-o", options.last_o + ".tmp"]

    if options.last_db:
       args += [ options.last_db ]
    
    if options.last_query:
       args += [ options.last_query ]

    try:
       result = getstatusoutput(' '.join(args) )
       rename(options.last_o + ".tmp", options.last_o) 
    except:
       message = "Could not run lastal correctly"
       if result and len(result) > 1:
          message = result[1]
       if logger:
          logger.printf("ERROR\t%s\n", message)
       return '1'

    return result[0]
    

def  _execute_BLAST(options, logger = None):
    args= [ ]

    if options.blast_executable :
       args.append( options.blast_executable )

    if options.blast_max_target_seqs:
       args +=["-max_target_seqs", options.blast_max_target_seqs]

    if options.blast_num_threads:
       args += [ "-num_threads", options.blast_num_threads ]
    
    if options.blast_outfmt:
       args += [ "-outfmt", options.blast_outfmt ]
    
    if options.blast_db:
       args += [ "-db", options.blast_db ]
    
    if options.blast_query:
       args += [ "-query", options.blast_query ]
    
    if options.blast_evalue:
       args += [ "-evalue", options.blast_evalue ]

    if options.blast_out:
       args += [ "-out", options.blast_out + ".tmp" ]

    try:
       result = getstatusoutput(' '.join(args) )
       rename(options.blast_out + ".tmp", options.blast_out) 
    except:
       return '1'

    return result[0]



def MetaPathways_func_search(argv, extra_command = None, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\tFUNC_SEARCH\n")
    createParser()
    main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (0,'')

if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

