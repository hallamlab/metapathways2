#!/usr/bin/python

"""This script run the orf prediction """

try:
   import  sys, re, csv, traceback
   from os import path, _exit
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


usage = sys.argv[0] + """ --algorithm <algorithm> [algorithm dependent options]"""

parser = None
def createParser():
    global parser

    epilog = """The preprocessed nucleotide sequences (contigs) are used as input to a gene prediction algorithm, currently prodigal, to detect the gene coding regions.  The output of the prodigal run is a set of untranslated ORFs and the same ORFs translated (into amino acid sequences). The resulting files are available in the 'orf_prediction' folder. The translation is done based on the translation table id provided by the user, by default it 11"""

    epilog = re.sub(r'\s+',' ', epilog)
 
    parser = OptionParser(usage=usage, epilog=epilog)

    # Input options

    parser.add_option('--algorithm', dest='algorithm', default='prodigal',
                           help='ORF prediction algorithm')

    prodigal_group =  OptionGroup(parser, 'Prodigal parameters')

    prodigal_group.add_option('--prod_input', dest='prod_input', default=None,
                           help='the input sequences  <inputfile>')

    prodigal_group.add_option('--prod_output', dest='prod_output', default=None,
                           help='the output <outfile>')

    prodigal_group.add_option('--prod_p', dest='prod_p', default=None,
                           help='Select procedure (single or meta).  Default is single')

    prodigal_group.add_option('--prod_f', dest='prod_f', default='gff',
                           help='Select output format (gbk, gff, or sco).  default is gff')

    prodigal_group.add_option('--prod_g', dest='prod_g', default='11',
                           help='Specify a translation table to use (default 11)')

    prodigal_group.add_option('--prod_m', dest='prod_m', action='store_true', default=False,
                           help='Treat runs of n\'s as masked sequence and do not build genes across them')

    prodigal_group.add_option('--prod_exec', dest='prod_exec', default=None,
                           help='prodigal executable')

    parser.add_option_group(prodigal_group)



def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser

    options, args = parser.parse_args(argv)

    if options.algorithm == 'prodigal':
        _execute_prodigal(options)
       


def  _execute_prodigal(options):
    args= [ ]

    if options.prod_exec :
       args.append( options.prod_exec )

    if options.prod_m:
       args.append("-m")

    if options.prod_p:
       args += [ "-p", options.prod_p ]
    
    if options.prod_f:
       args += [ "-f", options.prod_f ]
    
    if options.prod_g:
       args += [ "-g", options.prod_g ]
    
    if options.prod_input:
       args += [ "-i", options.prod_input ]

    if options.prod_output:
       args += [ "-o", options.prod_output ]

    result = getstatusoutput(' '.join(args) )
    return result[0]



def MetaPathways_orf_prediction(argv, extra_command = None, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\tORF_PREDICTION\n")
    createParser()
    main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (0,'')

if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

