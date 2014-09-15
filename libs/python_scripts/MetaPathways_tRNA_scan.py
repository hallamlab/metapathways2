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


help = sys.argv[0] + """ -i input -o output [algorithm dependent options]"""

parser = None
def createParser():
    global parser
    epilog = """This script is used for scanning for tRNA,  using tRNA-Scan 1.4, 
              on the set of metagenomics sample sequences """
    epilog = re.sub(r'\s+',' ', epilog)

    parser = OptionParser(usage=help, epilog=epilog)

    # Input options


    parser.add_option('-o', dest='trna_o', default=None,
                           help='Output from the tRNA-Scan 1.4 into <outfile>')

    parser.add_option('-i', dest='trna_i', default=None,
                           help='reads the sequences from <input> file')

    parser.add_option('-T', dest='trna_T', default='6',
                           help='reads the Tsignal from <TPCsignal>')

    parser.add_option('-D', dest='trna_D', default=None,
                           help='reads the Dsignal from <Dsignal>')

    parser.add_option('-F', dest='trna_F', default=None,
                           help='write predicted tRNA genes in fasta format<outfile>')

    parser.add_option('--executable', dest='trna_executable', default=None,
                           help='The tRNA-SCAN 1.4 executable')



def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser
    options, args = parser.parse_args(argv)
    return _execute_tRNA_Scan(options)


def  _execute_tRNA_Scan(options):
    args= [ ]

    if options.trna_executable :
       args.append( options.trna_executable )

    if options.trna_i:
       args += [ "-i", options.trna_i ]
    
    if options.trna_o:
       args += [ "-o", options.trna_o ]

    if options.trna_D:
       args += [ "-D",  options.trna_D ]
    
    if options.trna_T:
       args += [ "-T",  options.trna_T ]
    
    if options.trna_F:
       args += [ "-F",  options.trna_F]
    result = getstatusoutput(' '.join(args) )

    return result
    

def MetaPathways_tRNA_scan(argv, extra_command = None, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\ttRNA_SCAN\n")
    createParser()
    result = main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (result[0],'')

if __name__ == '__main__':
    createParser()
    result = main(sys.argv[1:])

