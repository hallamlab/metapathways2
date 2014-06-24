#!/usr/bin/python
# File created on 27 Jan 2012.
from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     from os import makedirs, sys, remove
     from sys import path
     import re
     
     from optparse import OptionParser, OptionGroup
     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, fprintf, printf
     from libs.python_modules.utils.sysutil import getstatusoutput
     from libs.python_modules.annotate.sequence import genbank, fasta
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)


usage= """./MetapathWays_copy_gff.py --source  source --target target [--source  source --target target ]"""
parser = OptionParser(usage)

parser.add_option("-s", "--source", dest="source", default=[], action='append', 
                  help='the source file')

parser.add_option("-t", "--target", dest="target", default=[], action='append', 
                  help='the target file')



def check_arguments(opts, args):
    if len(opts.source) ==0  or len(opts.source)%2 ==  1:
        return False
      
    if len(opts.source) != len(opts.target):
         print "The nuber of sources and targets should be the same"
         return False
    return True



# the main function
def copy_faa_gff_orf_prediction( source_files, target_files) :
      for source, target in zip(source_files, target_files):
         #print source + ' ' + target
         sourcefile = open(source, 'r')
         targetfile = open(target, 'w')
         sourcelines = sourcefile.readlines()
         for line in sourcelines:
            fprintf(targetfile, "%s\n", line.strip())
 
         sourcefile.close()
         targetfile.close()

def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)
    copy_faa_gff_orf_prediction( opts.source, opts.target) 



# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

