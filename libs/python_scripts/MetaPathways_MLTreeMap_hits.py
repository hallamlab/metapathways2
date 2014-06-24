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
     import os
     import re
     from os import makedirs, sys, remove, listdir
     from sys import path

     from optparse import OptionParser
     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, fprintf
     from libs.python_modules.utils.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)






usage= """./MetaPathway_MLTreeMap_hits.py -i input_folder -o output_file """
parser = OptionParser(usage)
parser.add_option("-i", "--input-folder", dest="input_folder",
                  help='the input mltreemap output folder [REQUIRED]')
parser.add_option("-o", "--output-file", dest="output_file",
                  help='the output file for COG hits [REQUIRED]')


def check_arguments(opts, args):
    if opts.input_folder == None or opts.output_file == None:
       return True
    else:
       return False


def main(argv): 
    (opts, args) = parser.parse_args()
    if check_arguments(opts, args):
       print usage
       sys.exit(0)

    input_folder = opts.input_folder
    output_file = opts.output_file

    filePATTERN = re.compile(r'.*COG[0-9]*.*\.fa');
    cogSeqMatchesPATTERN = re.compile(r'[a-zA-Z]*_(.*)__[0-9]*__*(COG[0-9]*).*.fa');
    list= []
    for file in  listdir(input_folder):
      if filePATTERN.match(file):
         hits =  cogSeqMatchesPATTERN.search( file) 
         if hits:
             list.append( (hits.group(1), hits.group(2)) )
         

    try:
        outputfile  = open(output_file, 'w')
    except:
        print "Cannot open file to MLTreeMap hits"
        sys.exit(0)




    fprintf(outputfile, "Sequences\tCOG\n")
    for seq, cog in list:
        fprintf(outputfile, "%s\t%s\n",seq, cog)

    outputfile.close()

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

