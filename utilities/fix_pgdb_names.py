#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar, Niels W. Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     import re
     import sys
     from optparse import OptionParser, OptionGroup
     from glob import glob
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwayrc\""""
     print """ """
     sys.exit(3)

script_name = sys.argv[0]
usage= script_name + """ --pgdb-folder --specific-pgdb [ all if not specified ]"""
parser = OptionParser(usage)

parser.add_option( "--pgdb-folder", dest="pgdb_folder",  
                  help='the pgdb folder')

parser.add_option( "--pgdb", dest="pgdb", default = [], 
                  help='a specific pgdb name [default: selects all in the --pgdb-folder]')

def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)



def check_arguments(opts, args):
    if opts.pgdb_folder == None:
         return False
    return True

def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

     return attributes


def  fix_pgdb_input_files(pgdb_folder, pgdbs = []):
     pgdb_list = glob(pgdb_folder + '*cyc/*/input/organism.dat')     

     for pgdb_organism_file in pgdb_list:
        process_organism_file(pgdb_organism_file)


def fixLine(line, id):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[0]+'\t' + id
     

def getID(line):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[1]
     
def process_organism_file(filel):
     patternsToFix = [ re.compile(r'NAME\tunclassified sequences'), re.compile(r'ABBREV-NAME\tu. sequences') ]
     patternID =  re.compile(r'^ID\t.*')
     try:
         orgfile = open(filel,'r')
     except IOError:
         print "ERROR : Cannot open organism file" + str(filel)
         return 

     lines = orgfile.readlines()
     newlines = []

     needsFixing = False

     id = None
     for line in lines:
         line = line.strip()
         if len(line)==0:
            continue
         flag = False

         result = patternID.search(line)
         if result:   
             id = getID(line)
          
         for patternToFix in patternsToFix:
             result = patternToFix.search(line)
             if result and id:
                 newline = fixLine(line, id)
                 newlines.append(newline)
                 flag= True
                 needsFixing = True

         if flag==False:
            newlines.append(line)

     orgfile.close()
     if needsFixing:
       write_new_file(newlines, filel)


def write_new_file(lines, output_file):
    
    print "Fixing file " + output_file 
    try:
       outputfile = open(output_file,'w')
       pass
    except IOError:
         print "ERROR :Cannot open output file "  + output_file
   
    for line in lines:
       fprintf(outputfile, "%s\n", line)

    outputfile.close()

# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)

    pgdb_folder = opts.pgdb_folder
    pgdbs = opts.pgdb

    fix_pgdb_input_files(pgdb_folder, pgdbs = pgdbs)

    
    

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

