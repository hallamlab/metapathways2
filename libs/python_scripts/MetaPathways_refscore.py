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
     import os, re
     from os import makedirs, sys, remove
     from sys import path
     from optparse import OptionParser

     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, fprintf
     from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)

PATHDELIM = pathDelim()

usage= """./run_pgdb_pipeline.py -i input_fasta_file -o output_file [-B blast_executable -F formatdb_executable] """
parser = None
def createParser():
    global parser
    parser = OptionParser(usage)
    parser.add_option("-i", "--input_file", dest="input_fasta",
                      help='the input fasta file [REQUIRED]')
    parser.add_option("-o", "--output_file", dest="output_file",
                      help='the output fasta file [REQUIRED]')
    parser.add_option("-B", "--BLAST_EXEUTABLE", dest="blast_executable",
                      help='the BLAST executable  [REQUIRED]')
    parser.add_option("-F", "--FORMAT_EXECUTABLE", dest="formatdb_executable",
                      help='the FORMATDB executable file [REQUIRED]')
    parser.add_option("-a", "--algorithm", dest="algorithm", choices = ['BLAST', 'LAST'], default = "BLAST", 
                      help='the algorithm used for computing homology [DEFAULT: BLAST]')


def check_arguments(opts, args):
    if opts.input_fasta == None or opts.output_file == None:
       return True
    else:
       return False

class FastaRecord(object):
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

#    return FastaRecord(title, sequence)

def read_fasta_records(input_file):
    records = []
    sequence=""
    name=""
    while 1:
         line = input_file.readline()
         if line == "": 
            if sequence!="" and name!="":
               records.append(FastaRecord(name, sequence))
            return  records

         if line=='\n':
            continue

         line = line.rstrip()
         if  line.startswith(">") :
            if sequence!="" and name!="":
               records.append(FastaRecord(name, sequence))

            name = line.rstrip()
            sequence =""
         else:
            sequence = sequence + line.rstrip()
    return records

def format_db_blast(formatdb_executable, seq_subset_file):
    cmd='%s -dbtype prot -in %s' %(formatdb_executable, seq_subset_file.name)
    result= getstatusoutput(cmd)
    

def format_db_last(formatdb_executable, seq_subset_file):
    dirname = os.path.dirname(seq_subset_file.name)     
    cmd='%s -p -c %s  %s' %(formatdb_executable, dirname + PATHDELIM + 'subset_db', seq_subset_file.name)
    result= getstatusoutput(cmd)
    

def blast_against_itself(blast_executable, seq_subset_file, blast_table_out):
    cmd='%s -outfmt 6 -db  %s -query %s -out  %s' %(blast_executable,  seq_subset_file.name, seq_subset_file.name, blast_table_out)
    result= getstatusoutput(cmd)

def last_against_itself(last_executable, seq_subset_file, last_table_out):
    dirname = os.path.dirname(seq_subset_file.name)     
    cmd='%s -o %s -f 0 %s %s' %(last_executable,  last_table_out, dirname + PATHDELIM + 'subset_db',  seq_subset_file.name)
    result= getstatusoutput(cmd)


def add_last_refscore_to_file(blast_table_out, refscore_file, allNames):
    commentPATTERN = re.compile(r'^#')

    infile = open( blast_table_out,'r')
    refscores = {}
    lines = infile.readlines()
    for line in lines:
       if commentPATTERN.match(line):
          continue
       line=line.rstrip()
       fields = line.split('\t')
       if len(fields) != 12:
          print 'Error in the blastout file'
          sys.exit(1)
       if fields[6].rstrip()==fields[1].rstrip():
      #    fprintf(refscore_file, "%s\t%s\n",fields[0], fields[11])
          refscores[fields[1]]=fields[0]

    for key, value in refscores.iteritems():
       allNames[key] = True
       fprintf(refscore_file, "%s\t%s\n",key, value)

    infile.close()



def add_blast_refscore_to_file(blast_table_out, refscore_file, allNames):
    infile = open( blast_table_out,'r')
    refscores = {}
    lines = infile.readlines()
    for line in lines:
       line=line.rstrip()
       fields = line.split('\t')
       if len(fields) != 12:
          print 'Error in the blastout file'
          sys.exit(1)
       if fields[0].rstrip()==fields[1].rstrip():
      #    fprintf(refscore_file, "%s\t%s\n",fields[0], fields[11])
          refscores[fields[0]]=fields[11]

    for key, value in refscores.iteritems():
       allNames[key] = True
       fprintf(refscore_file, "%s\t%s\n",key, value)

    infile.close()
        

# compute the refscores
def compute_refscores(formatdb_executable, blast_executable,seq_subset_file, refscore_file, allNames, algorithm):
    if algorithm =='LAST':
        format_db_last(formatdb_executable, seq_subset_file)
        last_table_out = seq_subset_file.name + ".lastout"
        last_against_itself(blast_executable, seq_subset_file, last_table_out)
        add_last_refscore_to_file(last_table_out,refscore_file, allNames)

    if algorithm =='BLAST':
       format_db_blast(formatdb_executable, seq_subset_file)
       blast_table_out = seq_subset_file.name + ".blastout"
       blast_against_itself(blast_executable, seq_subset_file, blast_table_out)
       add_blast_refscore_to_file(blast_table_out,refscore_file, allNames)


    return None

def remove_blast_index_files(filename):
    prefixes = [ 'blastout', 'phr', 'pin', 'psq' ] 
    for prefix in prefixes:
       try:
          remove(filename +"." + prefix)
       except IOError:
          pass


def remove_last_index_files(filename):
    prefixes = [ 'des', 'sds', 'suf', 'bck',  'ssp', 'tis' ]

    remove( filename+ '.lastout')
    dirname = os.path.dirname(filename)     
    remove( dirname + PATHDELIM + 'subset_db0' +'.prj')
    for prefix in prefixes:
       try:
          remove(dirname + PATHDELIM + 'subset_db0.' + prefix)
       except IOError:
          pass



# the main function
SIZE = 1000

def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)
    if check_arguments(opts, args):
       print usage
       sys.exit(0)

    input_fasta = opts.input_fasta
    output_file = opts.output_file
    blast_executable = opts.blast_executable
    formatdb_executable = opts.formatdb_executable
    algorithm = opts.algorithm
 
    # input file to blast with itself to commpute refscore
    infile = open(input_fasta,'r')
   
    #this file has the refscores of the entire file
    outfile = open(output_file, 'w') 

    count = 0

    allNames= dict()
    for record in read_fasta_records(infile):
        if count % SIZE == 0:
            if count > 0:
              seq_subset_file.close()
              compute_refscores(formatdb_executable, blast_executable,seq_subset_file, outfile, allNames, algorithm);

              # now remove the old file
              if algorithm == 'BLAST' :
                 remove_blast_index_files(seq_subset_file.name)

              if algorithm == 'LAST' :
                 remove_last_index_files(seq_subset_file.name)

              remove(seq_subset_file.name)

            seq_subset_file = open(output_file +'.tmp.'+ str(count) +'.fasta','w')
        allNames[record.name.replace(">","")] = False;    
        fprintf(seq_subset_file, "%s\n", record.name)
        fprintf(seq_subset_file, "%s\n", record.sequence)

        count = count + 1

    #print str(count) + "   "  + "going to blast last sequence "
    if (count) % SIZE != 0:
       #print str(count) + "   "  + "last sequence "
       seq_subset_file.close()
       compute_refscores(formatdb_executable, blast_executable,seq_subset_file, outfile, allNames, algorithm);
       remove(seq_subset_file.name)
       if algorithm == 'BLAST' :
          remove_blast_index_files(seq_subset_file.name)
       if algorithm == 'LAST' :
          remove_last_index_files(seq_subset_file.name)


    #print count
    for key in allNames:
        if allNames[key] ==False:
           fprintf(outfile, "%s\t%s\n",key, 1000000)

    outfile.close()

def MetaPathways_refscore(argv, errorlogger = None, runstatslogger = None):
    createParser( )
    errorlogger.write("#STEP\tCOMPUTE_REFSCORE\n")
    main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger)
    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

