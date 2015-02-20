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
     import os, re, traceback
     from os import makedirs, sys, remove, rename
     from sys import path
     from optparse import OptionParser

     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, fprintf
     from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     sys.exit(3)

PATHDELIM = pathDelim()

usage = __file__ + """ -i input_fasta_file -o output_file [-B blast_executable -F formatdb_executable] """
parser = None
def createParser():
    global parser
    epilog = """The amino acid sequences in the orf_prediction folder are used to do a self BLAST/LAST, which will be used to compute the bit score ratio (BSR) for the hits.  The BSR ratio can be defined at the ratio of a the bit-score between a query and a target sequence to the bitcore when both the query and target sequenes are the query sequence. Usually, a BSR ratio of 0.4 or more is considered as a good hit for protein sequences. Note that BSR ratio is designed in some sense to have a normalized value for the bit-score  since the score is also influenced by the length of the query. 
The results are written to a file  (usually in a folder called blast_results in the MetaPathway pipeline,  into a file named <samplename>.refscore.<algorithm> (where <algorithm> refers to the BLAST or LAST in the context of the pipeline) extension This script can be extended to add other sequence homology search algorithms."""

    epilog = re.sub(r'[ \t\f\v]+',' ', epilog)
    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-i", "--input_file", dest="input_fasta",
                      help='the input fasta file [REQUIRED]')
    parser.add_option("-o", "--output_file", dest="output_file",
                      help='the output fasta file [REQUIRED]')
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

class FastaReader():
    """Parses a fasta record from a string or file."""
    stop = False
    START_PATTERN = re.compile(r'^>')
    name = None
    future_name =None
    sequence=""

    def __init__(self, fasta_filename):
        try:
            self.file = open(fasta_filename, 'r')
        except IOError:
            print "Cannot open fasta file " + fasta_filename

    def __iter__(self):
        return self

    def close(self):
         self.file.close()

    def next(self):
        if self.stop:
          raise StopIteration

        try:
           if not self.name:
               self.name = self.file.readline().strip()
           line = self.file.readline()
        except:
           line = None

        if not line:
           self.stop = True
           raise StopIteration

        fragments = []
        while line and not self.START_PATTERN.search(line):
            fragments.append(line.strip())
            line = self.file.readline()

       # print line
        if self.future_name:
            self.name = re.sub('>','',self.future_name)

        if line:
          self.future_name = line.strip()

        self.sequence =''.join(fragments)
        self.seqname =  self.name
    
        return FastaRecord( re.sub('>', '', self.name), self.sequence)


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
        

# write the refscores
def write_refscores(refscore_file, refscores):
    for key, value in refscores.iteritems():
       fprintf(refscore_file, "%s\t%s\n",key, value)



SCORES= {
'A':  4,
'R':  5,
'N':  6,
'D':  6,
'C':  9,
'Q':  5,
'E':  5,
'G':  6,
'H':  8,
'I':  4,
'L':  4,
'K':  5,
'M':  5,
'F':  6,
'P':  7 ,
'S':  4 ,
'T':  5 ,
'W':  11,
'Y':  7 ,
'V':  4,
'B':  4,
'J':  3,
'Z':  4,
'X': -1,
'*': 1,
}

def getrefscore(seq):
    score =0
    for c in seq:
      try:
        score += SCORES[c]
      except:
        score = 0
    return score 

def compute_refscores(sequences_subset, refscore_file):
    refscores ={} 
    for key, value in sequences_subset.iteritems():
       refscores[key] = getrefscore(value)
    write_refscores(refscore_file, refscores)


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
def _old_compute_refscores(formatdb_executable, blast_executable,seq_subset_file, refscore_file, allNames, algorithm):
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
    suffixes = [ 'prj', 'des', 'sds', 'suf', 'bck',  'ssp', 'tis' ]

    remove( filename+ '.lastout')
    dirname = os.path.dirname(filename)     
    for suffix in suffixes:
       try:
          remove(dirname + PATHDELIM + 'subset_db.' + suffix)
       except IOError:
          pass



# the main function
SIZE = 10000

def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)
    if check_arguments(opts, args):
       print usage
       sys.exit(0)

    input_fasta = opts.input_fasta
    output_file = opts.output_file
    algorithm = opts.algorithm
 
    # input file to blast with itself to commpute refscore
   
    #this file has the refscores of the entire file
    outfile = open(output_file + ".tmp", 'w') 

    count = 0

    allNames= dict()
    sequence_subset = dict() 
    refscores = dict() 

    fastaReader = FastaReader(input_fasta)

    for record in fastaReader:
       count = count + 1
       sequence_subset[record.name] = record.sequence
       if count % SIZE == 0:
          compute_refscores(sequence_subset, outfile)
          count  = 0
          sequence_subset = dict() 

    if count % SIZE != 0:
        compute_refscores(sequence_subset, outfile)
        sequence_subset = dict() 
    #print count

    fastaReader.close()
    outfile.close()
    rename(output_file + ".tmp", output_file) 

def MetaPathways_refscore(argv, errorlogger = None, runstatslogger = None):
    createParser( )
    if errorlogger:
       errorlogger.write("#STEP\tCOMPUTE_REFSCORE\n")
    try:
       main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger)
    except:
       print 'error'
       print traceback.format_exc(10)

       return (1,traceback.format_exc(10))

    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

