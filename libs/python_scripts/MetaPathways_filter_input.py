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
     from os import makedirs, sys, remove, rename
     from sys import path
     from optparse import OptionParser

     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, fprintf
     from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
     from libs.python_modules.parsers.fastareader  import FastaReader
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     sys.exit(3)

PATHDELIM = pathDelim()



usage= sys.argv[0] + """ -i file.fna  --min_length N --log_file logfile.log """ +\
              """ -o outfasta  [ -M map file ]"""

parser = None
def createParser():
    global parser
    epilog = """
This script filters both nucleotide and amino acid sequences.   In the case
of nucleotide sequences it filters out the contigs based on the minimum
length, presence of ambiguous nucleotide positions. Contigs with ambiguous
sequences are split into, separate contigs, at the regions where ambiguous
sequences appear. The resulting conigs are filtered out if they are below the
minimum length threshold.  The resulting filtered sequences are made available
into the folder 'preprocessed'   In the case of amino acid sequences in the
folder orf_prediction produced by the GFF to Amino step are filtered by
length. Unless the user changes, the default minimum length is 60 amino acid
sequences. The resulting sequences are put in the folder orf_prediction in a
file names "samplename".qced.faa. These amino acid sequences can be used in
 the downstream process of annotating them functionally and taxonomically
(using LCA rule)"""

    epilog = re.sub(r'[ \t\f\v]+',' ', epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-i", "--input_file", dest="input_fasta",
                      help='the input fasta filename [REQUIRED]')
    parser.add_option("-o", "--output_file", dest="output_fasta",
                      help='the output fasta filename [REQUIRED]')
    parser.add_option("-L", "--lengths_file", dest="lengths_file",
                      help='the output file that the lengths [REQUIRED]')
    parser.add_option("-m", "--min_length", type='int', dest="min_length",
                      help='minmum length of the sequence allowed')
    parser.add_option("-l", "--log_file", dest="log_file",  
                      help='file name to write the statsitics and log into')
    parser.add_option("-M", "--map_file", dest="map_file", type='str',  
                      help='file name to store the sequence name maps')
    parser.add_option("-t", "--type", dest="seqtype", type='str', default ='nucleotide',
                      help='the type of sequences,  choices are [ nucleotide, amino]')


def valid_arguments(opts, args):
    state = True
    if opts.input_fasta == None :
        print 'ERROR: Missing input fasta file'
        state = False

    if opts.output_fasta == None :
        print 'ERROR: Missing output fasta file'
        state = False

    if opts.min_length == None :
        print 'ERROR: Missing minimum sequence length'
        state = False

    if opts.log_file == None :
        print 'ERROR: Missing a log filename'
        state = False

    return state




def isAminoAcidSequence(sequence):
    if sequence:
        count = 0 

        list= {
           'A':  0,  'R':  0,  'N':  0,  'D':  0,  'C':  0,  'Q':  0,  'E':  0,  'G':  0,  
           'H':  0,  'I':  0,  'L':  0,  'K':  0,  'M':  0,  'F':  0,  'P':  0, 'S':  0,  
           'T':  0, 'W':  0, 'Y':  0, 'V':  0,  'B':  0,  'J':  0,  'Z':  0,  }

        for x in sequence:
            if x.upper() in list:
               list[x.upper()]=1

        count = 0 
        for x in list:
           count += list[x]

        if count > 10: 
            return True
        else:
             return False
    return True
       

def filter_sequence(sequence):
   if isAminoAcidSequence(sequence):
       return sequence

   sequence = re.sub(r'[^atcgATCG]','-', sequence.strip())
   subsequences =  sequence.split('-')

   max_length = 0;
   longest_sequence = ""; 
   for seq  in subsequences: 
      if len(seq) > max_length :
          longest_sequence = seq
          max_length = len(seq)

   return  longest_sequence



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

        

# the main function
SIZE = 1000

def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    if not valid_arguments(opts, args):
       print usage
       sys.exit(0)

    min_length = opts.min_length
    outfile = open(opts.output_fasta + '.tmp', 'w') 
    logfile = open(opts.log_file, 'w') 
    lengthsfile = open(opts.lengths_file + '.tmp', 'w') 
     

    if opts.map_file:
       mapfile = open(opts.map_file, 'w') 
    else:
       mapfile = None

    
    sample_name = opts.input_fasta;
    sample_name = re.sub(r'^.*/','',sample_name, re.I)
    sample_name = re.sub(r'^.*\\','',sample_name, re.I)
    sample_name = re.sub(r'\.fasta$','',sample_name, re.I)
    sample_name = re.sub(r'\.fna$','',sample_name, re.I)
    sample_name = re.sub(r'\.faa$','',sample_name, re.I)
    sample_name = re.sub(r'\.fas$','',sample_name, re.I)
    sample_name = re.sub(r'\.fa$','',sample_name, re.I)

    BEFORE = 'BEFORE'
    AFTER = 'AFTER'
    NUMSEQ = "#INFO\tNumber of sequences :"   
    NUMSEQ_SHORTER = "@INFO\tNumber of sequences shorter than minimum length of sequences"
    AVG_LENGTH= "@INFO\tAverage length of sequences:"
    MIN_LENGTH= "@INFO\tMinimum length of sequences:"
    MAX_LENGTH= "@INFO\tMaximum length of sequences:" 

    _MAX = 1000000000000
    stats = { 
              MIN_LENGTH: { 'BEFORE':_MAX, 'AFTER':_MAX },  
              MAX_LENGTH: { 'BEFORE': 0, 'AFTER':0 },  
              NUMSEQ : { 'BEFORE' :0, 'AFTER':0},   
              NUMSEQ_SHORTER : { 'BEFORE':0, 'AFTER':0 },
              AVG_LENGTH : { 'BEFORE':0, 'AFTER':0 },
            }  

    length_distribution = {}
    length_cumulative_distribution = {}

    for i in range(0,31):
        length_distribution[i]= 0
        length_cumulative_distribution[i]= 0

    seq_count = 0
    allNames= dict()
    outputStr = ""
    outputLines = []
    fastareader= FastaReader(opts.input_fasta)

    """ process one fasta sequence at a time """
    lengths_str=""
    for record in fastareader:
        seqname = record.name
        seq = record.sequence
        length = len(seq)
        
        index = int(len(seq) / 50);
        if index >= 30:
            index = 30

        length_distribution[index] += 1
        if length < stats[MIN_LENGTH][BEFORE] :
            stats[MIN_LENGTH][BEFORE] = length

        if length > stats[MAX_LENGTH][BEFORE] : 
            stats[MAX_LENGTH][BEFORE] = length

        if length < MIN_LENGTH:
            stats[NUMSEQ_SHORTER][BEFORE] += 1

        stats[AVG_LENGTH][BEFORE]  =  stats[AVG_LENGTH][BEFORE] + length

        seqvalue = filter_sequence(seq)
    
        stats[NUMSEQ][BEFORE] += 1
        
        seqlen = len(seqvalue)
        if seqlen>= min_length :

           if len(lengths_str) > 100: 
              fprintf(lengthsfile,"%s\n",lengths_str);
              lengths_str = str(seqlen)
           else:
              lengths_str += '\t' + str(seqlen)

           stats[NUMSEQ][AFTER] += 1
           stats[AVG_LENGTH][AFTER]  =  stats[AVG_LENGTH][AFTER] + seqlen
           if mapfile==None:
              fprintf(outfile, "%s\n", seqname)
           else:
               fprintf(outfile, ">%s\n",  sample_name + '_' + str(seq_count) )
               key = re.sub(r'^>','',seqname)
               fprintf(mapfile, "%s\n", sample_name+ '_' + str(seq_count) + '\t' + key + '\t' + str(seqlen))
               seq_count += 1

           fprintf(outfile, "%s\n",seqvalue)

           if  seqlen < stats[MIN_LENGTH][AFTER] :
               stats[MIN_LENGTH][AFTER] = seqlen
             
           if  seqlen > stats[MAX_LENGTH][AFTER] :
               stats[MAX_LENGTH][AFTER] = seqlen

    fprintf(lengthsfile,"%s\n",lengths_str);
    
    if stats[NUMSEQ][BEFORE] > 0 :
      stats[AVG_LENGTH][BEFORE]  = stats[AVG_LENGTH][BEFORE]/stats[NUMSEQ][BEFORE]
    else:
      stats[AVG_LENGTH][BEFORE]  = 0
    if stats[NUMSEQ][AFTER] > 0 :
       stats[AVG_LENGTH][AFTER]  = stats[AVG_LENGTH][AFTER]/stats[NUMSEQ][AFTER]
    else :
       stats[AVG_LENGTH][AFTER]  = 0

    lengthsfile.close()
    outfile.close()

    rename(opts.output_fasta + ".tmp", opts.output_fasta)
    rename(opts.lengths_file + ".tmp", opts.lengths_file)

    #inputfile.close()
    if mapfile != None:
       mapfile.close()

    """ min length """
    if stats[MIN_LENGTH][BEFORE] == _MAX:
       stats[MIN_LENGTH][BEFORE] = 0
    if stats[MIN_LENGTH][AFTER] == _MAX:
       stats[MIN_LENGTH][AFTER] = 0

    fprintf(logfile, "@INFO\tBEFORE\tAFTER\n");
    fprintf(logfile, "%s\n", NUMSEQ +'\t' + str(stats[NUMSEQ][BEFORE]) + '\t' + str(stats[NUMSEQ][AFTER]));
    fprintf(logfile, "%s\n", NUMSEQ_SHORTER   + '\t'+ str(stats[NUMSEQ_SHORTER][BEFORE]) + '\t' + str(stats[NUMSEQ_SHORTER][AFTER]))
    fprintf(logfile, "%s\n", AVG_LENGTH +'\t' + str(stats[AVG_LENGTH][BEFORE]) + '\t'+ str(stats[AVG_LENGTH][AFTER]))
    fprintf(logfile, "%s\n", MIN_LENGTH + '\t' + str(stats[MIN_LENGTH][BEFORE]) +'\t'+ str(stats[MIN_LENGTH][AFTER]))
    fprintf(logfile, "%s\n", MAX_LENGTH +'\t'+ str(stats[MAX_LENGTH][BEFORE]) + '\t' +  str(stats[MAX_LENGTH][AFTER]))
    fprintf(logfile, "@INFO\tLOW\tHIGH\tFREQUENCY\tCUMULATIVE_FREQUENCY\n");
#    fprintf(logfile, "#   ---\t-----\t--------\t---------\t----------\n");

    i  = 30
    length_cumulative_distribution[i] = length_cumulative_distribution[i];
    i  -= 1
    while i >= 0:
       length_cumulative_distribution[i] = length_cumulative_distribution[i+1] + length_distribution[i];
       i -= 1

    for i in range(0,31):
       fprintf(logfile, "   %s\n", str(i*50) + '\t' + str((i+1)*50) + '\t' +\
                str(length_distribution[i]) +'\t' + str(length_cumulative_distribution[i]) )

    logfile.close()


    if opts.seqtype=='nucleotide':
       priority = 1000
    else:
       priority = 2000

    if runstatslogger != None:
       if opts.seqtype=='nucleotide':
         runstatslogger.write("%s\tNumber of sequences in input file BEFORE QC (%s)\t%s\n" %(str(priority), opts.seqtype,  str(stats[NUMSEQ][BEFORE])) )
         runstatslogger.write("%s\t-min length\t%s\n" %(str(priority + 1), str(stats[MIN_LENGTH][BEFORE])) )
         runstatslogger.write("%s\t-avg length\t%s\n" %( str(priority + 2), str(int(stats[AVG_LENGTH][BEFORE]))))
         runstatslogger.write("%s\t-max length\t%s\n" %(str(priority + 3), str(stats[MAX_LENGTH][BEFORE])) )
         runstatslogger.write("%s\t-total base pairs (bp)\t%s\n" %(str(priority + 4), str(int(stats[AVG_LENGTH][BEFORE]* stats[NUMSEQ][BEFORE]))))

         runstatslogger.write("%s\tNumber of sequences AFTER QC (%s)\t%s\n" %(str(priority + 5), opts.seqtype, str(stats[NUMSEQ][AFTER])))
         runstatslogger.write("%s\t-min length\t%s\n" %(str(priority + 6), str(stats[MIN_LENGTH][AFTER])) )
         runstatslogger.write("%s\t-avg length\t%s\n" %( str(priority + 7), str(int(stats[AVG_LENGTH][AFTER]))))
         runstatslogger.write("%s\t-max length\t%s\n" %( str(priority + 8), str(stats[MAX_LENGTH][AFTER])) )
         runstatslogger.write("%s\t-total base pairs (bp)\t%s\n" %( str(priority + 9), str(int(stats[AVG_LENGTH][AFTER]* stats[NUMSEQ][AFTER])) ))
       else:
         runstatslogger.write("%s\tNumber of translated ORFs BEFORE QC (%s)\t%s\n" %(str(priority), opts.seqtype,  str(stats[NUMSEQ][BEFORE])) )
         runstatslogger.write("%s\t-min length\t%s\n" %(str(priority + 1), str(stats[MIN_LENGTH][BEFORE])) )
         runstatslogger.write("%s\t-avg length\t%s\n" %( str(priority + 2), str(int(stats[AVG_LENGTH][BEFORE]))))
         runstatslogger.write("%s\t-max length\t%s\n" %(str(priority + 3), str(stats[MAX_LENGTH][BEFORE])) )
         runstatslogger.write("%s\t-total base pairs (bp)\t%s\n" %(str(priority + 4), str(int(stats[AVG_LENGTH][BEFORE]* stats[NUMSEQ][BEFORE]))))
         runstatslogger.write("%s\tNumber of tranlated ORFs AFTER QC (%s)\t%s\n" %(str(priority + 5), opts.seqtype, str(stats[NUMSEQ][AFTER])))
         runstatslogger.write("%s\t-min length\t%s\n" %(str(priority + 6), str(stats[MIN_LENGTH][AFTER])) )
         runstatslogger.write("%s\t-avg length\t%s\n" %( str(priority + 7), str(int(stats[AVG_LENGTH][AFTER]))))
         runstatslogger.write("%s\t-max length\t%s\n" %( str(priority + 8), str(stats[MAX_LENGTH][AFTER])) )
         runstatslogger.write("%s\t-total base pairs (bp)\t%s\n" %( str(priority + 9), str(int(stats[AVG_LENGTH][AFTER]* stats[NUMSEQ][AFTER])) ))



def MetaPathways_filter_input(argv, errorlogger = None, runstatslogger = None):
    createParser()
    main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger) 
    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

