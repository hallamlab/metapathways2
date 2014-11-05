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

     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, fprintf, read_fasta_records
     from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
     from libs.python_modules.parsers.fastareader  import FastaReader
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     sys.exit(3)

PATHDELIM = pathDelim()



usage= sys.argv[0] + """ -i <file-amino>  --min_length N --log_file <logfile> """ +\
              """ --fasta <dummycontigs>  --faa <amino> --fna <dummyorfs> -L <lengthsfiles> [ -M <mapfile> ]"""

parser = None
def createParser():
    global parser
    epilog = """
This script processes the input when it is provided as amino acid sequences 
instead of nucloetide sequences. From the input it creates the following 
files to conform  to the format as in the case of nucleotide sequences input:

 1. preprocess/sample.fasta

 2. preprocess/sample.mapping.txt

 3. run_statistics/orf_prediction

"""

    epilog = re.sub(r'[ \t\f\v]+',' ', epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-i", "--input_file", dest="input_fasta",
                      help='the input fasta filename [REQUIRED]')

    parser.add_option("--fasta", dest="output_fasta",
                      help='the dummy output fasta filename [REQUIRED]')
    parser.add_option("--faa", dest="output_faa",
                      help='the amino acid  fasta filename [REQUIRED]')
    parser.add_option("--fna", dest="output_fna",
                      help='the amino acid  fna filename [REQUIRED]')
    parser.add_option("--gff", dest="output_gff",
                      help='the gff file for the amino acid sequences filename [REQUIRED]')


    parser.add_option("-L", "--lengths_file", dest="lengths_file",
                      help='the output file that the lengths [REQUIRED]')
    parser.add_option("-l", "--log_file", dest="log_file",  
                      help='file name to write the statsitics and log into')
    parser.add_option("-M", "--map_file", dest="map_file", type='str',  
                      help='file name to store the sequence name maps')


def valid_arguments(opts, args):
    state = True
    if opts.input_fasta == None :
        print 'ERROR: Missing input fasta file'
        state = False

    if opts.output_fasta == None :
        print 'ERROR: Missing dummy output fasta file'
        state = False

    if opts.output_faa == None :
        print 'ERROR: Missing amino acid sequences  fasta file'
        state = False

    if opts.output_fna == None :
        print 'ERROR: Missing ORF sequences  fasta file'
        state = False

    if opts.log_file == None :
        print 'ERROR: Missing a log filename'
        state = False

    return state



def isAminoAcidSequence(sequence):
    if sequence:
        count = 0 
        list = [ 'a', 't', 'c', 'g', 'A', 'T', 'C', 'G']
        for x in sequence:
            if x in list:
               count+=1
        if count/len(sequence) < 0.80:
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


# the main function
SIZE = 1000

def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    if not valid_arguments(opts, args):
       print usage
       sys.exit(0)

    min_length = 0
    #inputfile = open(opts.input_fasta,'r')
    outfile = open(opts.output_fasta, 'w') 
    outfilefna = open(opts.output_fna, 'w') 
    outfilefaa = open(opts.output_faa, 'w') 
    outfilegff = open(opts.output_gff, 'w') 


    logfile = open(opts.log_file, 'w') 
    lengthsfile = open(opts.lengths_file, 'w') 
     
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
               contigID =  sample_name + '_' + str(seq_count) 
               orfID =  sample_name + '_' + str(seq_count) + "_0" 

               fprintf(outfile, ">%s\n",  contigID )
               fprintf(outfilefna, ">%s\n",  orfID )
               fprintf(outfilefaa, ">%s\n",  orfID )

               gffString =  sample_name + '_' + str(seq_count)
               gffString +=  "\t" + "AMINO_ACID_SEQ"
               gffString +=  "\t" + "CDS"
               gffString +=  "\t" + "0"
               gffString +=  "\t" + str(3*seqlen)
               gffString +=  "\t" + "0"
               gffString +=  "\t" + "+"
               gffString +=  "\t" + "0"
               gffString +=  "\t" + "ID=" + orfID + ";" 
               gffString +=  "locus_tag=" + orfID + ";" 
               gffString +=  "partial=00;" 
               gffString +=  "orf_length="+ str(seqlen)+";" 
               gffString +=  "contig_length="+ str(3*seqlen)

               fprintf(outfilegff, "%s\n", gffString)

               key = re.sub(r'^>','',seqname)
               fprintf(mapfile, "%s\n", sample_name+ '_' + str(seq_count) + '\t' + key + '\t' + str(seqlen))
               seq_count += 1

           fprintf(outfile, "%s\n","DUMMY CONTIGS FOR AMINO ACID SEQUENCES")
           fprintf(outfilefna, "%s\n","DUMMY ORFS FOR AMINO ACID SEQUENCES")
           fprintf(outfilefaa, "%s\n",seqvalue)

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
    outfilefna.close()
    outfilefaa.close()
    outfilegff.close()

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


    seqtype='amino'
    """priority is used to sort the output to print in the right order"""
    priority = 2000

    if runstatslogger != None:
         runstatslogger.write("%s\tSequences BEFORE Filtering (%s)\t%s\n" %(str(priority), seqtype,  str(stats[NUMSEQ][BEFORE])) )
         runstatslogger.write("%s\tmin length\t%s\n" %(str(priority + 1), str(stats[MIN_LENGTH][BEFORE])) )
         runstatslogger.write("%s\tavg length\t%s\n" %( str(priority + 2), str(int(stats[AVG_LENGTH][BEFORE]))))
         runstatslogger.write("%s\tmax length\t%s\n" %(str(priority + 3), str(stats[MAX_LENGTH][BEFORE])) )
         runstatslogger.write("%s\ttot length\t%s\n" %(str(priority + 4), str(int(stats[AVG_LENGTH][BEFORE]* stats[NUMSEQ][BEFORE]))))
         runstatslogger.write("%s\tSequences AFTER Filtering (%s)\t%s\n" %(str(priority + 5), seqtype, str(stats[NUMSEQ][AFTER])))
         runstatslogger.write("%s\tmin length\t%s\n" %(str(priority + 6), str(stats[MIN_LENGTH][AFTER])) )
         runstatslogger.write("%s\tavg length\t%s\n" %( str(priority + 7), str(int(stats[AVG_LENGTH][AFTER]))))
         runstatslogger.write("%s\tmax length\t%s\n" %( str(priority + 8), str(stats[MAX_LENGTH][AFTER])) )
         runstatslogger.write("%s\ttot length\t%s\n" %( str(priority + 9), str(int(stats[AVG_LENGTH][AFTER]* stats[NUMSEQ][AFTER])) ))

def MetaPathways_preprocess_amino_input(argv, errorlogger = None, runstatslogger = None):
    createParser()
    main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger) 

    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

