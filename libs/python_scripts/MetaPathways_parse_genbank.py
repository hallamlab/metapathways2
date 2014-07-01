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
     from os import makedirs, sys, remove, _exit
     from sys import path
     import re
     from optparse import OptionParser, OptionGroup
     from libs.python_modules.utils.utils  import *
     from libs.python_modules.annotate.sequence import  genbank, fasta
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)
usage= sys.argv[0] + """ -g/--genbank  genbank --output-gff output_gff_file --output-fna output_fna_file --output-faa output_faa_file """
parser = None

def createParser():
    global parser
    parser = OptionParser(usage)

    parser.add_option("-g", "--genbank", dest="genbank", 
                      help='the genbank file [ REQUIRED]')
    
    parser.add_option("--output-fna", dest="output_fna",
                      help='the output nucleotide sequences [REQUIRED]')
    
    parser.add_option("--output-faa", dest="output_faa",
                      help='the output amino acid sequences [REQUIRED]')
    
    parser.add_option("--output-gff", dest="output_gff",
                      help='the output gff file for the orfs [REQUIRED]')
    
    parser.add_option("-M", "--map_file",  dest="map_file",
                      help='file name to store the sequences name maps [REQUIRED]')

    parser.add_option("-L", "--lengths_file",  dest="lengths_file",
                      help='file name to store the sequences name maps [REQUIRED]')

    parser.add_option( "--create-functional-table",  dest="functional_table",
                      help='file name to store the functional table  [OPTIONAL]')
    

def check_arguments(opts, args):

    if opts.genbank==None:
         print "A genbank file should be specified"  
         return False

    if opts.output_fna==None:
         print "The output nucloetide sequences"  
         return False

    if opts.output_faa==None:
         print "The output amino acid sequences"  
         return False

    if opts.output_gff==None:
         print "The output amino acid sequences"  
         return False
    return True

def create_dictionary(databasemapfile, annot_map):
       seq_beg_pattern = re.compile(">")

       dbmapfile = open( databasemapfile,'r')
       lines=dbmapfile.readlines()
       dbmapfile.close()
       for line in lines:
          if seq_beg_pattern.search(line):
              words = line.rstrip().split()
              name = words[0].replace('>','',1)
               
              words.pop(0)
              annotation = ' '.join(words)
              annot_map[name]= annotation
       
def process_file(genbank_filename, output_fna, output_faa, output_gff,\
                 map_file, functional_table, lengths_file):
   """This script converts GenBank entries to PathoLogic files, 
   creating a .pf and .fna file for each GenBank entry."""


   try:
      functional_table_fp = None
      if functional_table!=None:
        functional_table_fp = open(functional_table,"w")
        fprintf(functional_table_fp, "ORF_ID\tORF_length\tstart\tend\tContig_Name\t" +\
                                  "Contig_length\tstrand\tec\ttaxonomy\tproduct\n")
   except:
       eprintf("ERROR:\tCannot open functional table file %s\n", functional_table)
       functional_table_fp = None

   try:
      lengths_file_fp = None
      if lengths_file!=None:
        lengths_file_fp = open(lengths_file,"w")
   except:
       eprintf("ERROR:\tCannot open contigs lengths file %s\n", lengths_file)
       functional_table_fp = None




   map_file_fp = open(map_file, "w");
   
   with open(genbank_filename, 'r') as genbank_file:

      outputfnafile = open(output_fna,'w')
      outputfaafile = open(output_faa,'w')
      outputgfffile = open(output_gff,'w')
      record_count = 0;
      sample_name = extractSampleName(genbank_filename, 'gbk-unannotated')

      for record in genbank.GenBankRecordParser(genbank_file.read()):
        _locus = sample_name + '_' + str(record_count)

        fprintf(outputfnafile, ">%s\n%s\n", _locus, record.sequence)
        feature_count = 0
        contig_length = len(record.sequence)
        fprintf(lengths_file_fp, "%s\n", contig_length)

        fprintf(map_file_fp, "%s\t%s\t%s\n",  _locus, record.locus, str(contig_length))

        for feature in record.features:
           if feature.type=='CDS':
              _locus_tag = _locus + '_' + str(feature_count)

              fprintf(outputfaafile, ">%s\n%s\n",  _locus_tag, feature.translation)

              start, end, strand = feature.coordinates 
              gff_Str = _locus 
              gff_Str += '\t' + 'Genbank file' 
              gff_Str += '\t' + 'CDS' 
              gff_Str += '\t' + start 
              gff_Str += '\t' + end 
              gff_Str += '\t' + '0'
              gff_Str += '\t' + strand
              gff_Str += '\t' + '0'

              try:
                 orf_length = int(end.strip()) - int(start.strip())
                 if orf_length < 0:
                    orf_length = - orf_length
              except:
                 orf_length = 0

 
              gff_Str +=  '\tID=' + _locus_tag
              gff_Str +=  ';locus_tag=' + _locus_tag
              gff_Str +=  ';orf_length=' +  str(orf_length)
              gff_Str +=  ';contig_length=' + str(contig_length)
              gff_Str +=  ';partial=00'
              if feature.product:  
                 gff_Str +=  ';product=' + feature.product
              else:
                 gff_Str +=  ';product=' + ''

              fprintf(outputgfffile, "%s\n",  gff_Str)

              ''' now prepare the functional stuff '''
              if functional_table_fp !=None:
                  funStr  = _locus_tag + "\t"
                  funStr += str(orf_length) + "\t"
                  funStr += start + "\t"
                  funStr += end + "\t"
                  funStr += _locus + "\t"
                  funStr += str(contig_length) + "\t"
                  funStr += strand + "\t"
                  funStr += " " + "\t"
                  funStr += " " + "\t"
                  funStr += feature.product
                  fprintf(functional_table_fp, funStr +"\n")

              feature_count+= 1
        record_count+= 1
            

      outputfnafile.close()
      outputgfffile.close()
      outputfaafile.close()
      map_file_fp.close()

      if functional_table_fp !=None:
         functional_table_fp.close()


def MetaPathways_parse_genbank(argv, errorlogger = None, runstatslogger = None):
     """ the main function of metapaths """
     createParser()
     main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger)
     return (0,'')


def main(argv, errorlogger = None, runstatslogger = None): 
    """ the main function """
    global parser 

    (opts, args) = parser.parse_args(argv )

    if not check_arguments(opts, args):
       print usage
       sys.exit(0)

    genbank  = opts.genbank
    output_fna =  opts.output_fna
    output_faa = opts.output_faa
    output_gff = opts.output_gff
    map_file = opts.map_file
    lengths_file = opts.lengths_file
    functional_table = opts.functional_table
    process_file(genbank, output_fna, output_faa, output_gff, map_file, functional_table, lengths_file)


# the main function of metapaths
if __name__ == "__main__":
    """ the main function of metapaths """
    createParser()
    main(sys.argv[1:])

