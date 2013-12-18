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
     from python_modules.metapaths_utils  import parse_command_line_parameters, fprintf, printf
     from python_modules.sysutil import getstatusoutput
    # from python_modules.annotate.sequence import  fasta
#     import python_modules.annotate.sequence.genbank
     from python_modules.annotate.sequence import  genbank, fasta
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)


usage= """./MetapathWays_parse_genbank.py -g/--genbank  genbank --output-gff output_gff_file --output-fna output_fna_file --output-faa output_faa_file """
parser = OptionParser(usage)
parser.add_option("-g", "--genbank", dest="genbank", 
                  help='the genbank file [ REQUIRED]')

parser.add_option("--output-fna", dest="output_fna",
                  help='the output nucleotide sequences [REQUIRED]')

parser.add_option("--output-faa", dest="output_faa",
                  help='the output amino acid sequences [REQUIRED]')

parser.add_option("--output-gff", dest="output_gff",
                  help='the output gff file for the orfs [REQUIRED]')



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
       


# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)
    process_file(opts.genbank, opts.output_fna, opts.output_faa, opts.output_gff)






"""This script converts GenBank entries to PathoLogic files, creating a .pf and
.fna file for each GenBank entry."""

def process_file(genbank_filename, output_fna, output_faa, output_gff):
   
   with open(genbank_filename, 'r') as genbank_file:
      outputfnafile = open(output_fna,'w')
      outputfaafile = open(output_faa,'w')
      outputgfffile = open(output_gff,'w')
      for record in genbank.GenBankRecordParser(genbank_file.read()):
        fprintf(outputfnafile, ">%s\n%s\n", record.locus, record.sequence)
        for feature in record.features:
           if feature.type=='CDS':
              fprintf(outputfaafile, ">%s\n%s\n",  feature.locus_tag, feature.translation)

              start, end, strand = feature.coordinates 
              gff_Str = record.locus 
              gff_Str += '\t' + 'Genbank file' 
              gff_Str += '\t' + 'CDS' 
              gff_Str += '\t' + start 
              gff_Str += '\t' + end 
              gff_Str += '\t' + '0'
              gff_Str += '\t' + strand
              gff_Str += '\t' + '0'
 
              gff_Str +=  '\tID=' + feature.locus_tag
              gff_Str +=  ';locus_tag=' + feature.locus_tag
              if feature.product:  
                 gff_Str +=  ';product=' + feature.product
              else:
                 gff_Str +=  ';product=' + ''

              fprintf(outputgfffile, "%s\n",  gff_Str)
            

      outputfnafile.close()
      outputgfffile.close()
      outputfaafile.close()



# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

