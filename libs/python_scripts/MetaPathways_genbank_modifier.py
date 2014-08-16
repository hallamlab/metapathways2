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
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     sys.exit(3)


usage= """./MetapathWays_parse_genbank.py --input-gbk genbank --output-gbk output 
                               [ --template-header template_file --input-gff input_gff] """
parser = OptionParser(usage)
parser.add_option( "--input-gbk", dest="input_gbk", 
                  help='the genbank file [ REQUIRED]')

parser.add_option("--output-gbk", dest="output_gbk", 
                  help='the genbank file [ REQUIRED]')

parser.add_option("--template-header", dest="template_header",
                  help='the template header')

parser.add_option("--input-gff", dest="input_gff",
                  help='the input gff for annotation')


def check_arguments(opts, args):

    if opts.input_gbk==None:
         print "A genbank file should be specified"  
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
       

    

def insert_attribute(attributes, attribStr):
     rawfields = re.split('=', attribStr)
     if len(rawfields) == 2:
       attributes[rawfields[0].strip().lower()] = rawfields[1].strip()



def split_attributes(str, attributes):        
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr) 

     return attributes
   
def insert_orf_into_dict(line, contig_dict):
     rawfields = re.split('\t', line)
     fields = [] 
     for field in rawfields:
        fields.append(field.strip());
     
     if( len(fields) != 9):
       return

     attributes = {}
     try:
         attributes['seqname'] =  fields[0]   # this is a bit of a  duplication  
         attributes['source'] =  fields[1] 
         attributes['feature'] =  fields[2] 
         attributes['start'] =  int(fields[3])
         attributes['end'] =  int(fields[4])
     except:
         print line
         print fields
         print attributes
         sys.exit(0)

     try:
        attributes['score'] =  float(fields[5])
     except:
        attributes['score'] =  fields[5]

     attributes['strand'] =  fields[6] 
     attributes['frame'] =  fields[7] 
     split_attributes(fields[8], attributes)
    

     if 'id' in attributes:
         if not attributes['id'] in contig_dict :
            contig_dict[attributes['id']] = {}
         contig_dict[attributes['id']]=attributes.copy()



def load_gff_file_to_dictionary(gff_file_name, gff_dict):
     try:
        gfffile = open(gff_file_name, 'r')
     except IOError:
        print "Cannot read file " + gff_file_name + " !"

     gff_lines = gfffile.readlines()
     gff_beg_pattern = re.compile("^#")
     gfffile.close()
     
     contig_dict={} 
     count = 0
     for line in gff_lines:
        line = line.strip() 
        if gff_beg_pattern.search(line):
          continue
        insert_orf_into_dict(line, gff_dict)


def read_template_header(headers,  template_header):
     template_header_file = open(template_header, 'r') 
     templatelines = template_header_file.readlines()

     headers['REFERENCES'] = []
     reference={}

     for line in templatelines:
         fields = [ x.strip() for x in line.split()]
         if len(fields)> 1:
            if fields[0] in ['REFERENCE', 'TITLE', 'JOURNAL', 'AUTHORS', 'PUBMED']:
                if fields[0]=='REFERENCE':
                  if not reference=={}:
                     headers['REFERENCES'].append(reference) 
                     reference={}
                  reference={'id':'.', 'consortium':'.', 'title':'.', 'journal':'.', 'authors':'.', 'pubmed':'.'}

                 
                secondToLastFields = fields[1:]
                if fields[0].lower()=='authors':
                   authorstring =' '.join(secondToLastFields) 
                   reference[fields[0].lower()] = authorstring.split(', ') 
                else:
                   reference[fields[0].lower()] = ' '.join(secondToLastFields) 
            else:
                secondToLastFields = fields[1:]
                headers[fields[0]] = ' '.join(secondToLastFields) 
           
     if not reference=={}:
         headers['REFERENCES'].append(reference) 
             


# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)

    headers={} 
    if opts.template_header!=None:
       read_template_header(headers, opts.template_header)
    gff_dictionary={}    


    if opts.input_gff!=None:
        load_gff_file_to_dictionary(opts.input_gff, gff_dictionary)

    process_gbk_file(opts.input_gbk, opts.output_gbk, headers, gff_dictionary)


"""This script converts GenBank entries to PathoLogic files, creating a .pf and
.fna file for each GenBank entry."""

def process_gbk_file(input_gbk, output_gbk, headers, gff_dictionary):
   
  tag = re.sub(r'[.]gbk','', input_gbk)
  tag = re.sub(r'.*/','', tag)

  output_gbk_file = open(output_gbk,'w') 
  serializer = genbank.GenBankRecordSerializer()
  with open(input_gbk, 'r') as genbank_file:
     out_list=[]
     count = 0
     for record in genbank.GenBankRecordParser(genbank_file.read()):
        count+=1
        
        record.locus = tag +  str(count)
        if count%1000==0:
           print 'Count = ' + str(count)
        
        if headers and 'REFERENCES' in headers:
           record.references_ = headers['REFERENCES']

        i = 0
        for feature in record.features:
           if feature.type =="CDS":
             if  feature.locus_tag in gff_dictionary:
                record.features[i].product= 'aaaaa ' + gff_dictionary[feature.locus_tag]['product']
           i+=1

        #record.locus = "hello"

        out_list.append(serializer.serialize(record))
        if count%1000 == 0:
           output_str = '\n'.join(out_list)
           out_list=[]
           fprintf(output_gbk_file,'%s\n',output_str)

     output_str = '\n'.join(out_list)
     fprintf(output_gbk_file,'%s\n',output_str)

     output_gbk_file.close()

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

