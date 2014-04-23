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
     from os import makedirs, sys, remove
     from sys import path
     from glob import glob
     
     from optparse import OptionParser
     from python_modules.metapaths_utils  import parse_command_line_parameters, fprintf, eprintf
     from python_modules.sysutil import pathDelim
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)

PATHDELIM = pathDelim()

usage=  sys.argv[0] + """ -s sample_name  -f folder_path """ 

parser = OptionParser(usage)
parser.add_option("-s", "--sample_name", dest="sample_name",
                  help='the sample name [REQUIRED]')
parser.add_option("-f", "--folder_path", dest="folder_path",
                  help='the folder path [REQUIRED]')


def valid_arguments(opts, args):
    state = True
    if opts.sample_name == None :
        print 'Must have a sample name'
        state = False

    if opts.folder_path == None :
        print 'Must have a folder path'
        state = False

    return state



class FastaRecord(object):
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence


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

def get_number_of_BLAST_LAST_hits(file_name):
    commentPATTERN = re.compile(r'^#')
    count = 0
    try:
       inputfilename = open(file_name, 'r')
    except:
       return results 

    line =  inputfilename.readline()
    while line: 
       if commentPATTERN.search(line):
          line =  inputfilename.readline()
          continue
       fields= [x.strip() for x in line.split('\t') ] 
       if len(fields) < 12:
          line =  inputfilename.readline()
          continue
       count += 1
       line =  inputfilename.readline()
    inputfilename.close()

    
     
    return count

def get_number_of_rRNA_hits(file_name):
    commentPATTERN = re.compile(r'similarity')
    count = 0
    try:
       inputfilename = open(file_name, 'r')
    except:
       return results 

    line =  inputfilename.readline()
    while line: 
       if commentPATTERN.search(line):
          line =  inputfilename.readline()
          continue

       fields= [x.strip() for x in line.split('\t') ] 
       if len(fields) < 7:
          line =  inputfilename.readline()
          continue
       count += 1
       line =  inputfilename.readline()

    inputfilename.close()
    return count

def get_number_of_tRNA_hits(file_name):
    dataPATTERN = re.compile(r'number of predicted tRNA=(.*)')
    count = 0
    try:
       inputfilename = open(file_name, 'r')
    except:
       return count 

    line =  inputfilename.readline()
    while line: 
       if not dataPATTERN.search(line):
          line =  inputfilename.readline()
          continue

       result = dataPATTERN.search(line)
       if result:
          if len(result.groups())==1:
            count = result.group(1)
            return count
       line =  inputfilename.readline()

    inputfilename.close()
    return count


# get the rRNA_hits
def get_rRNA_hits(sample_name, folder_path):
    results = []
    regPattern = re.compile(r'.rRNA.stats.txt')
    input_dir = folder_path +  PATHDELIM + 'results' + PATHDELIM + 'rRNA' 
    files = [ re.sub(r'.*\/','',f) for f in glob(input_dir + PATHDELIM + sample_name + '*')  if regPattern.search(f) ]

    regPattern = re.compile(r'[.](.*)[.]rRNA.stats.txt')
    
    for file in files:
      result = regPattern.search(file)
      if result:
         database = result.group(1)
         file_name = input_dir + PATHDELIM + sample_name + '.' +  result.group(1) + '.rRNA.stats.txt'
         count  =  get_number_of_rRNA_hits(file_name)
         results.append( (  'Number of rRNA hits in ' + database, count ) )
    return results

# get the rRNA_hits
def get_tRNA_hits(sample_name, folder_path):
    results = []
    regPattern = re.compile(r'.tRNA.stats.txt')
    input_dir = folder_path +  PATHDELIM + 'results' + PATHDELIM + 'tRNA' 
    file_name = input_dir + PATHDELIM + sample_name + '.tRNA.stats.txt'
    count  =  get_number_of_tRNA_hits(file_name)
    results.append( ('Number of tRNA hits in ', count ) )
    return results


def get_number_of_uncommented_lines(file_name):
    commentPATTERN = re.compile(r'^#')
    count = 0
    try:
       inputfilename = open(file_name, 'r')
    except:
       return count 

    line =  inputfilename.readline()
    while line: 
       if commentPATTERN.search(line):
          line =  inputfilename.readline()
          continue
       fields= [x.strip() for x in line.split('\t') ] 
       count += 1
       line =  inputfilename.readline()

    inputfilename.close()
    return count

#counts the number of taxonomic and annotated ORFs
def get_functional_taxonomic_hits(sample_name, folder_path):
    results = []
    # for the LAST algorithm
    regPattern = re.compile(r'.annot.gff$', flags=re.IGNORECASE)
    input_dir = folder_path +  PATHDELIM + 'results' + PATHDELIM + 'annotation_table' 
    file_name = input_dir + PATHDELIM +  'functional_and_taxonomic_table.txt'

    eprintf("\nCounting number of functionally and taxonomically ORFs ...")
    count  =  get_number_of_uncommented_lines(file_name)
    eprintf("done\n")
    results.append( ('Total functionally and taxonmically  annotated ORFs', count ) )
    return results

#counts the number of ORFs in the table ORF_annotation_table
def get_ORF_annotations_hits(sample_name, folder_path):
    results = []
    # for the LAST algorithm
    regPattern = re.compile(r'.annot.gff$', flags=re.IGNORECASE)
    input_dir = folder_path +  PATHDELIM + 'results' + PATHDELIM + 'annotation_table' 
    file_name = input_dir + PATHDELIM +  'ORF_annotation_table.txt'

    eprintf("\nCounting number of ORFs for mapping to functional classification ...")
    count  =  get_number_of_uncommented_lines(file_name)
    eprintf("done\n")
    results.append( ('Total orfs count for functional classification', count ) )
    return results


#counts the number of annotatations generated
def get_annotation_hits(sample_name, folder_path):
    results = []
    # for the LAST algorithm
    regPattern = re.compile(r'.annot.gff$', flags=re.IGNORECASE)
    input_dir = folder_path +  PATHDELIM + 'genbank' 
    files = [ re.sub(r'.*\/','',f) for f in glob(input_dir + PATHDELIM + sample_name + '*')  if regPattern.search(f) ]
    regPattern = re.compile(r'(.*)[.]annot.gff$', flags=re.IGNORECASE)
    
    for file in files:
      result = regPattern.search(file)
      if result:
         file_name = input_dir + PATHDELIM + sample_name +    '.annot.gff'
         eprintf("\nCounting number of annotations...")
         count  =  get_number_of_uncommented_lines(file_name)
         eprintf("done\n")
         results.append( ('Total number of valid annotations', count ) )
    return results

# counts the number of parsed BLAST or LAST hits
def get_BLAST_LAST_parsed_hits(sample_name, folder_path):
    results = []
    # for the LAST algorithm

    regPattern = re.compile(r'.LASTout.parsed.txt$', flags= re.IGNORECASE)
    input_dir = folder_path +  PATHDELIM + 'blast_results' 
    files = [ re.sub(r'.*\/','',f) for f in glob(input_dir + PATHDELIM + sample_name + '*')  if regPattern.search(f) ]
    regPattern = re.compile(r'[.](.*)[.]LASTout.parsed.txt$', flags=re.IGNORECASE)
    
    for file in files:
      result = regPattern.search(file)
      if result:
         database = result.group(1)
         file_name = input_dir + PATHDELIM + sample_name + '.' +  result.group(1) + '.LASTout.parsed.txt'
         eprintf("\nParse LAST hits for : %s...", database)
         count  =  get_number_of_uncommented_lines(file_name)
         results.append(('Total number of selected hits in ' + database + ' with LAST ', count ) )
         
    # now for the BLAST algorithm
    regPattern = re.compile(r'.BLASTout.parsed.txt')
    input_dir = folder_path +  PATHDELIM + 'blast_results' 
    files = [ re.sub(r'.*\/','',f) for f in glob(input_dir + PATHDELIM + sample_name + '*')  if regPattern.search(f) ]
    regPattern = re.compile(r'[.](.*)[.]BLASTout')
    
    for file in files:
      result = regPattern.search(file)
      if result:
         database = result.group(1)
         file_name = input_dir + PATHDELIM + sample_name + '.' +  result.group(1) + '.BLASTout.parsed.txt'
         eprintf("\nParse BLAST hits for : %s...", database)
         count  =  get_number_of_uncommented_lines(file_name)
         results.append(('Total number of selected hits in ' + database + ' with BLAST ', count ) )

    return results



# counts the number of BLAST or LAST hits
def get_BLAST_LAST_hits(sample_name, folder_path):
    results = []
    # for the LAST algorithm

    regPattern = re.compile(r'.LASTout$')
    input_dir = folder_path +  PATHDELIM + 'blast_results' 
    files = [ re.sub(r'.*\/','',f) for f in glob(input_dir + PATHDELIM + sample_name + '*')  if regPattern.search(f) ]
    regPattern = re.compile(r'[.](.*)[.]LASTout$')
    
    for file in files:
      result = regPattern.search(file)
      if result:
         database = result.group(1)
         file_name = input_dir + PATHDELIM + sample_name + '.' +  result.group(1) + '.LASTout'
         eprintf("\nProcess LAST hits for : %s...", database)
         count  =  get_number_of_BLAST_LAST_hits(file_name)
         eprintf("done")
         results.append(('Total number of hits in ' + database + ' with LAST ', count ) )
         
    # now for the BLAST algorithm
    regPattern = re.compile(r'.BLASTout')
    input_dir = folder_path +  PATHDELIM + 'blast_results' 
    files = [ re.sub(r'.*\/','',f) for f in glob(input_dir + PATHDELIM + sample_name + '*')  if regPattern.search(f) ]
    regPattern = re.compile(r'[.](.*)[.]BLASTout')
    
    for file in files:
      result = regPattern.search(file)
      if result:
         database = result.group(1)
         file_name = input_dir + PATHDELIM + sample_name + '.' +  result.group(1) + '.BLASTout'
         eprintf("\nProcess BLAST hits for : %s...", database)
         count  =  get_number_of_BLAST_LAST_hits(file_name)
         results.append( (  'Total number of hits in ' + database + ' with BLAST ', count ) )

    return results



def get_stats_from_stats_file(sample_name, folder_path, type):
    sequencesPATTERN = re.compile(r'\t([^\t]* of sequences[^\t]*)\t([^\t]*)\t([^\t]*)$')
    input_file_name = folder_path + PATHDELIM  + 'run_statistics' + PATHDELIM + sample_name +  '.' + type + '.stats'
    results = []

    if type=='nuc':
       tag = ' (nucleotide) '
    else:
       tag = ' (amino) '

    try:
       inputfilename = open(input_file_name, 'r')
    except:
       return results 

    lines = inputfilename.readlines()
    inputfilename.close()

    for line in lines:
       line = re.sub(r':', '', line)
       result = sequencesPATTERN.search(line.strip())
       if result:
          try:
             num2 = '%d' %(float(result.group(2)) )
          except:
             num2 = 0
          try:
             num3 = '%d' %(float(result.group(3)) )
          except:
             num3 = 0
          results.append( (result.group(1) + tag + 'BEFORE filtering ', num2 ) )
          results.append( (result.group(1) + tag + 'AFTER filtering ',  num3 ) )

    return results
     

def main(argv): 
    (opts, args) = parser.parse_args(argv)

    if not valid_arguments(opts, args):
       print usage
       sys.exit(0)
    
    sample_name = opts.sample_name
    folder_path = opts.folder_path
    results = []

    # read the nucleotide seequences
    status = get_stats_from_stats_file(sample_name, folder_path, 'nuc')
    results += status

    # read the nucleotide seequences
    status = get_stats_from_stats_file(sample_name, folder_path, 'amino')
    results += status

    # read the blast/last hits
    status = get_BLAST_LAST_hits(sample_name, folder_path)
    results += status

    # read the selected parsed blast/last hits
    status = get_BLAST_LAST_parsed_hits(sample_name, folder_path)
    results += status

    # read the annotated gff hits
    status = get_annotation_hits(sample_name, folder_path)
    results += status

    # read the annotated gff hits
    status = get_functional_taxonomic_hits(sample_name, folder_path)
    results += status

    # read the number of ORFs that are used for mapping to functional categories
    status =  get_ORF_annotations_hits(sample_name, folder_path)
    results += status

    # get the rRNA hits
    status = get_rRNA_hits(sample_name, folder_path)
    results += status

    # get the tRNA hits
    status = get_tRNA_hits(sample_name, folder_path)
    results += status

    stats_file_name = folder_path + PATHDELIM + 'run_statistics' + PATHDELIM + sample_name + '.run.stats.txt' 

    try:
       statsfilename = open(stats_file_name, 'w')
    except:
       print "ERRROR : Cannot open stats file format " + stats_file_name 
       sys.exit(0)
      
    for pair in results:
       fprintf(statsfilename, '%s\t%s\n', pair[0], pair[1])
    statsfilename.close()


def MetaPathways_gather_run_stats(argv):
    main(argv) 
    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

