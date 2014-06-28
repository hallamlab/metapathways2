#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     from os import makedirs, sys, remove, path, _exit
     import re
     from optparse import OptionParser, OptionGroup

     from libs.python_modules.taxonomy.LCAComputation import *
     from libs.python_modules.taxonomy.MeganTree import *
     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, fprintf, printf, eprintf,  GffFileParser, exit_process
     from libs.python_modules.utils.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)


usage= """./MetapathWays_annotate.py -d dbname1 -b parsed_blastout_for_database1 [-d dbname2 -b parsed_blastout_for_database2 ] --input-annotated-gff input.gff  """
parser=None
def createParser():
     global parser
     parser = OptionParser(usage)
     parser.add_option("-b", "--blastoutput", dest="input_blastout", action='append', default=[],
                       help='blastout files in TSV format [at least 1 REQUIRED]')
     
     parser.add_option("-d", "--dbasename", dest="database_name", action='append', default=[],
                       help='the database names [at least 1 REQUIRED]')
     
     cutoffs_group =  OptionGroup(parser, 'Cuttoff Related Options')
     
     cutoffs_group.add_option("--min_score", dest="min_score", type='float', default=20,
                       help='the minimum bit score cutoff [default = 20 ] ')
     
     cutoffs_group.add_option("--max_evalue", dest="max_evalue", type='float', default=1e-6,
                       help='the maximum E-value cutoff [ default = 1e-6 ] ')
     cutoffs_group.add_option("--min_length", dest="min_length", type='float', default=30,
                       help='the minimum length of query cutoff [default = 30 ] ')
     cutoffs_group.add_option("--max_length", dest="max_length", type='float', default=10000,
                       help='the maximum length of query cutoff [default = 10000 ] ')
     
     
     cutoffs_group.add_option("--min_identity", dest="min_identity", type='float', default=20,
                       help='the minimum identity of query cutoff [default 30 ] ')
     cutoffs_group.add_option("--max_identity", dest="max_identity", type='float', default=100,
                       help='the maximum identity of query cutoff [default = 100 ] ')
     
     cutoffs_group.add_option("--limit", dest="limit", type='float', default=5,
                       help='max number of hits per query cutoff [default = 5 ] ')
     
     cutoffs_group.add_option("--min_bsr", dest="min_bsr", type='float', default=0.0,
                       help='minimum BIT SCORE RATIO [default = 0.30 ] ')
     parser.add_option_group(cutoffs_group)
     
     
     output_options_group =  OptionGroup(parser, 'Output table Options')
     output_options_group.add_option("--ncbi-taxonomy-map", dest="ncbi_taxonomy_map",  default=False,
                       help='add the ncbi taxonomy map ')
     
     output_options_group.add_option( "--input-cog-maps", dest="input_cog_maps",
                      help='input cog maps file')
     
     output_options_group.add_option( "--subsystems2peg-file", dest="subsystems2peg_file", default = False,
                      help='the subsystems to peg file from fpt.theseed.org')
     
     output_options_group.add_option( "--input-kegg-maps", dest="input_kegg_maps",
                      help='input kegg maps file')
     
     output_options_group.add_option( "--input-seed-maps", dest="input_seed_maps",
                      help='input seed maps file')
     
     output_options_group.add_option('--input-annotated-gff', dest='input_annotated_gff',
                     metavar='INPUT', help='Annotated gff file [REQUIRED]')
     
     output_options_group.add_option('--output-dir', dest='output_dir',
                     metavar='INPUT', help='Output directory [REQUIRED]')

     output_options_group.add_option("-v", "--verbose",
                     action="store_true", dest="verbose", default=False,
                     help="print lots of information to the stdout [default Off]")

     parser.add_option_group(output_options_group)
     
     lca_options_group =  OptionGroup(parser, 'LCA algorithm Options')
     lca_options_group.add_option("--lca-min-score", dest="lca_min_score",  type='float', default=50,
                       help='minimum BLAST/LAST score to consider as for LCA rule')
     lca_options_group.add_option("--lca-top-percent", dest="lca_top_percent",  type='float', default=10,
                       help='set of considered matches are within this percent of the highest score hit')
     lca_options_group.add_option("--lca-min-support", dest="lca_min_support",  type='int', default=2,
                       help='minimum number of reads that must be assigned to a taxon for ' +\
                            'that taxon to be present otherwise move up the tree until there ' + 
                            'is a taxon that meets the requirement')
     parser.add_option_group(lca_options_group)


def printlist(list, lim):
    i = 0;
    for item in list:
       print item
       i += 1
       if i > lim:
           break

def check_arguments(opts, args):
    return True

    if len(opts.input_blastout) == 0:
         print "There sould be at least one blastoutput file"  
         return False

    if len(opts.database_name) == 0:
         print "There sould be at least one database name"  
         return False

    if len(opts.input_blastout) != len(opts.database_name)  :
         print "The number of database names, blastoutputs files should be equal"
         return False

    if opts.input_annotated_gff == None:
       print "Must specify the input annotated gff file"
       return False

    if opts.output_dir == None:
       print "Must specify the output dir"
       return False

    return True

def process_gff_file(gff_file_name, orf_dictionary):
     try:
        gfffile = open(gff_file_name, 'r')
     except IOError:
        print "Cannot read file " + gff_file_name + " !"

     gff_lines = gfffile.readlines()
     gff_beg_pattern = re.compile("^#")
     gfffile.close()
     
     for line in gff_lines:
        line = line.strip() 
        if gff_beg_pattern.search(line):
          continue
        insert_orf_into_dict(line, orf_dictionary)


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
           

def copyList(a, b):
    [ b.append(x) for x in a ]

    
def get_species(hit):
    if not 'product' in hit: 
        return None

    species = []
    try:
        m = re.findall(r'\[([^\[]+)\]', hit['product'])
        if m != None:
          copyList(m,species)
    except:
          return None

    if species:
       return species
    else:
       return None


def create_annotation(results_dictionary, dbname,  annotated_gff,  output_dir, Taxons, orfsPicked, orfToContig, lca):

    meganTree = None
    #lca.set_results_dictionary(results_dictionary)
    if not path.exists(output_dir):
       makedirs(output_dir)

    orf_dictionary={}
    #process_gff_file(annotated_gff, orf_dictionary)
    gffreader = GffFileParser(annotated_gff)
    output_table_file = open(output_dir + '/functional_and_taxonomic_table.txt', 'a')

    count = 0
    for contig in  gffreader:
       for orf in  gffreader.orf_dictionary[contig]:
          if orf['id'] not in orfsPicked:
            continue
          
          orfToContig[orf['id']] = contig
          
          taxonomy = None
          # if count%10000==0 :
          #    pass 

          #_results = re.search(r'refseq', opts_global.database_name, re.I)
          if orf['id'] in Taxons:
              taxonomy1=Taxons[orf['id']]
              taxonomy=lca.get_supported_taxon(taxonomy1)
          else:
              taxonomy = 'root'

          fprintf(output_table_file, "%s", orf['id'])
          fprintf(output_table_file, "\t%s", orf['orf_length'])
          fprintf(output_table_file, "\t%s", orf['start'])
          fprintf(output_table_file, "\t%s", orf['end'])
          fprintf(output_table_file, "\t%s", orf['seqname'])
          fprintf(output_table_file, "\t%s", orf['contig_length'])
          fprintf(output_table_file, "\t%s", orf['strand'])
          fprintf(output_table_file, "\t%s", orf['ec'])
          # fprintf(output_table_file, "\t%s", str(species))
          fprintf(output_table_file, "\t%s", taxonomy)
          fprintf(output_table_file, "\t%s\n", orf['product'])

          # adding taxons to the megan tree
          #if meganTree and taxonomy != '':
          #    meganTree.insertTaxon(taxonomy)
          #print meganTree.getChildToParentMap()

    output_table_file.close()

    # this prints out the megan tree
#    if meganTree:
#        megan_tree_file = open(output_dir + '/megan_tree.tre', 'w')
#        fprintf(megan_tree_file,  "%s;", meganTree.printTree('1'))
#        megan_tree_file.close()
    

            #write_annotation_for_orf(outputgff_file, candidatedbname, dbname_weight, results_dictionary, orf_dictionary, contig, candidate_orf_pos,  orf['id']) 


def remove_repeats(filtered_words):
    word_dict = {}
    newlist = []
    for word in filtered_words:
       if not word in word_dict:
          if not word in ['', 'is', 'have', 'has', 'will', 'can', 'should',  'in', 'at', 'upon', 'the', 'a', 'an', 'on', 'for', 'of', 'by', 'with' ,'and',  '>' ]:
             word_dict[word]=1
             newlist.append(word)
    return ' '.join(newlist)


class BlastOutputTsvParser(object):

    def __init__(self, dbname,  blastoutput):
        self.lineToProcess = ""
        self.dbname = dbname
        self.blastoutput = blastoutput
        self.i=0
        self.SIZE = 10000
        self.data = {}
        self.fieldmap={}
        self.seq_beg_pattern = re.compile("^#")
        self.lines = []
        self.headerline = None

        self.MAX_READ_ERRORS_ALLOWED = 0
        self.ERROR_COUNT = 0
        self.STEP_NAME = 'CREATE_REPORT_FILES' #PARSE_BLAST'
        self.error_and_warning_logger = None

        try:
           self.blastoutputfile = open( blastoutput,'r')
           line = self.blastoutputfile.readline()
           if not self.seq_beg_pattern.search(line) :
              eprintf("First line must have field header names and begin with \"#\"\n")
              exit_process()

           self.headerline = line.strip()
           self.lineToProcess = self.headerline
           header = re.sub('^#','',line)
           fields = [ x.strip()  for x in header.rstrip().split('\t')]
           k = 0 
           for x in fields:
             self.fieldmap[x] = k 
             k += 1

        except AttributeError:
           print "Cannot read the map file for database :" + dbname
           sys.exit(0)

    def setMaxErrorsLimit(self, max):
       self.MAX_READ_ERRORS_ALLOWED = max

    def setErrorAndWarningLogger(self, logger):
       self.error_and_warning_logger = logger

    def setSTEP_NAME(self, step_name):
        self.STEP_NAME  = step_name

    def getHeaderLine(self):
       return self.headerline

    def getProcessedLine(self):
       return self.lineToProcess

    def refillBuffer(self):
       i = 0 
       self.lines = []
       while i < self.SIZE:
         line = self.blastoutputfile.readline().strip()
         if not line:
           break
         if self.seq_beg_pattern.search(line):
            continue

         self.lines.append(line)
         i += 1
       self.size = len(self.lines)
 
    def rewind(self): 
        self.i = self.i - 1

    def __iter__(self):
        return self
 
    def next(self):
        if self.i % self.SIZE == 0:
           self.refillBuffer()
           if len(self.lines)==0:
              raise StopIteration()

        if self.i % self.SIZE < self.size:
           fields = [ x.strip()  for x in self.lines[self.i % self.SIZE].split('\t')]
           try:
              self.data = {}
              self.data['query'] = fields[self.fieldmap['query']]
              self.data['q_length'] = int(fields[self.fieldmap['q_length']])
              self.data['bitscore'] = float(fields[self.fieldmap['bitscore']])
              self.data['bsr'] = float(fields[self.fieldmap['bsr']])
              self.data['target'] = fields[self.fieldmap['target']]
              self.data['aln_length'] = float(fields[self.fieldmap['aln_length']])
              self.data['expect'] = float(fields[self.fieldmap['expect']])
              self.data['identity'] = float(fields[self.fieldmap['identity']])
              self.data['ec'] = fields[self.fieldmap['ec']]
              self.data['product'] = re.sub(r'=',' ',fields[self.fieldmap['product']])
              self.lineToProcess = self.lines[self.i % self.SIZE]
           except:
              self.ERROR_COUNT += 1
              if self.MAX_READ_ERRORS_ALLOWED > self.ERROR_COUNT:
                 eprintf("%s\tWARNING\till-formatted line \"%s\" \t %s\n", self.STEP_NAME,  self.lines[self.i % self.SIZE], self.blastoutput)
                 if self.error_and_warning_logger != None:
                     self.error_and_warning_logger.write("%s\tWARNING\till-formatted line :\"%s\" \t source : %s\n" %(self.STEP_NAME,  re.sub(r'\t', '<tab>', self.lines[self.i % self.SIZE]) , self.blastoutput))
                 self.i = self.i + 1
                 self.next()
              else:
                 if self.error_and_warning_logger != None:
                     self.error_and_warning_logger.write("%s\tERROR\tThe number of lines in file %s exceeded the max tolerance %d\n" %(self.blastoutput,  self.MAX_READ_ERRORS_ALLOWED) )
                 exit_process() 


#              print "<<<<<<-------"
#              print 'self size ' + str(self.size)
#              print 'line ' + self.lines[self.i % self.SIZE]
#              print 'num fields ' + str(len(fields))
#              fields = [ x  for x in self.lines[self.i % self.SIZE].split('\t')]
#              for field in fields:
#                 print field
#              print 'next line ' + self.lines[(self.i + 1) % self.SIZE]
#              print ' field map ' + str(self.fieldmap)
#              print 'index ' + str(self.i)
#              print 'data ' + str(self.data)
#              print 'fields ' + str(fields)
#              print ' while processing file ' + self.blastoutput
#              print ">>>>>>-------"
#              import traceback 
#              print traceback.print_exc()
           
           self.i = self.i + 1
           return self.data
        else:
           self.lineToProcess = None
           self.blastoutputfile.close()
           raise StopIteration()
              
def isWithinCutoffs(data, cutoffs):
  import traceback

  try:
    if data['q_length'] < cutoffs.min_length:
       return False

    if data['bitscore'] < cutoffs.min_score:
       return False

    if data['expect'] > cutoffs.max_evalue:
       return False

    if data['identity'] < cutoffs.min_identity:
       return False

    if data['bsr'] < cutoffs.min_bsr:
       return False
  except:
     print traceback.print_exc()
    # print cutoffs
     sys.exit(0)

  return True


def process_parsed_blastoutput(dbname, blastparser, cutoffs, annotation_results, pickorfs):
    fields = ['target', 'q_length', 'bitscore', 'bsr', 'expect', 'identity', 'ec', 'query' ]
    fields.append('product')
    for data in blastparser:
        if  data!=None and isWithinCutoffs(data, cutoffs) :
           # if dbname=='refseq':
           # print data['query'] + '\t' + str(data['q_length']) +'\t' + str(data['bitscore']) +'\t' + str(data['expect']) +'\t' + str(data['identity']) + '\t' + str(data['bsr']) + '\t' + data['ec'] + '\t' + data['product']
           annotation = {}
           for field in fields:
             if field in data:
                annotation[field] = data[field] 

           if not data['query'] in pickorfs: 
              blastparser.rewind()
              return None
          
           annotation['dbname'] = dbname

           if not data['query'] in annotation_results:
               annotation_results[data['query']] = []

           annotation_results[data['query']].append(annotation)

    return None

def beginning_valid_field(line):
    fields = [ x.strip() for x in line.split('\t') ]
    count =0
    for field in fields:
       if len(field) > 0:
         return count
       count+=1

    return -1

# creates an empty hierarchical tree with zeros at the lowest count
def read_map_file(dbname_map_filename, field_to_description, hierarchical_map) :
    try:
       map_file = open(dbname_map_filename, 'r')

       map_filelines = map_file.readlines()
    except:
       eprintf("ERROR: Cannot open file %s\n", dbname_map_filename)
       exit_process()
    
    tempfields = [ '', '', '', '', '', '', '' ]
    for line in map_filelines:
       pos = beginning_valid_field(line)
       if pos==-1: 
          continue

       fields = [ x.strip() for x in line.split('\t') ]
       
       tempfields[pos] = fields[pos]
       if len(fields) > pos + 1:
          field_to_description[fields[pos]] = fields[pos+1]
       else:
          field_to_description[fields[pos]] = fields[pos]
       
       i=0
       temp_hierarchical_map = hierarchical_map
       while i < pos :
          temp_hierarchical_map = temp_hierarchical_map[ tempfields[i] ]
          i+=1
     
       temp_hierarchical_map[ tempfields[i] ] = {}
    fill_hierarchy_with_zeroes(hierarchical_map)



def fill_hierarchy_with_zeroes(dictionary):
     for key in dictionary.keys():
        if len(dictionary[key]) ==0 :
           dictionary[key] = 0
        else:
           fill_hierarchy_with_zeroes(dictionary[key])


def cog_id(product):
    results = re.search(r'COG[0-9][0-9][0-9][0-9]', product)
    cog_id = ''
    if results:
       cog_id=results.group(0)
    return cog_id
    

def kegg_id(product):
    results = re.search(r'K[0-9][0-9][0-9][0-9][0-9]', product)
    kegg_id = ''
    if results:
       kegg_id=results.group(0)
    return kegg_id

def create_table(results,  std_dbname, output_dir, hierarchical_map, field_to_description):
    if not path.exists(output_dir):
       makedirs(output_dir)

    #print field_to_description
    orthology_count = {}
    for key in field_to_description[std_dbname]:
       orthology_count[key] = 0 
    
    #print hierarchical_map 
    for seqname in results:    
       for orf in results[seqname]:
           if std_dbname =='seed': 
              seed = re.sub(r'\[.*\]','', orf['product']).strip()
              if seed in orthology_count:
                 orthology_count[seed]+=1

           if std_dbname =='cog': 
              cog =  cog_id(orf['product'])
              if cog in orthology_count:
                 orthology_count[cog]+=1

           if std_dbname =='kegg': 
              kegg =  kegg_id(orf['product'])
              if kegg in orthology_count:
                 orthology_count[kegg]+=1

   # print orthology_count.keys()
    add_counts_to_hierarchical_map(hierarchical_map[std_dbname], orthology_count)
    

def print_kegg_cog_tables(dbname, output_dir, hierarchical_map, field_to_description, filePermType = 'w'):

    if dbname=='cog':
       outputfile = open( output_dir +'/COG_stats_1.txt', filePermType)
       print_counts_at_level(hierarchical_map[dbname], field_to_description[dbname],  0, 0, outputfile, printKey=False,\
          header="Functional Category\tGene Count") 
       outputfile.close()
       outputfile = open( output_dir +'/COG_stats_2.txt', filePermType)
       print_counts_at_level(hierarchical_map[dbname], field_to_description[dbname],  0, 1, outputfile,\
          header="Function Abbr\tFunctional Category\tGene Count") 
       outputfile.close()
       outputfile = open( output_dir +'/COG_stats_3.txt', filePermType)
       print_counts_at_level(hierarchical_map[dbname], field_to_description[dbname],  0, 2, outputfile,\
          header="COGID\tFunction\tGene Count") 
       outputfile.close()

    if dbname=='kegg':
       outputfile = open( output_dir +'/KEGG_stats_1.txt', filePermType)
       print_counts_at_level(hierarchical_map[dbname], field_to_description[dbname],  0, 0, outputfile, printKey=False,\
          header="Function Category Level 1\tGene Count") 
       outputfile.close()
       outputfile = open( output_dir +'/KEGG_stats_2.txt', filePermType)
       print_counts_at_level(hierarchical_map[dbname], field_to_description[dbname],  0, 1, outputfile, printKey=False,\
          header="Function Category Level 2a\tGene Count") 
       outputfile.close()
       outputfile = open( output_dir +'/KEGG_stats_3.txt', filePermType)
       print_counts_at_level(hierarchical_map[dbname], field_to_description[dbname],  0, 2, outputfile,\
         header="ID\tFunction Category Level 3\tGene Count" ) 
       outputfile.close()
       outputfile = open( output_dir +'/KEGG_stats_4.txt', filePermType)
       print_counts_at_level(hierarchical_map[dbname], field_to_description[dbname],  0, 3, outputfile,\
         header="KO\tFunction Category Level 4\tGene Count") 
       outputfile.close()


def print_counts_at_level(hierarchical_map, field_to_description,  depth, level, outputfile, printKey=True, header=None): 
    
    if type(hierarchical_map) is type(0):
       return hierarchical_map
    if header:
       fprintf(outputfile, "%s\n",header )

    count = 0
    for key in hierarchical_map:  
       tempcount = print_counts_at_level(hierarchical_map[key],field_to_description, depth+1, level, outputfile, printKey=printKey)
       if depth==level:
          if key in field_to_description:
              if printKey:
                 fprintf(outputfile, "%s\n", key + '\t' + field_to_description[key] + '\t' +  str(tempcount) )
              else:
                 fprintf(outputfile, "%s\n",  field_to_description[key] + '\t' +  str(tempcount) )
          else:
              if printKey:
                 fprintf(outputfile, "%s\n", key + '\t' + ' ' + '\t' + str(tempcount))
              else:
                 fprintf(outputfile, "%s\n", key +  '\t' + str(tempcount))
       count+=tempcount
    return count


# this function adds the count into the hierarchical map at the lowest level
def  add_counts_to_hierarchical_map(hierarchical_map, orthology_count):
     
  try:
    for key in hierarchical_map:  
      if type(hierarchical_map[key])==int:  
         if key in orthology_count:
             hierarchical_map[key]+=int(orthology_count[key])
      else:
          add_counts_to_hierarchical_map(hierarchical_map[key], orthology_count)
  except:

      import traceback
      traceback.print_exc()
      print len(hierarchical_map[key])
      sys.exit(0)
   
def get_list_of_queries(annotated_gff):
    orfList = {}
    gffreader = GffFileParser(annotated_gff)
    count = 0
    for contig in  gffreader:
       for orf in  gffreader.orf_dictionary[contig]:
          orfList[orf['id']]  = 1
          count += 1
    #      if count%500000==0:
    #         print count
            
    return orfList.keys()


def Heapify(A, i, S):
    while True:
       l = 2*i + 1
       r = l + 1
       max = i
       if l < S and A[l][1] < A[max][1]:  # was > 
          max = l
       if r < S and A[r][1] < A[max][1]:  # was >
          max = r

       if max != i and i < S:
          temp = A[i]
          A[i] =  A[max]
          A[max] =  temp
       else:
          break
       i = max

def BuildHeap(S, A):
    i = int(S/2)
    while  i >= 0:
        Heapify(A, i, S)
        i = i - 1

def writeParsedLines(fieldmapHeaderline, parsedLines, list, names, outputfilename):
    try:
      outputfile = open(outputfilename, 'w')
    except OSError:
      print "ERROR: Cannot create sequence file : " + faa_file
      sys.exit(0)

    outputStr=fieldmapHeaderline + "\n"
    fprintf(outputfile, "%s", outputStr)

    outputStr=""
    i = 0
    for item in list:
       outputStr += parsedLines[item[0]]+'\n'
       if i% 1000==0 and i > 0:
          fprintf(outputfile, "%s", outputStr)
          outputStr=""
       i += 1

    if len(outputStr) > 0:
      fprintf(outputfile, "%s", outputStr)

    outputfile.close()

def merge_sorted_parsed_files(dbname, filenames, outputfilename, orfRanks, verbose=False, errorlogger = None):
    linecount = 0
    readerhandles = []

    if verbose:
       eprintf("Processing for database  : %s\n", dbname)

    if len(filenames)==0:
       eprintf("WARNING : Cannot find any B/LAST output file for database : %\n", dbname)
       exit_process()
    
    try:
       for i in range(len(filenames)):
         #print filenames
         readerhandles.append(BlastOutputTsvParser(dbname, filenames[i]) )
    except OSError:
      eprintf("ERROR: Cannot read sequence file : %s\n", filenames[i])
      exit_process()

    # set error and warning parameters 
    for readerhandle in readerhandles:
        readerhandle.setMaxErrorsLimit(5)
        readerhandle.setErrorAndWarningLogger(errorlogger)
        readerhandle.setSTEP_NAME('PARSE BLAST')
    
    try:
       outputfile = open(outputfilename, 'w')
       fieldmapHeaderLine = readerhandles[0].getHeaderLine()
       fprintf(outputfile, "%s\n",fieldmapHeaderLine) 
    except OSError: 
       eprintf("ERROR: Cannot create sequence file : %s\n", outputfilename)
       exit_process()

    values = []    
    for i in range(len(filenames)):
       iterate = iter(readerhandles[i])
       try :
          next(iterate)
          line = readerhandles[i].getProcessedLine()  
          fields  = [ x.strip() for x in line.split('\t') ]
          values.append( (i, orfRanks[fields[0]], line) )
       except:
          outputfile.close()
          return

    S = len(filenames) 
    BuildHeap(S, values)
    
    while S>0:
       try:
          iterate = iter(readerhandles[values[0][0]])
          line = readerhandles[values[0][0]].getProcessedLine()  
          fields  = [ x.strip() for x in line.split('\t') ]
          #print fields[0], orfRanks[fields[0]] 
          fprintf(outputfile, "%s\n",line) 
          next(iterate)

          line = readerhandles[values[0][0]].getProcessedLine()  
          fields  = [ x.strip() for x in line.split('\t') ]
          values[0] = (values[0][0], orfRanks[fields[0]], line) 
       except:
          #import traceback
          #traceback.print_exc()
          #print 'finished ' + str(S)
          values[0] = values[S-1] 
          S = S - 1

       if S>0:
          Heapify(values, 0, S)

    #print 'line count ' + str(linecount)
    outputfile.close()


def  create_sorted_parse_blast_files(dbname, blastoutput, listOfOrfs, size = 100000, verbose=False, errorlogger= None):
    orfRanks = {}
    count = 0
    for orf in listOfOrfs:
       orfRanks[orf] = count
       count += 1

    sorted_parse_file =  blastoutput + ".tmp"

    currSize = 0 
    parsedLines={}
    list = []
    names = {}
    seqid =0
    batch = 0 
    filenames = []

    if verbose:
       eprintf("\n\n\n")
       eprintf("dbname                   : %s\n", dbname)
       eprintf("Parsed file              : %s\n",  blastoutput)

    blastparser =  BlastOutputTsvParser(dbname, blastoutput)
    blastparser.setMaxErrorsLimit(5)
    blastparser.setErrorAndWarningLogger(errorlogger)

    fieldmapHeaderLine = blastparser.getHeaderLine()
    
    for data in blastparser:
       names[seqid] = data['query']
       parsedLines[seqid] = blastparser.getProcessedLine()
       list.append( (seqid, names[seqid]) )
         
       seqid +=1 
       currSize += 1
       if currSize % size ==0: 
           list.sort(key=lambda tup: tup[1], reverse=False)
           #print "Num of lines writing to file " + sorted_parse_file + "." + str(batch) + "  :  " + str(len(list))
           writeParsedLines(fieldmapHeaderLine, parsedLines, list, names, sorted_parse_file + "." + str(batch))
           filenames.append(sorted_parse_file + "." + str(batch))
           batch += 1
           list = []
           names = {}
           seqid =0
           parsedLines = {}

    
    if currSize == 0:
        list.sort(key=lambda tup: tup[1], reverse=False)
        writeParsedLines(fieldmapHeaderLine, parsedLines, list, names, sorted_parse_file + "." + str(batch))
        filenames.append(sorted_parse_file + "." + str(batch))
    else :
       if currSize % size !=0 : 
          list.sort(key=lambda tup: tup[1], reverse=False)
          writeParsedLines(fieldmapHeaderLine, parsedLines, list, names, sorted_parse_file + "." + str(batch))
          filenames.append(sorted_parse_file + "." + str(batch))


    if verbose:
       eprintf("Number of lines          : %s\n", str(currSize))
       eprintf("Merge the files          : %s\n" , sorted_parse_file)
       eprintf("Number of files to merge : %s\n", str(len(filenames)))

    merge_sorted_parsed_files(dbname, filenames, sorted_parse_file, orfRanks, verbose = verbose, errorlogger = errorlogger)

    #remove the split files
    for file in filenames:
       remove(file)


import gc
import resource

opts_global = ""

# the main function
def main(argv, errorlogger = None,  runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)
    global opts_global
    opts_global = opts
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)


    db_to_map_Maps =  {'cog':opts.input_cog_maps, 'seed':opts.input_seed_maps, 'kegg':opts.input_kegg_maps}


    results_dictionary={}
    dbname_weight={}

#    print "memory used  = %s" %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss /1000000)
    listOfOrfs =  get_list_of_queries(opts.input_annotated_gff)       
    listOfOrfs.sort(key=lambda tup: tup, reverse=False)
    
    #printlist(listOfOrfs,5)
    #sys.exit(0)

##### uncomment the following lines
    for dbname, blastoutput in zip(opts.database_name, opts.input_blastout):
      create_sorted_parse_blast_files(dbname, blastoutput, listOfOrfs, verbose= opts.verbose, errorlogger = errorlogger) 
#####

    # process in blocks of size _stride
    output_table_file = open(opts.output_dir + '/functional_and_taxonomic_table.txt', 'w')
    fprintf(output_table_file, "ORF_ID\tORF_length\tstart\tend\tContig_Name\tContig_length\tstrand\tec\ttaxonomy\tproduct\n")
    output_table_file.close()

    lca = LCAComputation(opts.ncbi_taxonomy_map)
    lca.setParameters(opts.lca_min_score, opts.lca_top_percent, opts.lca_min_support)

    blastParsers={}
    for dbname, blastoutput in zip( opts.database_name, opts.input_blastout):
        blastParsers[dbname] =  BlastOutputTsvParser(dbname, blastoutput + '.tmp')
        blastParsers[dbname].setMaxErrorsLimit(5)
        blastParsers[dbname].setErrorAndWarningLogger(errorlogger)

    #this part of the code computes the occurence of each of the taxons
    # which is use in the later stage is used to evaluate the min support
    # as used in the MEGAN software

    start = 0
    Length = len(listOfOrfs) 
    _stride = 100000
    Taxons = {}
    while start < Length:
       pickorfs= {}
       last =  min(Length, start + _stride)
       for  i in range(start, last): 
          pickorfs[listOfOrfs[i]]= 'root'
       start = last
       #print 'Num of Min support orfs ' + str(start)

       results_dictionary={}
       for dbname, blastoutput in zip( opts.database_name, opts.input_blastout):
          results = re.search(r'refseq', dbname, re.I)
          if results:
            try:
               results_dictionary[dbname]={}
               process_parsed_blastoutput(dbname, blastParsers[dbname], opts, results_dictionary[dbname], pickorfs)
               lca.set_results_dictionary(results_dictionary)
               lca.compute_min_support_tree(opts.input_annotated_gff, pickorfs, dbname = dbname )
               for key, taxon  in pickorfs.iteritems():
                   Taxons[key] = taxon
            except:
               eprintf("ERROR: while training for min support tree %s\n", dbname)
               import traceback
               traceback.print_exc()

    blastParsers={}
    for dbname, blastoutput in zip( opts.database_name, opts.input_blastout):
        blastParsers[dbname] =  BlastOutputTsvParser(dbname, blastoutput + '.tmp')

    # this loop determines the actual/final taxonomy of each of the ORFs 
    # taking into consideration the min support
    filePermTypes= {}
    start = 0
    outputfile = open( opts.output_dir +'/ORF_annotation_table.txt', 'w')

    
    short_to_long_dbnames = {}
    for dbname in opts.database_name: 
      results = re.search(r'^seed', dbname,  re.IGNORECASE)
      if results:
          short_to_long_dbnames['seed'] = dbname

      results = re.search(r'^cog', dbname,  re.IGNORECASE)
      if results:
          short_to_long_dbnames['cog'] = dbname

      results = re.search(r'^kegg', dbname, re.IGNORECASE)
      if results:
          short_to_long_dbnames['kegg'] = dbname

    standard_dbs = ['cog', 'seed', 'kegg']
    standard_db_maps = [opts.input_cog_maps, opts.input_seed_maps, opts.input_kegg_maps]
    field_to_description = {}
    hierarchical_map = {}

    for db in standard_dbs: 
      if db in short_to_long_dbnames:
        field_to_description[db] = {}
        hierarchical_map[db] = {}

    for dbname in standard_dbs:
       if dbname in short_to_long_dbnames:
          try:
            read_map_file(db_to_map_Maps[dbname], field_to_description[dbname], hierarchical_map[dbname])
          except:
            raise 
            pass

    while start < Length:
       pickorfs= {}
       last =  min(Length, start + _stride)
       for  i in range(start, last): 
          pickorfs[listOfOrfs[i]]= True
       start = last
       gc.collect()
       eprintf("\nMemory used  = %s MB\n", str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000000))
       results_dictionary={}
       for dbname, blastoutput in zip( opts.database_name, opts.input_blastout):
            try:
               results_dictionary[dbname]={}
               eprintf("Processing database %s...", dbname)
               process_parsed_blastoutput(dbname, blastParsers[dbname], opts, results_dictionary[dbname], pickorfs)
               eprintf("done\n")
            except:
               import traceback
               traceback.print_exc()
               eprintf("ERROR: %s\n", dbname)
               pass
            print dbname + ' ' + str(len(results_dictionary[dbname]))

       eprintf("Num orfs processed  : %s\n", str(start))

       # create the annotations now
       orfToContig = {}
       create_annotation(results_dictionary, opts.database_name,  opts.input_annotated_gff, opts.output_dir, Taxons, pickorfs, orfToContig, lca)
        
       for std_dbname, db_map_filename in zip(standard_dbs, standard_db_maps):
         if std_dbname in short_to_long_dbnames:
              create_table(results_dictionary[short_to_long_dbnames[std_dbname]], std_dbname,  opts.output_dir, hierarchical_map, field_to_description)

#             create_table(results_dictionary[dbname], opts.input_kegg_maps, 'kegg', opts.output_dir, filePermType)

       print_orf_table(results_dictionary, orfToContig, opts.output_dir, outputfile)

    for std_dbname, db_map_filename in zip(standard_dbs, standard_db_maps):
       if std_dbname in short_to_long_dbnames:
          print_kegg_cog_tables(std_dbname, opts.output_dir, hierarchical_map, field_to_description,  filePermType = 'w')

    outputfile.close() 
    # now remove the temporary files
    for dbname, blastoutput in zip( opts.database_name, opts.input_blastout):
        try:
           remove( blastoutput + '.tmp')
        except:
           pass

def refseq_id(product):
    results = re.search(r'gi\|[0-9.]*', product)
    refseq_id = ''
    if results:
       refseq_id=results.group(0)
    return refseq_id

def process_subsys2peg_file(subsystems2peg, subsystems2peg_file):
     try:
         orgfile = open(subsystems2peg_file,'r')
     except IOError:
         print "Cannot open " + str(org_file)
     lines = orgfile.readlines()
     orgfile.close()
     for line in lines:
        hits = line.split('\t')
        if len(hits) > 2:
           subsystems2peg[hits[2]]=hits[1]
     try:
        orgfile.close()
     except:
         print "Cannot close " + str(org_file)
 
def print_orf_table(results, orfToContig,  output_dir,  outputfile):
    if not path.exists(output_dir):
       makedirs(output_dir)
   

    orf_dict = {}
    for dbname in results.keys():
      for orfname in results[dbname]:
         for orf in results[dbname][orfname]:
           if not orf['query'] in orf_dict:
               orf_dict[orf['query']] = {}

           _results = re.search(r'cog', dbname, re.I)
           if _results:
              orf_dict[orf['query']][dbname] = cog_id(orf['product'])

           _results = re.search(r'kegg', dbname, re.I)
           if _results:
              orf_dict[orf['query']][dbname] =  kegg_id(orf['product'])

           _results = re.search(r'seed', dbname, re.I)
           if _results:
              orf_dict[orf['query']][dbname] = re.sub(r'\[.*\]','', orf['product']).strip()
             
           _results = re.search(r'metacyc', dbname, re.I)
           if _results:
              orf_dict[orf['query']][dbname] =  orf['product']

           orf_dict[orf['query']]['contig'] = orfToContig[orfname]

    # compute the databases
    database_maps = {}
    for dbname in results.keys():
       _results = re.search(r'cog', dbname, re.I)
       if _results:
         database_maps['cog'] = dbname

       _results = re.search(r'kegg', dbname, re.I)
       if _results:
         database_maps['kegg'] = dbname

       _results = re.search(r'seed', dbname, re.I)
       if _results:
         database_maps['seed'] = dbname

       _results = re.search(r'metacyc', dbname, re.I)
       if _results:
         database_maps['metacyc'] = dbname

    
    for orfn in orf_dict:
#       print orfn, '<<',  orf_dict[orfn], ' >> xxxx'
       #_keys =  orf_dict[orfn].keys()
       #_results = re.search(r'cog', dbname, re.I)

       if 'cog' in database_maps and  database_maps['cog'] in orf_dict[orfn]:
          cogFn = orf_dict[orfn][database_maps['cog']]
       else:
          cogFn = ""

       if 'kegg' in database_maps and database_maps['kegg'] in orf_dict[orfn]:
          keggFn = orf_dict[orfn][database_maps['kegg']]
       else:
          keggFn = ""

       if 'metacyc' in database_maps and database_maps['metacyc'] in orf_dict[orfn]:
          metacycPwy = orf_dict[orfn][database_maps['metacyc']]
       else:
          metacycPwy = ""

       if 'seed' in database_maps and database_maps['seed'] in orf_dict[orfn]:
          seedFn = orf_dict[orfn][database_maps['seed']]
       else:
          seedFn = ""

       fprintf(outputfile, "%s\n", orfn + "\t" + orf_dict[orfn]['contig'] + '\t' + cogFn + '\t' + keggFn +'\t' + seedFn + '\t' + metacycPwy)


def MetaPathways_create_reports_fast(argv, errorlogger =  None, runstatslogger = None):       
    createParser()
    errorlogger.write("#STEP\tPARSE BLAST\n")
    main(argv,errorlogger= errorlogger, runstatslogger = runstatslogger )
    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])
