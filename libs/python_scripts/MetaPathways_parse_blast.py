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
     from os import makedirs, sys, remove, rename
     from sys import path
     import re
     from copy import copy
     
     from optparse import OptionParser, OptionGroup
     from python_modules.metapaths_utils  import parse_command_line_parameters, fprintf, printf
     from python_modules.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)


usage= """./MetapathWays_parse_blast.py -d dbname1 -b blastout_for_database1 -m map_for_database1 [-d dbname2 -b blastout_for_database2 -m map_for_database2 ] """


parser = None

def createParser():
    global parser
    parser = OptionParser(usage)
    parser.add_option("-b", "--blastoutput", dest="input_blastout", action='append', default=[],
                      help='the input blastout files [at least 1 REQUIRED]')
    parser.add_option("-d", "--dbasename", dest="database_name", action='append', default=[],
                      help='the database names [at least 1 REQUIRED]')
    
    parser.add_option("-r", "--ref_score", dest="refscore_file", 
                      help='the refscore  table [REQUIRED]')
    
    parser.add_option("-m", "--map_for_database", dest="database_map", action='append', default=[],
                      help='the map file for the database  [at least 1 REQUIRED]')
    
    parser.add_option("-a", "--algorithm", dest="algorithm", choices = ['BLAST', 'LAST'], default = "BLAST",
                       help='the algorithm used for computing homology [DEFAULT: BLAST]')
    
    cutoffs_group =  OptionGroup(parser, 'Cuttoff Related Options')
    
    cutoffs_group.add_option("--min_score", dest="min_score", type='float', default=20,
                      help='the minimum bit score cutoff [default = 20 ] ')
    cutoffs_group.add_option("--min_query_coverage", dest="min_query_coverage", type='float', default=0,
                      help='the minimum bit query_coverage cutoff [default = 0 ] ')
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
    
    cutoffs_group.add_option("--max_gaps", dest="max_gaps", type='float', default=1000,
                      help='the maximum gaps of query cutoff [default = 1000] ')
    cutoffs_group.add_option("--limit", dest="limit", type='float', default=5,
                      help='max number of hits per query cutoff [default = 5 ] ')
    
    cutoffs_group.add_option("--min_bsr", dest="min_bsr", type='float', default=0.30,
                      help='minimum BIT SCORE RATIO [default = 0.30 ] ')
    parser.add_option_group(cutoffs_group)
    
    output_options_group =  OptionGroup(parser, 'Output Options')
    output_options_group.add_option("--tax", dest="taxonomy", action='store_true', default=False,
                      help='add the taxonomy info [useful for refseq] ')
    output_options_group.add_option("--remove_tax", dest="remove_taxonomy", action='store_true', default=False,
                      help='removes the taxonomy from product [useful for refseq] ')
    output_options_group.add_option("--remove_ec", dest="remove_ec", action='store_true', default=False,
                      help='removes the EC number from product [useful for kegg/metacyc] ')
    parser.add_option_group(output_options_group)




def check_arguments(opts, args):
    if len(opts.input_blastout) == 0:
         print "There sould be at least one blastoutput file"  
         return False

    if len(opts.database_name) == 0:
         print "There sould be at least one database name"  
         return False

    if len(opts.database_map) == 0:
         print "There sould be at least one database map file name"  
         return False

    if len(opts.input_blastout) != len(opts.database_name) or len(opts.input_blastout) !=  len(opts.database_map) :
         print "The number of database names, blastoutputs and database map file should be equal"
         return False


    if opts.refscore_file == None:
       print "Must specify the refscore"
       return False

    return True



def create_query_dictionary(blastoutputfile, query_dictionary, algorithm ):
       seq_beg_pattern = re.compile("^#")

       try:
          blastoutfh = open( blastoutputfile,'r')
       except:
          print ""
          print "ERROR : cannot open B/LAST output file " + blastoutputfile + " to parse "
          pass
  
       try:
          for line in blastoutfh:
             if not seq_beg_pattern.search(line):
                 words = line.rstrip().split('\t')
                 if len(words) != 12: 
                     continue
   
                 if algorithm =='BLAST': 
                    query_dictionary[words[0]] = 1
   
                 if algorithm =='LAST': 
                    query_dictionary[words[1]]= 1
          blastoutfh.close()
       except:
          print "ERROR : while reading  B/LAST output file " + blastoutputfile + " to partse "
          print "      : make sure B/LAST ing was done for the particular database"
          pass 

def create_dictionary(databasemapfile, annot_map, query_dictionary):
       if not query_dictionary:
          print "WARNING : empty query dictionary in parse B/LAST"
          return 

       seq_beg_pattern = re.compile(">")
       dbmapfile = open( databasemapfile,'r')

       for line in dbmapfile:
          if seq_beg_pattern.search(line):
              words = line.rstrip().split()
              name = words[0].replace('>','',1)
              
              if not name in query_dictionary: 
                 continue
              words.pop(0)
              if len(words)==0:
                 annotation = 'hypothetical protein'
              else:
                 annotation = ' '.join(words)
              annot_map[name] = annotation
       dbmapfile.close()

       if len(annot_map) ==0:
           sys.exit( "File " + databasemapfile + " seems to be empty!" ) 
        
def create_refscores(refscores_file, refscore_map):
#       print 'in refscores ' + refscores_file
       refscorefile = open(refscores_file,'r')
       lines=refscorefile.readlines()
       refscorefile.close()
       for line in lines:
           words =[ x.strip()  for x in  line.split('\t') ]
           if len(words) == 2:
              try:
                refscore_map[words[0]]= float(words[1])
              except:
                refscore_map[words[0]]= 1

class BlastOutputParser(object):
    commentPATTERN = re.compile(r'^#')

    def __init__(self, dbname,  blastoutput, database_mapfile, refscore_file, opts):
        self.Size = 10000
        self.dbname = dbname
        self.blastoutput = blastoutput
        self.database_mapfile =database_mapfile
        self.refscore_file = refscore_file
        self.annot_map = {} 
        self.i=0
        self.opts = opts
        self.hits_counts = {}
        self.data = {}
        self.refscores = {}

        #print "trying to open blastoutput file " + blastoutput
        try:
           query_dictionary = {}
           create_query_dictionary(self.blastoutput, query_dictionary, self.opts.algorithm) 
           try:
              self.blastoutputfile = open(self.blastoutput,'r')
           except:
              print ""
              print "ERROR : cannot open B/LAST output file " + blastoutput + " to parse "
              print "      : make sure \"B/LAST\"ing was done for the particular database"
              sys.exit(0)

           create_refscores(refscore_file, self.refscores)

#           print "Going to read the dictionary\n"
           create_dictionary(database_mapfile, self.annot_map, query_dictionary)
           query_dictionary = {}
#           print "Doing reading  dictionary\n"
        except AttributeError:
           print "Cannot read the map file for database :" + dbname
           sys.exit(0)
  
    def __iter__(self):
        return self
 
    def permuteForLAST(self, words):
        try :
           temp = copy(words)
           words[0] = temp[6] # query
           words[1] = temp[1] # target
           words[2] = 100.0 # percent id
           words[3] = temp[3]  #aln length
           words[6] = temp[2]
           words[7] = int(temp[2]) + int(temp[3]) - 1
           words[10] = 0.0   # evalue
           words[11] = temp[0]
        except:
           print "ERROR : Invalid B/LAST output file"
           sys.exit(0) 

    def refillBuffer(self):
        i = 0
        self.lines = []
        line = self.blastoutputfile.readline()
        while line and i < self.Size:
          line=self.blastoutputfile.readline()
          if self.commentPATTERN.match(line):
             continue
          self.lines.append(line)
          if not line:
            break
          i += 1
        self.size = len(self.lines)
       
    def next(self):
        if self.i % self.Size ==0:
           self.refillBuffer()

        if  self.i % self.Size < self.size:
           words = [ x.strip()  for x in self.lines[self.i % self.Size].rstrip().split('\t')]

           if len(words) != 12:
               self.i = self.i + 1
               return None
            
           if  self.opts.algorithm =='LAST':
                self.permuteForLAST(words)

           if not words[0] in self.hits_counts:
              self.hits_counts[words[0]] = 0

           if self.hits_counts[words[0]] >= self.opts.limit:
              self.i = self.i + 1
              return None 

           if len(words) != 12 or not isWithinCutoffs(words, self.data, self.opts, self.annot_map, self.refscores):
             self.i = self.i + 1
             return None 

           self.hits_counts[words[0]] += 1
           self.i = self.i + 1
              
           try:
              return self.data
           except:
              return None
        else:
           self.blastoutputfile.close()
           raise StopIteration()
              
def isWithinCutoffs(words, data, cutoffs, annot_map, refscores):
    data['query'] = words[0]

    try:
       data['target'] = words[1]
    except:
       data['target'] = 0

    try:
       data['q_length'] = int(words[7]) - int(words[6]) + 1
    except:
       data['q_length'] = 0

    try:
       data['bitscore'] = int(words[11])
    except:
       data['bitscore'] = 0

    try:
       data['bsr'] = float(words[11])/refscores[words[0]]
    except:
       #print "words 0 " + str(refscores[words[0]])
       #print "words 11 " + str( words[11])
       data['bsr'] = 0

    try:
       data['expect'] = float(words[10])
    except:
       data['expect'] = 0

    try:
       data['aln_length'] = float(words[3])
    except:
       data['aln_length'] = 0

    try:
       data['identity'] = float(words[2])
    except:
       data['identity'] = 0

    try:
       data['product'] = annot_map[words[1]]
    except:
       print("Sequence with name \"" + words[1] + "\" is not present in map file ")
       data['product'] = 'hypothetical protein'

    try:
       m = re.search(r'(\d+[.]\d+[.]\d+[.]\d+)', data['product'])
       if m != None:
         data['ec'] = m.group(0)
       else:
         data['ec'] = ''
    except:
        data['ec'] = ''

    if cutoffs.taxonomy:
       try:
          m = re.search(r'\[([^\[]+)\]', data['product'])
          if m != None:
            data['taxonomy'] = m.group(1)
          else:
            data['taxonomy'] = ''
       except:
            data['taxonomy'] = ''

    
    if cutoffs.remove_taxonomy:
       try:
          data['product'] = re.sub(r'\[([^\[]+)\]','', data['product'])
       except:
          data['product'] = ''

    if cutoffs.remove_ec:
       try:
          data['product'] = re.sub(r'\([Ee][Ce][:]\d+[.]\d+[.]\d+[.]\d+\)', '', data['product'])
          data['product'] = re.sub(r'\[[Ee][Ce][:]\d+[.]\d+[.]\d+[.]\d+\]', '', data['product'])
          data['product'] = re.sub(r'\[[Ee][Ce][:]\d+[.]\d+[.]\d+[.-]\]', '', data['product'])
          data['product'] = re.sub(r'\[[Ee][Ce][:]\d+[.]\d+[.-.-]\]', '', data['product'])
          data['product'] = re.sub(r'\[[Ee][Ce][:]\d+[.-.-.-]\]', '', data['product'])
       except:
          data['product'] = ''


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

#min_length'
#'min_score'
#'max_evalue'
# 'min_identity'
#'limit'
#'max_length'
#'min_query_coverage'
#'max_gaps'
#min_bsr'


    return True

def add_refscore_to_file(blast_table_out, refscore_file, allNames):
    infile = open( blast_table_out,'r')

    refscores = {}
    lines = infile.readlines()
    for line in lines:
       line=line.rstrip()
       fields = line.split('\t')
       if len(fields) != 12:
          print 'Error in the blastout file'
          sys.exit(1)

    for key, value in refscores.iteritems():
       allNames[key] = True
       fprintf(refscore_file, "%s\t%s\n",key, value)

    infile.close()
        

# compute the refscores
def process_blastoutput(dbname, blastoutput,  mapfile, refscore_file, opts):
    blastparser =  BlastOutputParser(dbname, blastoutput, mapfile, refscore_file, opts)

    fields = ['target','q_length', 'bitscore', 'bsr', 'expect', 'aln_length', 'identity', 'ec' ]
    if opts.taxonomy:
       fields.append('taxonomy')
    fields.append('product')

    output_blastoutput_parsed = blastoutput + '.parsed.txt'
    
    # temporary file is used to deal with incomplete processing of the file
    output_blastoutput_parsed_tmp =  output_blastoutput_parsed + ".tmp"
    outputfile = open(output_blastoutput_parsed_tmp, 'w') 

    # write the headers out
    fprintf(outputfile, "#%s",'query')
    for field in fields:
         fprintf(outputfile,"\t%s",field)
    fprintf(outputfile, "\n")

    for data in blastparser:
        if not data:
          continue
        try:
          fprintf(outputfile, "%s",data['query'])
        except:
           print 'data is : ', data, '\n'
           sys.exit()
        for field in fields:
           fprintf(outputfile, "\t%s",data[field])
        fprintf(outputfile, "\n")

    outputfile.close()
    rename(output_blastoutput_parsed_tmp, output_blastoutput_parsed)


    return None

# the main function
def main(argv): 
    global parser
    (opts, args) = parser.parse_args(argv)
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)
    


    
  # input file to blast with itself to commpute refscore
#    infile = open(input_fasta,'r')

    dictionary={}
    for dbname, blastoutput, mapfile in zip( opts.database_name, opts.input_blastout, opts.database_map):
        temp_refscore = ""
        if opts.algorithm == "LAST":
            temp_refscore = opts.refscore_file + ".LAST"
        if opts.algorithm == "BLAST":
            temp_refscore = opts.refscore_file + ".BLAST"
        process_blastoutput( dbname, blastoutput,  mapfile, temp_refscore, opts)

def MetaPathways_parse_blast(argv):       
    createParser()
    main(argv)
    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

