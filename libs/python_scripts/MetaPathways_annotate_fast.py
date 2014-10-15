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
    from os import makedirs, sys, remove, rename
    from sys import path
    import re, traceback
    from optparse import OptionParser, OptionGroup
    from glob import glob
    from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, fprintf, printf, eprintf
    from libs.python_modules.utils.metapathways_utils  import strip_taxonomy
    from libs.python_modules.utils.sysutil import getstatusoutput

except:
    print """ Could not load some user defined  module functions"""
    print """ Make sure your typed 'source MetaPathwaysrc'"""
    print """ """
    sys.exit(3)


usage=  sys.argv[0] + """ -d dbname1 -b parsed_blastout_for_database1 -w weight_for_database1 [-d dbname2 -b parsed_blastout_for_database2 -w weight_for_database2 ] [ --rRNA_16S  16SrRNA-stats-table ] [ --tRNA tRNA-stats-table ]"""
parser = None

def createParser():
     global parser
     epilog = """Based on the parsed homology search results against reference protein databases, 
               tRNA scan results and matches against the LSU and SSU rRNA sequence databases
               annotations for individual ORFs are created.
               The functional annotation is done by synthesizing  the annotations from the 
               reference databases, after removing redundant texts and also making sure 
               the annotations from different reference databases are in agreement with each other.
               The taxonomic annotation is done by applying the LCA (lowest common ancestor rule) 
               on the hits from the RefSeq NR protein database."""

     epilog = re.sub(r'\s+',' ', epilog)
     parser = OptionParser(usage = usage, epilog = epilog)

     parser.add_option("-b", "--blastoutput", dest="input_blastout", action='append', default=[],
                       help='blastout files in TSV format [at least 1 REQUIRED]')
     parser.add_option("-a", "--algorithm", dest="algorithm", default="BLAST", help="algorithm BLAST or LAST" )

     parser.add_option("-m", "--contig-map", dest="contig_map_file", default=None, help="contig map file" )
     
     parser.add_option("-d", "--dbasename", dest="database_name", action='append', default=[],
                       help='the database names [at least 1 REQUIRED]')

     parser.add_option("-D", "--blastdir", dest="blastdir",  default=None,
                       help='the blast dir where all the BLAST outputs are located')
     
     parser.add_option("-s", "--samplename", dest="sample_name",  default=None,
                       help='the sample name')
     
     parser.add_option("-w", "--weight_for_database", dest="weight_db", action='append', default=[], type='float',
                       help='the map file for the database  [at least 1 REQUIRED]')
     
     parser.add_option( "--rRNA_16S", dest="rRNA_16S", action="append", default=[], 
                       help='the 16s rRNA stats file [OPTIONAL]')
     
     parser.add_option( "--tRNA", dest="tRNA", action="append", default=[], 
                       help='the tRNA stats file [OPTIONAL]')
     
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
     
     cutoffs_group.add_option("--min_bsr", dest="min_bsr", type='float', default=0.00,
                       help='minimum BIT SCORE RATIO [default = 0.00 ] ')
     parser.add_option_group(cutoffs_group)
     
     
     output_options_group =  OptionGroup(parser, 'Output Options')
     output_options_group.add_option("--tax", dest="taxonomy", action='store_true', default=False,
                       help='add the taxonomy info [useful for refseq] ')
     parser.add_option_group(output_options_group)
     
     parser.add_option('-o' , "--output_gff", dest="output_gff",
                      help='the output gff file [REQUIRED]')
     
     parser.add_option('--output-comparative-annotation', dest="output_comparative_annotation",
                      help='the comparative output table [REQUIRED]')
     
     parser.add_option('--input_gff', dest='input_gff',
                metavar='INPUT', help='Unannotated gff file [REQUIRED]')



def check_arguments(opts, args):

    return True 

    if len(opts.input_blastout) == 0:
         eprintf("There should be at least one blastoutput file\n")  
         return False

    if len(opts.database_name) == 0:
         eprintf("There should be at least one database name\n")  
         return False

    if len(opts.weight_db) == 0:
         eprint("There should be at least one weight\n")  
         return False

    if len(opts.input_blastout) != len(opts.database_name) or\
         len(opts.input_blastout) !=  len(opts.weight_db) :
         eprint("The num of database names, blastoutputs and database map file should be equal\n")
         return False

    if opts.output_gff == None:
       eprintf("Must specify the output gff file\n")
       return False

    if opts.output_comparative_annotation == None:
       eprintf("Must specify the output tables for comparative annotation\n")
       return False

    if opts.input_gff == None:
       eprintf("Must specify the input gff file\n")
       return False

    return True

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
     attributes['seqname'] =  fields[0]   # this is a bit of a  duplication  
     attributes['source'] =  fields[1]
     attributes['feature'] =  fields[2]
     attributes['start'] =  int(fields[3])
     attributes['end'] =  int(fields[4])

     try:
        attributes['score'] =  float(fields[5])
     except:
        attributes['score'] =  fields[5]

     attributes['strand'] =  fields[6]
     attributes['frame'] =  fields[7]
     
     split_attributes(fields[8], attributes)

     if not fields[0] in contig_dict :
       contig_dict[fields[0]] = []

     contig_dict[fields[0]].append(attributes)
  

class GffFileParser(object):

    def __init__(self, gff_filename):
        self.Size = 10000
        self.i=0
        self.orf_dictionary = {}
        self.gff_beg_pattern = re.compile("^#")
        self.lines= []
        self.size=0
        try:
           self.gff_file = open( gff_filename,'r')
        except AttributeError:
           eprintf("Cannot read the map file for database : %s\n", dbname)
           exit_process()
  
    def __iter__(self):
        return self
 
    def refillBuffer(self):
       self.orf_dictionary = {}
       i = 0
       while i < self.Size:
          line=self.gff_file.readline()
          if not line:
            break
          if self.gff_beg_pattern.search(line):
            continue
          insert_orf_into_dict(line, self.orf_dictionary)
          i += 1

       self.orfs = self.orf_dictionary.keys()
       self.size = len(self.orfs)
       self.i = 0

    def next(self):
        if self.i == self.size:
           self.refillBuffer()

        if self.size==0:
           self.gff_file.close()
           raise StopIteration()

        #print self.i
        if self.i < self.size:
           self.i = self.i + 1
           return self.orfs[self.i-1]
       


def process_gff_file(gff_file_name, orf_dictionary):
     try:
        gfffile = open(gff_file_name, 'r')
     except IOError:
        eprintf("Cannot read file %s!\n", gff_file_name)

     gff_lines = gfffile.readlines()
     gff_beg_pattern = re.compile("^#")
     gfffile.close()
     
     count = 0
     for line in gff_lines:
        line = line.strip() 
        if gff_beg_pattern.search(line):
          continue
        insert_orf_into_dict(line, orf_dictionary)
        count += 1
        #if count %10000 == 0:
        #   print count


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
           

def write_annotation_for_orf(outputgff_file, candidatedbname, dbname_weight, results_dictionary, orf_dictionary, contig, candidate_orf_pos,  orfid):
   try:
      fields = [  'source', 'feature', 'start', 'end', 'score', 'strand', 'frame' ]


      output_line= orf_dictionary[contig][candidate_orf_pos]['seqname']

      for field in fields:
        # printf("\t%s", orf_dictionary[contig][candidate_orf_pos][field])
         output_line += "\t"+ str(orf_dictionary[contig][candidate_orf_pos][field])

      attributes = "ID="+orf_dictionary[contig][candidate_orf_pos]['id']
      attributes += ";" + "locus_tag="+orf_dictionary[contig][candidate_orf_pos]['locus_tag']
      attributes += ";" + "contig_length="+orf_dictionary[contig][candidate_orf_pos]['contig_length']
      attributes += ";" + "orf_length="+orf_dictionary[contig][candidate_orf_pos]['orf_length']
      attributes += ";" + "partial="+orf_dictionary[contig][candidate_orf_pos]['partial']
      attributes += ";" + "sourcedb="+candidatedbname
     
      if candidatedbname in results_dictionary:
         attributes += ";" + "annotvalue="+str(results_dictionary[candidatedbname][orfid]['value'])
         attributes += ";" + "ec="+str(results_dictionary[candidatedbname][orfid]['ec'])
         attributes += ";" + "product="+results_dictionary[candidatedbname][orfid]['product']
      else:
         attributes += ";" + "annotvalue="+str('0')
         attributes += ";" + "ec="+str('')
         attributes += ";" + "product="+'hypothetical protein'

      output_line += '\t' + attributes
      fprintf(outputgff_file, "%s\n", output_line);
   except:
      eprintf("ERROR : Failure to annotate in contig %s\n", contig)
      #print orf_dictionary[contig]
      print traceback.print_exc(10)
      exit_process()


def  write_16S_tRNA_gene_info(rRNA_dictionary, outputgff_file, tag):
      fields = [  'source', 'feature', 'start', 'end', 'score', 'strand', 'frame' ]
      for rRNA in rRNA_dictionary:
          output_line= rRNA_dictionary[rRNA]['seqname']
          for field in fields:
             output_line += "\t"+ str(rRNA_dictionary[rRNA][field])

          attributes = "ID="+rRNA_dictionary[rRNA]['seqname'] + tag
          attributes += ";" + "locus_tag="+rRNA_dictionary[rRNA]['seqname'] + tag
          attributes += ";" + "orf_length=" + str(rRNA_dictionary[rRNA]['orf_length'])
          attributes += ";" + "contig_length=" + str(rRNA_dictionary[rRNA]['contig_length'])
          attributes += ";" + "ec="
          attributes += ";" + "product="+rRNA_dictionary[rRNA]['product']
          output_line += '\t' + attributes
          fprintf(outputgff_file, "%s\n", output_line);


def process_rRNA_16S_stats(rRNA_16S_file, rRNA_16S_dictionary):
     try:
        taxonomy_file = open(rRNA_16S_file, 'r')
     except IOError:
        eprintf("Cannot read file %s!\n", rRNA_16S_file)
        exit_process()

     tax_lines = taxonomy_file.readlines()
     similarity_pattern = re.compile("similarity")
     evalue_pattern = re.compile("evalue")
     bitscore_pattern = re.compile("bitscore")
     taxonomy_pattern = re.compile("taxonomy")
     headerScanned = False
     for line in tax_lines:
         if headerScanned == False:
            if similarity_pattern.search(line) and evalue_pattern.search(line) and bitscore_pattern.search(line) and  taxonomy_pattern.search(line):
                headerScanned = True
            continue
         fields = [ x.strip() for x in line.split('\t') ]
         if len(fields) >=6:
           if fields[1]!='-':
              rRNA_16S_dictionary[fields[0]] =  [ fields[1], fields[2], fields[5] ]
           else:
              if len(fields) >=12:
                 if fields[7]!='-':
                     rRNA_16S_dictionary[fields[0]] =  [ fields[7], fields[8], fields[11] ]

     taxonomy_file.close()

def process_tRNA_stats(tRNA_stats_file, tRNA_dictionary):
     try:
        tRNA_file = open(tRNA_stats_file, 'r')
     except IOError:
        eprintf("Cannot read file %s!\n", tRNA_stats_file)
        exit_process()
     tRNA_lines = tRNA_file.readlines()

     sequence_name_pattern = re.compile("sequence name", re.I)
     number_pattern = re.compile("number", re.I)

     headerScanned = False
     for line in tRNA_lines:
         if number_pattern.search(line):
            continue
         if headerScanned == False:
            if sequence_name_pattern.search(line):
                headerScanned = True
            continue
         fields = [ x.strip() for x in line.split('\t') ]
         if len(fields) >=6:
              tRNA_dictionary[fields[0]] =  [ fields[3], fields[4], fields[5], fields[1] ]

# this adds the features and attributes to  be added to the gff file format for the tRNA genes
def add_tRNA_genes(tRNA_dictionary, tRNA_gff_dictionary, contig_lengths) :

     for tRNA in tRNA_dictionary: 
        try:
           orf_length = abs(int( tRNA_dictionary[tRNA][1] )-int( tRNA_dictionary[tRNA][0] )) + 1
        except:
           orf_length = 0

        if tRNA in  contig_lengths: 
           contig_length = contig_lengths[tRNA]
        else:
           contig_length = 0

        dict = { 'id':tRNA, 'seqname': tRNA, 'start':str(tRNA_dictionary[tRNA][0]), 'end':str(tRNA_dictionary[tRNA][1]),\
                 'strand':tRNA_dictionary[tRNA][2], 'score':" ", 'orf_length':str(orf_length),\
                 'contig_length':str(contig_length),\
                 'feature':'tRNA', 'source':'tranScan-1.4', 'frame':0, 'product':'tRNA-' + tRNA_dictionary[tRNA][3], 'ec':'' }      
        tRNA_gff_dictionary[tRNA] = dict.copy() 


# this adds the features and attributes to  be added to the gff file format for the 16S rRNA genes
def add_16S_genes(rRNA_16S_dictionary, rRNA_dictionary, contig_lengths) :

    for rRNA in rRNA_16S_dictionary: 
        try:
           orf_length = abs(int( tRNA_dictionary[rRNA][1] )-int( tRNA_dictionary[rRNA][0] )) + 1
        except:
           orf_length = 0

        if rRNA in  contig_lengths: 
           contig_length = contig_lengths[rRNA]
        else:
           contig_length = 0


        dict = { 'id':rRNA, 'seqname': rRNA, 'start':str(rRNA_16S_dictionary[rRNA][0]), 'end':str(rRNA_16S_dictionary[rRNA][1]),\
                 'strand':'+', 'score':str(rRNA_16S_dictionary[rRNA][2]),  'orf_length':str(orf_length),\
                 'contig_length':str(contig_length),\
                 'feature':'CDS', 'source':'BLAST Search', 'frame':0, 'product':'16S rRNA', 'ec':'' }      
        rRNA_dictionary[rRNA] = dict.copy() 

    
def create_annotation(dbname_weight, results_dictionary, input_gff,  rRNA_16S_stats_files, tRNA_stats_files,  output_gff, output_comparative_annotation, contig_lengths):
    orf_dictionary={}
#    process_gff_file(input_gff, orf_dictionary)
    gffreader = GffFileParser(input_gff)

    output_gff_tmp = output_gff + ".tmp"
    outputgff_file = open( output_gff_tmp, 'w')
    output_comp_annot_file1 = open( output_comparative_annotation + '.1.txt', 'w')
    output_comp_annot_file2 = open( output_comparative_annotation + '.2.txt', 'w')

    output_comp_annot_file1_Str = 'orf_id\tref dbname\tEC\tproduct\tvalue'
    fprintf(output_comp_annot_file1,'%s\n', output_comp_annot_file1_Str)

    output_comp_annot_file2_Str = 'orf_id'
    dbnames = dbname_weight.keys()
    for dbname in dbnames:
         weight = dbname_weight[dbname]
         output_comp_annot_file2_Str += '\t{0}(EC) \t{0}(product)\t{0}(value)'.format(dbname)
    fprintf(output_comp_annot_file2,'%s\n', output_comp_annot_file2_Str)
       

#    gffreader = GffReader(input_gff)
    for contig in  gffreader:
       count = 0
       for orf in  gffreader.orf_dictionary[contig]:
         #print orf['id']
         value = 0.0001
         success =False
         output_comp_annot_file1_Str = ''
         output_comp_annot_file2_Str = ''
         for dbname in dbnames:
            weight = dbname_weight[dbname]
            value = 0
            if orf['id'] in results_dictionary[dbname]:
                if value < results_dictionary[dbname][orf['id']]['value']:
                    value = results_dictionary[dbname][orf['id']]['value']
                    candidatedbname=dbname
                    success =True
                    candidate_orf_pos = count 

                    if output_comp_annot_file1_Str:
                        output_comp_annot_file1_Str += '{0}\t{1}\t{2}\t{3}\t{4}\n'.format('', dbname,\
                               results_dictionary[dbname][orf['id']]['ec'],\
                               results_dictionary[dbname][orf['id']]['product'],\
                               str(results_dictionary[dbname][orf['id']]['value']*float(weight)))
                    else:
                        output_comp_annot_file1_Str += '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(orf['id'], dbname,\
                               results_dictionary[dbname][orf['id']]['ec'],\
                               results_dictionary[dbname][orf['id']]['product'],\
                               str(results_dictionary[dbname][orf['id']]['value']*float(weight)))


                    if output_comp_annot_file2_Str:
                        output_comp_annot_file2_Str += '\t{0}\t{1}\t{2}'.format(\
                               results_dictionary[dbname][orf['id']]['ec'],\
                               results_dictionary[dbname][orf['id']]['product'],\
                               str(results_dictionary[dbname][orf['id']]['value']*float(weight)))
                    else:
                        output_comp_annot_file2_Str += '{0}\t{1}\t{2}\t{3}'.format(orf['id'], 
                               results_dictionary[dbname][orf['id']]['ec'],\
                               results_dictionary[dbname][orf['id']]['product'],\
                               str(results_dictionary[dbname][orf['id']]['value']*float(weight)))

            else: 
                if not output_comp_annot_file1_Str:
                   output_comp_annot_file1_Str += '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(orf['id'], '','','','')

                if output_comp_annot_file2_Str:
                   output_comp_annot_file2_Str += '\t{0}\t{1}\t{2}'.format('', '','')
                else:
                   output_comp_annot_file2_Str += '{0}\t{1}\t{2}\t{3}'.format(orf['id'], '','','','')

         if success:  # there was a database hit
            fprintf(output_comp_annot_file1,'%s\n', output_comp_annot_file1_Str)
            fprintf(output_comp_annot_file2,'%s\n', output_comp_annot_file2_Str)
            write_annotation_for_orf(outputgff_file, candidatedbname, dbname_weight, results_dictionary, gffreader.orf_dictionary, contig, candidate_orf_pos,  orf['id']) 
         else:   # if it was not  a hit then it is a hypothetical protein
            #print gffreader.orf_dictionary
            write_annotation_for_orf(outputgff_file, 'None', '0', results_dictionary, gffreader.orf_dictionary, contig, count, orf['id']) 
         
         count +=1  #move to the next orf

       #del orf_dictionary[contig]   
    output_comp_annot_file1.close()
    output_comp_annot_file2.close()

    # now deal with the rRNA sequences  if there is rRNA stats file
    if len(rRNA_16S_stats_files) > 0 :
       rRNA_16S_dictionary={} 
       for rRNA_16S_stats_file in rRNA_16S_stats_files:
          process_rRNA_16S_stats(rRNA_16S_stats_file, rRNA_16S_dictionary)

       rRNA_dictionary = {}
       add_16S_genes(rRNA_16S_dictionary, rRNA_dictionary, contig_lengths) 
       write_16S_tRNA_gene_info(rRNA_dictionary, outputgff_file, '_rRNA')

    # now deal with the tRNA sequences  if there is tRNA stats file
    if len(tRNA_stats_files) > 0 :
       tRNA_dictionary={} 
       for tRNA_stats_file in tRNA_stats_files:
          process_tRNA_stats(tRNA_stats_file, tRNA_dictionary)

       tRNA_gff_dictionary = {}
       add_tRNA_genes(tRNA_dictionary, tRNA_gff_dictionary, contig_lengths) 
       write_16S_tRNA_gene_info(tRNA_gff_dictionary, outputgff_file, '_tRNA')
       #print tRNA_dictionary


    outputgff_file.close()     
    rename(output_gff_tmp, output_gff)


def process_product(product, database, similarity_threshold=0.9):
    """Returns the best set of products from the list of (*database*,
    *product*) tuples *products*.

    Each product in the set is first trimmed down, removing database-specific
    information.

    The set is then determined by first sorting the products by length
    (ascending), and then, for each product, sequentially applying the longest
    common substring algorithm to determine the similarity between the product
    and already determined products. If this similarity is greater than the
    specified *similarity_threshold*, the longer of the two products is chosen
    to be a determined product.
    """

    processed_product = ''

    # print 'dbase', database
    # COG
    if database == 'cog':
        results = re.search(r'Function: (.+?) #', product)
        if results:
           processed_product=results.group(1)

    # KEGG: split and process

    elif database == 'kegg':
        kegg_products = re.split(r'\s*;\s+', product)
        for kegg_product in kegg_products:
            # Toss out organism:ID pairs, gene names, and KO IDs
            kegg_product = re.sub(r'^lcl[|]', '', kegg_product)
            kegg_product = re.sub(r'[a-z]{3}:\S+', '', kegg_product)
            kegg_product = kegg_product.strip()
            kegg_product = re.sub(r'(, \b[a-z]{3}[A-Z]?\b)+', '', kegg_product)
            kegg_product = re.sub(r'^\b[a-z]{3}[A-Z]?\b', '', kegg_product)
            # get KO number 
            kegg_product = re.sub(r'\bK\d{5}\b', '', kegg_product)

            # Also toss out anything between square brackets
            kegg_product = re.sub(r'\[.*\]', '', kegg_product)

            if kegg_product.strip():
                processed_product=kegg_product.strip()
                

    # RefSeq: split and process

    elif database == 'refseq':
        for subproduct in product.split('; '):
            subproduct = re.sub(r'[a-z]{2,}\|(.+?)\|\S*', '', subproduct)
            subproduct = re.sub(r'\[.+?\]', '', subproduct)
            if subproduct.strip():
                processed_product=subproduct.strip()

    # MetaCyc: split and process

    elif database == 'metacyc':
        # Pull out first name after the accession code:
        product_name = product.split('#')[0].strip()
        product_name = re.sub(r'^[^ ]* ', '', product_name)
        product_name = re.sub(r' OS=.*', '', product_name)
        if product_name:
            processed_product=product_name

    # Seed: split and process

    elif database == 'seed':
        for subproduct in product.split('; '):
            #subproduct = re.sub(r'[a-z]{2,}\|(.+?)\|\S*', '', subproduct)
            subproduct = re.sub(r'\[.+?\]', '', subproduct)
            subproduct = re.sub(r'\(.+?\)', '', subproduct)
            if subproduct.strip():
                processed_product=subproduct.strip()

    elif database == 'cazy':
        for subproduct in product.split('; '):
            #subproduct = re.sub(r'[a-z]{2,}\|(.+?)\|\S*', '', subproduct)
            subproduct = re.sub(r'\[.+?\]', '', subproduct)
            subproduct = re.sub(r'\(.+?\)', '', subproduct)
            if subproduct.strip():
                processed_product=subproduct.strip()
                print processed_product

    # MetaCyc: split and process

    # Generic
    else:
        processed_product=strip_taxonomy(product)

    words = [ x.strip() for x in processed_product.split() ]
    filtered_words =[]
    underscore_pattern = re.compile("_")
    arrow_pattern = re.compile(">")
    for word in words:
       if not  underscore_pattern.search(word) and not arrow_pattern.search(word):
           filtered_words.append(word)
    

    #processed_product = ' '.join(filtered_words)
    # Chop out hypotheticals
    processed_product = remove_repeats(filtered_words)
    processed_product = re.sub(';','',processed_product)

    processed_product = re.sub(r'hypothetical protein','', processed_product)

    return processed_product

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
        self.dbname = dbname
        self.blastoutput = blastoutput
        self.i=1
        self.data = {}
        self.fieldmap={}
        self.seq_beg_pattern = re.compile("#")

        try:
           self.blastoutputfile = open( blastoutput,'r')
           self.lines=self.blastoutputfile.readlines()
           self.blastoutputfile.close()
           self.size = len(self.lines)
           if not self.seq_beg_pattern.search(self.lines[0]) :
              exit_process("First line must have field header names and begin with \"#\"")
           header = self.lines[0].replace('#','',1)
           fields = [ x.strip()  for x in header.rstrip().split('\t')]
           k = 0 
           for x in fields:
             self.fieldmap[x] = k 
             k+=1
           eprintf("\nProcessing database : %s\n", dbname)
           
        except AttributeError:
           eprintf("Cannot read the map file for database :%s\n", dbname)
           exit_process()
  
    def __iter__(self):
        return self
    count = 0 
    def next(self):
        if self.i < self.size:
           
           try:
              fields = [ x.strip()  for x in self.lines[self.i].split('\t')]
              #print self.fieldmap['ec'], fields, self.i,  self.blastoutput
              self.data['query'] = fields[self.fieldmap['query']]
              self.data['q_length'] = int(fields[self.fieldmap['q_length']])
              self.data['bitscore'] = float(fields[self.fieldmap['bitscore']])
              self.data['bsr'] = float(fields[self.fieldmap['bsr']])
              self.data['expect'] = float(fields[self.fieldmap['expect']])
              self.data['identity'] = float(fields[self.fieldmap['identity']])
              self.data['ec'] = fields[self.fieldmap['ec']]
              self.data['product'] = re.sub(r'=',' ',fields[self.fieldmap['product']])
              self.i = self.i + 1
              return self.data
           except:
              print self.lines[self.i]
              print data
              sys.exit(0)
              return None
        else:
           raise StopIteration()
              
def isWithinCutoffs(data, cutoffs):
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

    return True


def word_information(string_of_words):
    words = [ x.strip() for x in string_of_words.split() ]

    information = 0
    wordlist = {}
    underscore_pattern = re.compile("_")
    for word in words:
       if not word in ['', 'is', 'have', 'has', 'will', 'can', 'should',  'in', 'at', 'upon', 'the', 'a', 'an', 'on', 'for', 'of', 'by', 'with' ,'and',  '>', 'predicted', 'protein', 'conserved' ]:
          if not underscore_pattern.search(word):
             wordlist[word]=1

    #print string_of_words
    #print wordlist
    #print len(wordlist)
    return len(wordlist)


def compute_annotation_value(data):
    score = 0;
    if len(data['ec'] ) > 0:
       score += 10

    if not re.search(r'hypothetical protein', data['product']):
       score += word_information(data['product'])

    return score
        

# compute the refscores
def process_parsed_blastoutput(dbname, weight,  blastoutput, cutoffs, annotation_results):
    blastparser =  BlastOutputTsvParser(dbname, blastoutput)

    fields = ['q_length', 'bitscore', 'bsr', 'expect', 'aln_length', 'identity', 'ec' ]
    if cutoffs.taxonomy:
       fields.append('taxonomy')
    fields.append('product')

    annotation = {}
    for data in blastparser:
        #if count%10000==0:
        #   print count

        if isWithinCutoffs(data, cutoffs) :
           #print data['query'] + '\t' + str(data['q_length']) +'\t' + str(data['bitscore']) +'\t' + str(data['expect']) +'\t' + str(data['identity']) + '\t' + str(data['bsr']) + '\t' + data['ec'] + '\t' + data['product']
#           if data['query'] =='NapDC_illum_asm_188606_0':

    #       print dbname 
           annotation['bsr'] = data['bsr']
           annotation['ec'] = data['ec']
           annotation['product'] = strip_taxonomy(process_product(data['product'], dbname) )
           annotation['value'] = compute_annotation_value(annotation)*weight
         #  print annotation
           

              #sys.exit(0)

           if not data['query'] in annotation_results:
               annotation_results[data['query']] = {'value':0}

           if annotation_results[data['query']]['value'] <= annotation['value'] :
                 annotation_results[data['query']] = annotation.copy()

#    add_refscore_to_file(blastoutput,refscore_file, allNames)
    count =  len(annotation_results.keys())
    return count

def read_contig_lengths(contig_map_file, contig_lengths):
     try:
        mapfile = open(contig_map_file, 'r')
     except IOError:
        print "Cannot read file " + contig_map_file + " !"

     mapfile_lines = mapfile.readlines()
     mapfile.close()
     
     for line in mapfile_lines:
        line = line.strip() 
        fields = [ x.strip() for x in line.split('\t') ]
        if len(fields) != 3: 
            contig_lengths = {}
            return 
        contig_lengths[fields[0] ] = int(fields[2])
     
def getBlastFileNames(opts) :
    
    database_names = []
    parsed_blastouts = [] 
    weight_dbs = []

    dbnamePATT = re.compile(r'' + opts.blastdir + '*' + opts.sample_name + '*[.](.*)[.]' + opts.algorithm.upper() + 'out.parsed.txt')

    blastOutNames = glob(opts.blastdir + '*' + opts.algorithm.upper() + 'out.parsed.txt')  
    for blastoutname in blastOutNames :
        result = dbnamePATT.search(blastoutname)
        if result:
            dbname = result.group(1)
            database_names.append(dbname)
            parsed_blastouts.append(blastoutname)
            weight_dbs.append(1)

    return database_names, parsed_blastouts, weight_dbs


# the main function
def main(argv, errorlogger =None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    if not check_arguments(opts, args):
       print usage
       sys.exit(0)

    results_dictionary={}
    dbname_weight={}

    contig_lengths = {}     
    read_contig_lengths(opts.contig_map_file, contig_lengths) 

    
    if opts.blastdir !=None and opts.sample_name != None:
        try:
           database_names, input_blastouts, weight_dbs = getBlastFileNames(opts) 
        except:
           print traceback.print_exc(10)
           pass
    else:
        database_names = opts.database_name
        input_blastouts = opts.input_blastout 
        weight_dbs = opts.weight_db 


    priority = 6000
    count_annotations = {}
    for dbname, blastoutput, weight in zip( database_names, input_blastouts, weight_dbs): 
        results_dictionary[dbname]={}
        dbname_weight[dbname] = weight
        count = process_parsed_blastoutput( dbname, weight, blastoutput, opts, results_dictionary[dbname])
        if runstatslogger!=None:
           runstatslogger.write("%s\tProtein Annotations from %s\t%s\n" %( str(priority), dbname, str(count)))
        count_annotations 
        priority += 1

    for dbname in results_dictionary: 
      for seqname in results_dictionary[dbname]: 
         count_annotations[seqname] = True      
    count = len(count_annotations)
    if runstatslogger!=None:
       runstatslogger.write("%s\tTotal Protein Annotations\t%s\n" %( str(priority),  str(count)))
        

    #create the annotations from he results
    
    #print results_dictionary
    create_annotation(dbname_weight, results_dictionary, opts.input_gff, opts.rRNA_16S, opts.tRNA, opts.output_gff, opts.output_comparative_annotation, contig_lengths)


def MetaPathways_annotate_fast(argv, errorlogger = None, runstatslogger = None):       
    createParser()
    errorlogger.write("#STEP\tANNOTATE_ORFS\n")
    main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger)
    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

