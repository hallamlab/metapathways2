#!/usr/bin/python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import optparse
    import csv
    from os import makedirs, path, listdir, remove, rename
    import shutil
    import traceback
    import sys
    import logging.handlers
    import re
    from glob import glob
    from libs.python_modules.utils.utils import *
    from libs.python_modules.utils.sysutil import pathDelim, genbankDate, getstatusoutput
    from libs.python_modules.parsers.parse  import parse_parameter_file
except:
    print """ Could not load some user defined  module functions"""
    print """ Make sure your typed 'source MetaPathwaysrc'"""
    print """ """
    sys.exit(3)

PATHDELIM = pathDelim()

def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

def files_exist( files ):
     for file in files:
        if not path.exists(file):
           print 'Could not read File ' +  file
           print 'Please make sure these sequences are in the \"blastDB\" folder'
           sys.exit(3)
           return False
     return True


def removeDir(origFolderName):
    folderName = origFolderName + PATHDELIM + '*'
    files = glob(folderName)
    for  f in files:
       remove(f)
    if path.exists(origFolderName):
       shutil.rmtree(origFolderName)
    

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
    
     if not fields[0] in contig_dict :
       contig_dict[fields[0]] = []

 #    print attributes
     contig_dict[fields[0]].append(attributes)


def get_sequence_name(line): 
     fields = re.split(' ', line)
     name = re.sub('>','',fields[0])
     #print name
     return name
     


note = """GFF File Format
Fields

Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
"""

#def get_date():


def process_gff_file(gff_file_name, output_filenames, nucleotide_seq_dict, protein_seq_dict, input_filenames):
     #print output_filenames
     try:
        gfffile = open(gff_file_name, 'r')
     except IOError:
        print "Cannot read file " + gff_file_name + " !"

     sample_name= re.sub('.annotated.gff', '', gff_file_name)
     sample_name= re.sub('.annot.gff', '', gff_file_name) # Niels: somewhere we changed the file name?
     sample_name= re.sub(r'.*[/\\]', '', sample_name)

     gff_lines = gfffile.readlines()
     gff_beg_pattern = re.compile("^#")
     gfffile.close()
     
     contig_dict={} 
     count = 0
     for line in gff_lines:
        line = line.strip() 
        if gff_beg_pattern.search(line):
          continue
        insert_orf_into_dict(line, contig_dict)
        #if count > 100:
        #   break 
        #count = count + 1
#        print line


     if "gbk" in output_filenames:
       write_gbk_file(output_filenames['gbk'], contig_dict, sample_name, nucleotide_seq_dict, protein_seq_dict)

     if "sequin" in output_filenames:
       write_sequin_file(output_filenames['sequin'], contig_dict, sample_name, nucleotide_seq_dict, protein_seq_dict, input_filenames)

     if "ptinput" in output_filenames:
       write_ptinput_files(output_filenames['ptinput'], contig_dict, sample_name, nucleotide_seq_dict, protein_seq_dict)

# this function creates the pathway tools input files
def  write_ptinput_files(output_dir_name, contig_dict, sample_name, nucleotide_seq_dict, protein_seq_dict):

     try:
        #print output_dir_name
        removeDir(output_dir_name)
        #print output_dir_name
        makedirs(output_dir_name)
        genetic_elementsfile = open(output_dir_name + "/.tmp.genetic-elements.dat", 'w')
        zerofastafile = open(output_dir_name + "/.tmp.0.fasta", 'w')
        zeropffile = open(output_dir_name + "/tmp.0.pf", 'w')
        organism_paramsfile = open(output_dir_name + "/.tmp.organism-params.dat", 'w')
     except:
        print "cannot create the pathway tools files"
        print "perhaps there is already a folder " + output_dir_name
        traceback.print_exc(file=sys.stdout)

     count =0 
     outputStr=""
     # iterte over every sequence
     startbase = 0
     endbase = 0
     fastaStr =""
     cumulFastaLength=0
     fprintf(zerofastafile, ">0\n")
     for key in contig_dict:
        first = True
        if count %10000 == 0:
           #print "count " + str(count)
           #outputfile.write(outputStr)
           outputStr=""
        count+=1

        for attrib in contig_dict[key]:     
           id  = attrib['id']
           try:
              protein_seq = protein_seq_dict[id]
           except:
              protein_seq = ""
           try:   
              if attrib['product']=='hypothetical protein':
                 continue
           except:
              print attrib
              sys.exit(0)

           startbase = cumulFastaLength + attrib['start']
           endbase =  cumulFastaLength + attrib['end']

           try: 
              fprintf(zeropffile, "ID\t%s\n", attrib['id'])
           except:
              fprintf(zeropffile, "ID\t%s \n", attrib['id'])

           try: 
              fprintf(zeropffile, "NAME\t%s\n", attrib['id'])
           except:
              fprintf(zeropffile, "NAME\t%s \n", attrib['id'])

           try: 
              fprintf(zeropffile, "STARTBASE\t%s\n", startbase)
           except:
              fprintf(zeropffile, "STARTBASE\t%s\n", startbase)

           try: 
              fprintf(zeropffile, "ENDBASE\t%s\n", endbase)
           except:
              fprintf(zeropffile, "ENDBASE\t%s\n", endbase)

           try: 
              fprintf(zeropffile, "PRODUCT\t%s\n", attrib['product'])
           except:
              fprintf(zeropffile, "PRODUCT\t%s \n", 'hypothetical protein')

           fprintf(zeropffile, "PRODUCT-TYPE\tP\n")
           try:
             if  len(attrib['ec']) > 0:
                fprintf(zeropffile, "EC\t%s\n", attrib['ec'])
           except: 
                pass
 #             print attrib
 #             sys.exit(0)

           fprintf(zeropffile, "//\n")

        try:
           dna_seq =  nucleotide_seq_dict[key]
           dna_seq += 'NNNNNNNNNN'
        except:
           print "Key missing for " + key
           continue

        fastaStr+=(wrap("",0,62, dna_seq)+'\n')
        cumulFastaLength += len(dna_seq)

        fprintf(zerofastafile, "%s",  fastaStr)
        fastaStr=""

     zerofastafile.close()
     zeropffile.close()
     rename(output_dir_name + "/.tmp.0.fasta", output_dir_name + "/0.fasta")
     rename(output_dir_name + "/tmp.0.pf", output_dir_name + "/0.pf")

     # Niels: removing annotated.gff from sample_name
     sample_name = re.sub(".annot.gff", '', sample_name)
     sample_name = re.sub('.*/', '', sample_name)
     sample_name = re.sub(r'[\\].', '', sample_name)
     
     # Niels: trim sample_name to less than 35 characters 
     # as it causes PGDB creation to fail
     if (len(sample_name) > 35):
        sample_name = sample_name[0:35]
     
     
     if not sample_name[0].isalpha() :
        sample_name = 'E' + sample_name

     fprintf(organism_paramsfile,"ID\t%s\n",sample_name)
     fprintf(organism_paramsfile,"STORAGE\tFILE\n")
     fprintf(organism_paramsfile,"NAME\t%s\n",sample_name)
     fprintf(organism_paramsfile,"ABBREV-NAME\t%s\n",sample_name)
     fprintf(organism_paramsfile,"STRAIN\t1\n")
     fprintf(organism_paramsfile,"RANK\t|species|\n")
     fprintf(organism_paramsfile,"NCBI-TAXON-ID\t12908\n")
     organism_paramsfile.close()
     rename(output_dir_name + "/.tmp.organism-params.dat", output_dir_name + "/organism-params.dat")

     fprintf(genetic_elementsfile,"ID\t0\n")
     fprintf(genetic_elementsfile,"NAME\t0\n")
     fprintf(genetic_elementsfile,"TYPE\t:READ/CONTIG\n")
     fprintf(genetic_elementsfile,"ANNOT-FILE\t0.pf\n")
     fprintf(genetic_elementsfile,"SEQ-FILE\t0.fasta\n")
     fprintf(genetic_elementsfile,"//")
     genetic_elementsfile.close()
     rename(output_dir_name + "/.tmp.genetic-elements.dat", output_dir_name + "/genetic-elements.dat")

def get_parameter(config_params, category, field, default = None):
     if config_params == None:
       return default
 
     if category in config_params:
        if field in config_params[category]:
             if config_params[category][field]: 
                return config_params[category][field]
             else:
                return default
        else:
             return default
     return default

#this function creates the sequin  file from the gff, protein and nucleotide sequences  
def  write_sequin_file(tbl_file_name, contig_dict, sample_name, nucleotide_seq_dict, protein_seq_dict, sequin_input_files):
      

     sequin_src_filename = re.sub(r'tbl$', 'src', tbl_file_name)
     sequin_output_fasta = re.sub(r'tbl$', 'fasta', tbl_file_name)
     sequin_output_sbt = re.sub(r'tbl$', 'sbt', tbl_file_name)

     shutil.copy(sequin_input_files['sequin_fasta'], sequin_output_fasta)
     shutil.copy(sequin_input_files['sequin_sbt_file'], sequin_output_sbt)
     sequin_required_files = { 'fasta': sequin_output_fasta, 'tbl': tbl_file_name, 'src': sequin_src_filename, 'tbl2asn': sequin_input_files['sequin_tbl2asn'], 'sbt': sequin_output_sbt }

     outputfile = open(tbl_file_name, 'w')
     #print contig_dict
    
     count =0 
     outputStr=""
     for key in contig_dict:
        first = True
        if count %10000 == 0:
           #print "count " + str(count)
           outputfile.write(outputStr)
           outputStr=""
        count+=1

        for attrib in contig_dict[key]:     
           id  = attrib['id']
           try:
              protein_seq = protein_seq_dict[id]
           except:
              protein_seq = ""
              None
           
           definition = sample_name
           accession = '.'
           version = '.' +spaces(10) + "GI:."
           dblink = sample_name
           keywords = '.'
           source = sample_name
           organism = sample_name
           if first:   
              first = False
              try:
                dna_seq =  nucleotide_seq_dict[key]
                dna_seq_formatted =  format_sequence_origin(dna_seq)
                dna_length = len(dna_seq)
                sourceStr = "1.." + str(dna_length)
              except:
                dna_seq = ""
                dna_seq_formatted =  ""
                dna_length = 0
                sourceStr ="0..0"

              outputStr+=(">Feature %s\n" % (key))
              outputStr+=re.sub('\.\.','\t',sourceStr)+'\t'+"REFERENCE" + '\n'
            
           startPrefix = ''
           endPrefix = ''
           if 'partial' in attrib:
               if attrib['partial']=='10':
                 startPrefix = '<'
               if attrib['partial']=='01':
                 endPrefix = '>'
               if attrib['partial']=='11':
                 startPrefix = '<'
                 endPrefix = '>'


           if 'start' in attrib and 'end' in attrib:
              if 'strand' in attrib:
                 if attrib['strand']=='-':
                     geneLoc = str(attrib['end']) + endPrefix +'\t' + startPrefix +  str(attrib['start'])
                 else:
                     geneLoc = startPrefix + str(attrib['start']) +'\t' + str(attrib['end']) + endPrefix
              outputStr+=geneLoc + '\t' + "gene" + '\n'
 

           if 'locus_tag' in attrib:
               locus_tag = "gene" + '\t' + attrib['locus_tag'] 
               outputStr+='\t\t\t' + locus_tag +'\n'


           outputStr+=geneLoc + '\t' + "CDS" + '\n'

           if 'product' in attrib:
              product_tag = "product" + '\t' + attrib['product'] 
              outputStr+='\t\t\t' + product_tag +'\n'

     outputfile.write(outputStr)
     outputfile.close() 

     outputsrcfile = open(sequin_src_filename, 'w')
     ncbi_sequin_params = parse_parameter_file(sequin_input_files['sequin_params'])
       
     headers =  ['Collection_date', 'Country', 'isolation_source',  'Lat_Lon', 'Organism', 'environmental_sample']
     
     header_values = {}
     headerStr = 'Sequence_ID'
     for header_name in headers:
        headerStr += '\t' + header_name
        header_values[header_name]= get_parameter(ncbi_sequin_params, 'SequinHeader', header_name, default='__'+ header_name + '__')

    
     valueStr =""
     for header_name in headers:
         valueStr += "\t" + header_values[header_name]

     fprintf(outputsrcfile, "%s\n", key + headerStr)
     for key in contig_dict:
        fprintf(outputsrcfile, "%s\n", key + valueStr)
     outputsrcfile.close()

     # Now open a pipe process and run the tbl2asn script on the sequin input
       
     for file in sequin_required_files:
        if not path.exists(sequin_required_files[file]):
           print "Could not find file : " + sequin_required_files[file]
           print "Make sure all of the following files are present :"
           for file in sequin_required_files:
                print file
           sys.exit(0)

     args = [ sequin_required_files['tbl2asn'], '-t', sequin_required_files['sbt'] , '-i', sequin_required_files['fasta'], '-a', 's', '-V', 'v']  
     command = ' '.join(args)
     result = getstatusoutput(command)
     if result[0] == 0 :
         print "Successfully created the SEQUIN file"
       
     
     #print contig_dict




#this function creates the genbank file from the gff, protein and nucleotide sequences  
def  write_gbk_file(output_file_name, contig_dict, sample_name, nucleotide_seq_dict, protein_seq_dict):

     date = genbankDate()
     output_file_name_tmp = output_file_name + ".tmp"
     outputfile = open(output_file_name_tmp, 'w')
     #print contig_dict
    
     count =0 
     outputStr=""
     for key in contig_dict:
        first = True
        if count %10000 == 0:
           #print "count " + str(count)
           outputfile.write(outputStr)
           outputStr=""
        count+=1

        for attrib in contig_dict[key]:     
           id  = attrib['id']
           try:
              protein_seq = protein_seq_dict[id]
           except:
              protein_seq = ""
              None
           
           definition = sample_name
           accession = '.'
           version = '.' +spaces(10) + "GI:."
           dblink = sample_name
           keywords = '.'
           source = sample_name
           organism = sample_name
           if first:   
              first = False
              try:
                dna_seq =  nucleotide_seq_dict[key]
                dna_seq_formatted =  format_sequence_origin(dna_seq)
                dna_length = len(dna_seq)
                sourceStr = "1.." + str(dna_length)
              except:
                dna_seq = ""
                dna_seq_formatted =  ""
                dna_length = 0
                sourceStr ="0..0"

              outputStr+=("LOCUS       %-18s  %4d bp   DNA           BCT      %-11s\n" % (key, dna_length,  date))
              outputStr+=(wrap("DEFINITION  ",12,74, definition)+'\n')
              outputStr+=(wrap("ACCESSION   ", 12, 74, accession)+'\n')
              outputStr+=(wrap("VERSION     ", 12, 74, version)+'\n')
              outputStr+=(wrap("DBLINK      ", 12, 74, dblink)+'\n')
              outputStr+=(wrap("KEYWORDS    ", 12, 74,keywords)+'\n')
              outputStr+=(wrap("SOURCE    ", 12, 74, keywords)+'\n')
              outputStr+=(wrap("  ORGANISM  ",12, 74, organism)+'\n')
              outputStr+=(wrap("", 12, 74, "Metagenome")+'\n')
              outputStr+=( wrap("REFERENCE   ",12,74, "1  (bases 1 to XXXXX)")+'\n')
              outputStr+=( wrap("  AUTHORS   ",12,74, "YYYYYY,X.")+'\n')
              outputStr+=( wrap("  CONSRTM   ",12,74, "XXXXX")+'\n')
              outputStr+=( wrap("  TITLE     ",12,74, "XXXXX")+'\n')
              outputStr+=( wrap("  JOURNAL   ",12,74,"XXXXX")+'\n')
              outputStr+=( wrap("   PUBMED   ",12,74,"XXXXX")+'\n')
              outputStr+=( wrap("  REMARK   ",12,74, "XXXXX")+'\n')
              outputStr+=( wrap("COMMENT     ", 12, 74,"PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review   COMPLETENESS: XXXXX")+'\n')
             
              outputStr+=( wrap("FEATURES ",21,74,"Location/Qualifiers") +'\n')
              outputStr+=( wrap("     source",21,74,sourceStr) +'\n')
              outputStr+=( wrap("",21,74,"/organism=\"" + sourceStr +"\"") +'\n')
              outputStr+=( wrap("",21,74,"/strain=\"1\"")+'\n')
              outputStr+=( wrap("",21,74,"/chromosome=\"1\"") +'\n')
            



           if 'start' in attrib and 'end' in attrib:
               geneLoc = str(attrib['start']) +".." + str(attrib['end'])
           else:
               geneLoc = "0..0"

           if 'strand' in attrib:
              if attrib['strand']=='-':
                 geneLoc='complement' + '(' + geneLoc +')'

           outputStr+=( wrap("     gene",21,74,geneLoc) +'\n')
           if 'locus_tag' in attrib:
               locus_tag = "/locus_tag=" + "\"" + attrib['locus_tag'] + "\""
           else:
               locus_tag = "/locus_tag" + "\"\"" 
           outputStr+=( wrap("",21,74,locus_tag) +'\n')
           outputStr+=( wrap("     CDS",21,74,geneLoc) +'\n')
           if 'product' in attrib:
              product="/product=" + attrib['product']
           else:
              product="/product=\"\""
           outputStr+=( wrap("",21,74,product) +'\n')
           outputStr+=( wrap("",21,74,locus_tag) +'\n')

           codon_start="/codon_start=1"
           translation_table="/transl_table=11"
           outputStr+=( wrap("",21,74,codon_start) +'\n')
           outputStr+=( wrap("",21,74,translation_table) +'\n')
       
           translation= "/translation="+ protein_seq
           outputStr+=( wrap("",21,74,translation) +'\n')

        outputStr+=(wrap("ORIGIN", 21, 74, "")+'\n')
        outputStr+=(dna_seq_formatted +'\n')
        outputStr+=("//\n")

     outputfile.write(outputStr)
     outputfile.close() 
     rename(output_file_name_tmp, output_file_name)


def format_sequence_origin(dna_seq):
    output=""
    Len =  len(dna_seq)
    for i in range(0, Len):    
       if i==0:
          output+= '%9d' % (i+1)
       if i%10==0:
          output+=' '
       if i!=0  and i%60==0:
          output+= '\n%9d ' % (i+1)
       output +=dna_seq[i]
       i+=1  
    return output

 
def spaces(n):
    space=''
    for  i in  range(0, n):
      space+=' '
    return space
       

def wrap(prefix, start, end, string):
    output=''
    prefixLen = len(prefix)
    i = prefixLen
    output+=prefix
    while i< start:
      output += ' '
      i += 1

    for c in string:
       if i > end:
          output+='\n' 
          i = start
          output+=spaces(i) 
       i+=1
       output+=c 
    return output

    
    

def process_sequence_file(sequence_file_name,  seq_dictionary):
     try:
        sequencefile = open(sequence_file_name, 'r')
     except IOError:
        print "Cannot read file " + sequence_file_name + " !"

     sequence_lines = sequencefile.readlines()
     sequencefile.close()
     fragments= []
     name=""

#     count = 0
     seq_beg_pattern = re.compile(">")
     for line in sequence_lines:
        line = line.strip() 
        if seq_beg_pattern.search(line):
          if len(name) > 0:
             sequence=''.join(fragments)
             seq_dictionary[name]=sequence
             fragments = []
          name=get_sequence_name(line)
        else:
          fragments.append(line)

    #add the final sequence
     if len(name) > 0:
         sequence=''.join(fragments)
         seq_dictionary[name]=sequence
     
#        if count > 1000:
#           sys.exit(0)
#        count=count+1
#     fields = re.split('\t', line)
        #print table[str(fields[0].strip())]
     #print blast_file + ' ' + tax_maps + ' ' + database


usage =  sys.argv[0] + """ -g gff_files -n nucleotide_sequences -p protein_sequences [--out-gbk gbkfile --out-sequin sequinfile --out-ptinput ptinputdir]\n"""

parser = None

def createParser():
    global parser
    epilog = """This script has three functions : (i) The functional and taxonomic annotations created for the individual ORFs are used to create the inputs required by the Pathway-Tools's Pathologic algorithm to build the ePGDBs. The input consists of 4 files that contains functional annotations and sequences with relevant information, this information is used by Pathologic to create the ePGDBs. (ii) It can create a genbank file for the ORFs and their annotations, (iii) An option is added where it can create a sequin file, which is required for sometimes for sequence submission to NCBI data repository, such as trace archive"""
    epilog = re.sub(r'\s+', ' ', epilog)

    parser = optparse.OptionParser(usage=usage, epilog = epilog)

    # Input options

    input_group = optparse.OptionGroup(parser, 'input options')

    input_group.add_option('-g', '--gff', dest='gff_file',
                           metavar='GFF_FILE', 
                           help='GFF files, with annotations,  to convert to genbank format')

    input_group.add_option('-n', '--nucleotide', dest='nucleotide_sequences',
                           metavar='NUCLEOTIDE_SEQUENCE', 
                           help='Nucleotide sequences')

    input_group.add_option('-p', '--protein', dest='protein_sequences',
                           metavar='PROTEIN_SEQUENCE', 
                           help='Protein sequences')


    input_group.add_option('-o', '--output', dest='output',
                           metavar='OUTPUT', help='Genbank file')

    parser.add_option_group(input_group)

    output_options_group =  optparse.OptionGroup(parser, 'Output Options')

    output_options_group.add_option("--out-gbk", dest="gbk_file",  
                     help='option to create a genbank file')
    output_options_group.add_option("--out-sequin", dest="sequin_file", 
                     help='option to create a sequin file')
    output_options_group.add_option("--sequin-params-file", dest="sequin_params_file", 
                     help='NCBI sequin parameters  file')

    output_options_group.add_option("--sequin-tbl2asn", dest="sequin_tbl2asn", 
                     help='the executable for the NCBI tbl2asn')

    output_options_group.add_option("--ncbi-sbt-file", dest="ncbi_sbt_file", 
                     help='the NCBI sbt file location created by the \"Create Submission Template\" form: http://www.ncbi.nlm.nih.gov/WebSub/temp    late.cgi"the sbt file executable for the NCBI tbl2asn')

    output_options_group.add_option("--out-ptinput", dest="ptinput_file", 
                     help='option to create ptools input files')
    parser.add_option_group(output_options_group)


def main(argv, errorlogger = None, runstatslogger = None):
    # Parse options (many!)
    # TODO: Create option groups
    # filtering options


    global parser
    options, args = parser.parse_args(argv)

    if not(options.gff_file or options.nucleotide_sequences or options.protein_sequences or options.output):
      print help
      sys.exit(0)
    
    if not options.gff_file:
       eprintf("ERROR\tGFF file not specified\n")
       errorlogger.printf("ERROR\tGFF file not specified\n")

    if not options.nucleotide_sequences:
       eprintf("ERROR\tNucleotide sequences not specified\n")
       errorlogger.printf("ERROR\tNucleotide sequences not specified\n")

    if not options.protein_sequences:
       parser.error('Protein sequences')

    if not options.gbk_file and not options.sequin_file and not options.ptinput_file:
       #print 'No genbank or sequin or ptools input is specified'
       return (0,'')

    if options.sequin_file and  not options.sequin_params_file:
        parser.error('Cannot create NCBI Sequin input without a parameter file')
     
    #print options

    if not path.exists(options.gff_file):
        print "gff file does not exist"
        eprintf("ERROR\tGFF file %s  not found\n", options.gff_file)
        errorlogger.printf("ERROR\tGFF file %s  not found\n", options.gff_file)
        sys.exit(0)

    if not path.exists(options.nucleotide_sequences):
        print "nucloetide sequences file does not exist"
        sys.exit(0)

    if not path.exists(options.protein_sequences):
        print "protein sequences file does not exist"
        sys.exit(0)

    output_files = {}
    input_files = {}
    if  options.gbk_file:
       output_files['gbk'] = options.gbk_file

    if  options.sequin_file:
       output_files['sequin'] = options.sequin_file
       input_files['sequin_params'] = options.sequin_params_file
       input_files['sequin_fasta'] = options.nucleotide_sequences
       input_files['sequin_tbl2asn'] = options.sequin_tbl2asn
       input_files['sequin_sbt_file'] = options.ncbi_sbt_file


    if  options.ptinput_file:
       output_files['ptinput'] = options.ptinput_file
    

    nucleotide_seq_dict = {}
    process_sequence_file( options.nucleotide_sequences, nucleotide_seq_dict) 
    protein_seq_dict = {}
    process_sequence_file(options.protein_sequences, protein_seq_dict) 
    process_gff_file(options.gff_file, output_files, nucleotide_seq_dict, protein_seq_dict, input_files) 
    #print params['bitscore']

def MetaPathways_create_genbank_ptinput_sequin(argv, errorlogger = None, runstatslogger = None):
    createParser()
    main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger)
    return (0,'')
    
if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])


