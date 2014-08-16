#!/usr/bin/python
"""This script applies parsed BLAST results from a CSV file (produced by
parseblast.pl) to a series of GenBank files.

Oct 26, 2009 by Simon Eng
    Changed the method of outputting statistics.
"""

try:
     import optparse, csv, re, sys
     from os import path
     import logging.handlers
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     sys.exit(3)


def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

def get_translation_table(num):
   translation_table={}
   translation_table['11'] =\
        {
          'GCC' : 'A',
          'AGT' : 'S',
          'TGA' : '*',
          'TGT' : 'C',
          'CGA' : 'R',
          'ATC' : 'I',
          'AAC' : 'N',
          'AGC' : 'S',
          'TAC' : 'Y',
          'ACA' : 'T',
          'TCG' : 'S',
          'CCG' : 'P',
          'CTG' : 'L',
          'GCA' : 'A',
          'GTG' : 'V',
          'AAG' : 'K',
          'GTT' : 'V',
          'CAC' : 'H',
          'AGA' : 'R',
          'ACC' : 'T',
          'CCA' : 'P',
          'TGG' : 'W',
          'CTC' : 'L',
          'CGC' : 'R',
          'TTG' : 'L',
          'TAA' : '*',
          'CAG' : 'Q',
          'ACG' : 'T',
          'AAA' : 'K',
          'ATG' : 'M',
          'GTA' : 'V',
          'TAG' : '*',
          'CTT' : 'L',
          'GGA' : 'G',
          'GTC' : 'V',
          'TGC' : 'C',
          'TCA' : 'S',
          'ATT' : 'I',
          'TAT' : 'Y',
          'AAT' : 'N',
          'ACT' : 'T',
          'CAA' : 'Q',
          'GAC' : 'D',
          'GGT' : 'G',
          'TCC' : 'S',
          'TTT' : 'F',
          'AGG' : 'R',
          'CGT' : 'R',
          'CGG' : 'R',
          'CAT' : 'H',
          'ATA' : 'I',
          'CCC' : 'P',
          'GGG' : 'G',
          'TTA' : 'L',
          'GAG' : 'E',
          'CTA' : 'L',
          'GAT' : 'D',
          'TCT' : 'S',
          'TTC' : 'F',
          'GCG' : 'A',
          'GGC' : 'G',
          'GCT' : 'A',
          'GAA' : 'E',
          'CCT' : 'P'
        }

   return translation_table[str(num)] 

def files_exist( files ):
     for file in files:
        if not path.exists(file):
           print 'Could not read File ' +  file
           print 'Please make sure these sequences are in the \"blastDB\" folder'
           sys.exit(3)
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
     attributes['orf_length'] =  attributes['end'] - attributes['start']

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


def process_gff_file(gff_file_name, output_amino_file_name, output_nuc_file_name, output_gff_file_name,  nucleotide_seq_dict):
     try:
        gfffile = open(gff_file_name, 'r')
     except IOError:
        print "Cannot read file " + gff_file_name + " !"

     sample_name= re.sub('.gff', '', gff_file_name)
     gff_lines = gfffile.readlines()
     gff_beg_pattern = re.compile("^#")
     gfffile.close()
     
     contig_dict={} 
     for line in gff_lines:
        line = line.strip() 
        if gff_beg_pattern.search(line):
          continue
        insert_orf_into_dict(line, contig_dict)

     create_the_gene_IDs(contig_dict, nucleotide_seq_dict)
     create_sequnces(output_amino_file_name, output_nuc_file_name, contig_dict, nucleotide_seq_dict)
     write_gff_file(output_gff_file_name, contig_dict)

def write_gff_file(output_gff_file, contig_dict):
     outputfile = open(output_gff_file, 'w')
     
     for key in contig_dict:
        for elem in contig_dict[key]  :
            fprintf(outputfile,"%s",elem['seqname'])
            fprintf(outputfile,"\t%s",elem['source'])
            fprintf(outputfile,"\t%s",elem['feature'])
            fprintf(outputfile,"\t%s",elem['start'])
            fprintf(outputfile,"\t%s",elem['end'])
            fprintf(outputfile,"\t%s",elem['score'])
            fprintf(outputfile,"\t%s",elem['strand'])
            fprintf(outputfile,"\t%s",elem['frame'])
           
            if 'ID' in elem :
               line = 'ID=' + elem['ID']
               line = line +';' + 'locus_tag=' + elem['ID']
            else:
               line = 'ID= '

            if 'partial' in elem :
               line = line +';' + 'partial=' + elem['partial']

            if 'orf_length' in elem :
               line = line +';' + 'orf_length=' + str(elem['orf_length'])

            if 'contig_length' in elem :
               line = line +';' + 'contig_length=' + str(elem['contig_length'])


            fprintf(outputfile,"\t%s",line)
            fprintf(outputfile,"\n")

     outputfile.close()

def create_the_gene_IDs(contig_dict, nucleotide_seq_dict):
     for key in contig_dict:
        count = 0 
        for elem in contig_dict[key]  :
           name = elem['seqname'] + "_%d" % (count) 
           contig_dict[key][count]['ID'] = name
           contig_dict[key][count]['contig_length'] = len(nucleotide_seq_dict[key])
           count += 1

def create_sequnces(output_amino_file_name, output_nuc_file_name, contig_dict, nucleotide_seq_dict): 

     translation_table = get_translation_table(11)
     #print translation_table
     #print contig_dict
     
     aa_outputfile = open(output_amino_file_name, 'w')
     nucl_outputfile = open(output_nuc_file_name, 'w')

     for key in contig_dict:
        for elem in contig_dict[key]  :
          nuc_orf_sequence, aa_orf_sequence=get_amino_acid_sequence(nucleotide_seq_dict, translation_table,  elem['seqname'], elem['start'], elem['end'], elem['strand'])  
          name = '>' + elem['ID'] 
          if len(aa_orf_sequence) > 0 :
             fprintf(aa_outputfile, "%s\n", name)
             fprintf(aa_outputfile, "%s\n", aa_orf_sequence)

          if len(nuc_orf_sequence) > 0 :
             fprintf(nucl_outputfile, "%s\n", name)
             fprintf(nucl_outputfile, "%s\n", nuc_orf_sequence)

     aa_outputfile.close() 
     nucl_outputfile.close() 

def get_amino_acid_sequence(nucleotide_seq_dict, translation_table, seqname, start, end, strand):  

    try:
       if strand=='-':
#          nuc_orf_sequence = reverse(nucleotide_seq_dict[seqname])
#          nuc_orf_sequence = nuc_orf_sequence[start-1:end-1].upper()

          nuc_orf_sequence = nucleotide_seq_dict[seqname]

          nuc_orf_sequence = nuc_orf_sequence[start-1:end].upper()
          a = len(nuc_orf_sequence)
          nuc_orf_sequence = reverse(nuc_orf_sequence)
          b = len(nuc_orf_sequence)
          nuc_orf_sequence = complement(nuc_orf_sequence)
          c = len(nuc_orf_sequence)
          #print str(a) + ' ' + str(b) + ' ' + str(c)
       else:
          nuc_orf_sequence = nucleotide_seq_dict[seqname][start-1:end-1].upper()

       
       aa_orf_sequence = convert_to_amino_acid(translation_table, nuc_orf_sequence, strand) 
       return (nuc_orf_sequence, aa_orf_sequence)
    except:
       return ("XXXXX", "YYYYY")

def convert_to_amino_acid(translation_table, nuc_orf_sequence, strand) :
    #print nuc_orf_sequence
    letters = list(nuc_orf_sequence)
    aa_sequence =''
    codon = ''
    for i in range(len(nuc_orf_sequence)):
        codon = codon + letters[i] 
        if i%3==2 :
           if (codon in translation_table)  and not (translation_table[codon]=='*') :
              aa_sequence = aa_sequence + translation_table[codon]
           codon=''

    return aa_sequence



def reverse(nuc_orf_sequence):
    letters = list(nuc_orf_sequence)
    i = len(nuc_orf_sequence)-1
    reverselist = []
    while i >= 0:
        reverselist.append(letters[i])
        i=i-1

    result_string ="".join(reverselist)
    return result_string

    
def complement(nuc_orf_sequence):
    letters = list(nuc_orf_sequence)
    length = len(nuc_orf_sequence)
    complementlist = []
    i = 0
    while i < length:
        if letters[i] == 'A':
           complementlist.append('T')
        elif letters[i] == 'T':
           complementlist.append('A')
        elif letters[i] == 'C':
           complementlist.append('G')
        elif letters[i] == 'G':
           complementlist.append('C')
        elif letters[i] == 'N':
           complementlist.append('N')
        else:
           complementlist.append('*')
        i=i+1

    
    result_string ="".join(complementlist)
    return result_string




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
     if len(name) > 0:
        sequence=''.join(fragments)
        seq_dictionary[name]=sequence



usage = sys.argv[0] + """ -g gff_files -n nucleotide_sequences --output-nuc output_nuc --output_amino  output_amino  --output_gff output_gff_file"""

parser = None
def createParser():
    # Parse options (many!)
    # TODO: Create option groups
    global parser

    epilog = """The GFF files produced by the gene prediction algorithm, say prodigal, in the ORF Prediction step are used to translate into the amino acid sequences and gff files that contains more detailed information such as orf_length, strand, start and end positions, contig_length, etc. In the pipeline the resulting file is in the orf_prediction folder with the extension .faa """

    epilog = re.sub(r'\s+', ' ', epilog)
    parser = optparse.OptionParser(usage = usage, epilog = epilog)
    # Input options

    input_group = optparse.OptionGroup(parser, 'input options')

    input_group.add_option('-g', '--gff', dest='gff_file',
                           metavar='GFF_FILE', 
                           help='GFF files to convert to genbank format')

    input_group.add_option('-n', '--nucleotide', dest='nucleotide_sequences',
                           metavar='NUCLEOTIDE_SEQUENCE', 
                           help='nucleotide sequences input file ')

    input_group.add_option('--output_amino', dest='output_amino',
                           metavar='OUTPUT', help='nmino acid output file')

    input_group.add_option('--output_nuc', dest='output_nuc',
                           metavar='OUTPUT', help='nucleotide Sequence file')

    input_group.add_option('--output_gff', dest='output_gff',
                           metavar='OUTPUT', help='output file in gff format')


    parser.add_option_group(input_group)


def main(argv, errorlogger = None, runstatslogger = None):
    # filtering options
    global parser
    options, args = parser.parse_args(argv)

    if not(options.gff_file or options.nucleotide_sequences or options.output_amino or  options.output_nuc  or options.output_gff):
      print help
      sys.exit(0)
    
    if not options.gff_file:
       parser.error('No gff files are specified')

    if not options.nucleotide_sequences:
       parser.error('Nucleotide sequences')

    if not options.output_amino:
       parser.error('Output anino acid file must be specified')

    if not options.output_nuc:
       parser.error('Output nucloetide sequences file must be specified')

    if not options.output_gff:
       parser.error('Output gff file must be specified')
    #print options

    if not path.exists(options.gff_file):
        print "gff file does not exist"
        sys.exit(0)

    if not path.exists(options.nucleotide_sequences):
        print "nucloetide sequences file does not exist"
        sys.exit(0)

    nucleotide_seq_dict = {}
    process_sequence_file( options.nucleotide_sequences, nucleotide_seq_dict) 
    process_gff_file(options.gff_file, options.output_amino, options.output_nuc, options.output_gff, nucleotide_seq_dict) 


    #print params['bitscore']

def MetaPathways_create_amino_sequences(argv, errorlogger = None, runstatslogger = None):
    createParser()
    main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger)
    return (0,'')
    
if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

