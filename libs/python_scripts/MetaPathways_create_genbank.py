#!/usr/bin/python
__author__ = 'nielshanson'

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
usage =  sys.argv[0] + """ -g gff_files -n nucleotide_sequences -p protein_sequences [--out-gbk gbkfile --out-sequin sequinfile --out-ptinput ptinputdir]\n"""
parser = None


def createParser(argv, errorlogger = None, runstatslogger = None):
    global parser
    epilog = """This script has one function: It can create a genbank file for the ORFs and their annotations"""
    epilog = re.sub(r'\s+', ' ', epilog)

    parser = optparse.OptionParser(usage=usage, epilog = epilog)

    # Input options
    genbank_group = optparse.OptionGroup(parser, 'genbank options')

    genbank_group.add_option('-g', '--gff', dest='gff_file',
                             metavar='GFF_FILE',
                             help='GFF files, with annotations, to convert to genbank format')
    genbank_group.add_option('-n', '--nucleotide', dest='nucleotide_sequences',
                             metavar='NUCLEOTIDE_SEQUENCE',
                             help='Nucleotide sequences')
    genbank_group.add_option('-p', '--protein', dest='protein_sequences',
                             metavar='PROTEIN_SEQUENCE',
                             help='Protein sequences')
    genbank_group.add_option('-o', '--output', dest='gbk_file',
                             metavar='OUTPUT', help='Genbank file')
    parser.add_option_group(genbank_group)

def insertAttribute(attributes, attribStr):
    rawfields = re.split('=', attribStr)
    if len(rawfields) == 2:
        attributes[rawfields[0].strip().lower()] = rawfields[1].strip()

def splitAttributes(str, attributes):
    rawattributes = re.split(';', str)
    for attribStr in rawattributes:
        insertAttribute(attributes, attribStr)

    return attributes

def getSequenceFile(line):
    fields = re.split(' ', line)
    name = re.sub('>','',fields[0])

    return name

def processSequenceFile(sequence_file_name,  seq_dictionary):
    try:
        sequence_file = open(sequence_file_name, 'r')
    except IOError:
        print "Cannot read file " + sequence_file_name + " !"

    sequence_lines = sequence_file.readlines()
    sequence_file.close()
    fragments = []
    name = ""

    seq_beg_pattern = re.compile(">")
    for line in sequence_lines:
        line = line.strip()
        if seq_beg_pattern.search(line):
            if len(name) > 0:
                sequence=''.join(fragments)
                seq_dictionary[name]=sequence
                fragments = []
            name = getSequenceFile(line)
        else:
            fragments.append(line)

            # add the final sequence
    if len(name) > 0:
        sequence=''.join(fragments)
        seq_dictionary[name]=sequence


def formatSequenceOrigin(dna_seq):
    output=""
    l =  len(dna_seq)
    for i in range(0, l):
        if i == 0:
            output+= '%9d' % (i+1)
        if (i % 10) == 0:
            output+=' '
        if (i != 0) and ((i % 60) == 0):
            output+= '\n%9d ' % (i+1)
        output += dna_seq[i]

    return output


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
            output+= " "*i
        i+=1
        output+=c

    return output


def writeGbkFile(output_file_name, contig_dict, sample_name, nucleotide_seq_dict, protein_seq_dict):
    """This function creates the genbank file from the gff, protein and nucleotide sequences."""
    date = genbankDate()
    output_file_name_tmp = output_file_name + ".tmp"
    outputfile = open(output_file_name_tmp, 'w')

    count = 0
    outputStr = ""
    for key in contig_dict:
        first = True
        if count % 10000 == 0:
            outputfile.write(outputStr)
            outputStr = ""
        count+=1

        for attrib in contig_dict[key]:
            id  = attrib['id']
            try:
                protein_seq = protein_seq_dict[id]
            except:
                protein_seq = ""

            definition = sample_name
            accession = '.'
            version = '.' + " "*10 + "GI:."
            dblink = sample_name
            keywords = '.'
            source = sample_name
            organism = sample_name
            dna_seq_formatted = ""
            if first:
                first = False
                try:
                    dna_seq =  nucleotide_seq_dict[key]
                    dna_seq_formatted =  formatSequenceOrigin(dna_seq)
                    dna_length = len(dna_seq)
                    sourceStr = "1.." + str(dna_length)
                except:
                    dna_seq = ""
                    dna_seq_formatted = ""
                    dna_length = 0
                    sourceStr = "0..0"

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


def processGffFile(gff_file_name, output_filenames, nucleotide_seq_dict, protein_seq_dict, input_filenames):
    # print output_filenames
    try:
        gff_file = open(gff_file_name, 'r')
    except IOError:
        print "Cannot read file " + gff_file_name + " !"

    sample_name= re.sub(".annot.gff", '', gff_file_name)
    sample_name= re.sub(r'.*[/\\]', '', sample_name)

    gff_lines = gff_file.readlines()
    gff_beg_pattern = re.compile("^#")
    gff_file.close()

    contig_dict={}
    for line in gff_lines:
        line = line.strip()
        if gff_beg_pattern.search(line):
            continue
        insertOrfIntoDict(line, contig_dict)

    if "gbk" in output_filenames:
        writeGbkFile(output_filenames['gbk'], contig_dict, sample_name, nucleotide_seq_dict, protein_seq_dict)


def insertOrfIntoDict(line, contig_dict):
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
    splitAttributes(fields[8], attributes)

    if not fields[0] in contig_dict :
        contig_dict[fields[0]] = []

    contig_dict[fields[0]].append(attributes)

def main(argv, extra_command = None, errorlogger = None, runstatslogger =None):
    global parser
    options, args = parser.parse_args(argv)

    if not(options.gff_file or options.nucleotide_sequences or options.protein_sequences or options.gbk_file):
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

    if not path.exists(options.gff_file):
        print "gff file does not exist"
        eprintf("ERROR\tGFF file %s  not found\n", options.gff_file)
        errorlogger.printf("ERROR\tGFF file %s  not found\n", options.gff_file)
        sys.exit(0)

    if not path.exists(options.nucleotide_sequences):
        print "nucleotide sequences file does not exist"
        sys.exit(0)

    if not path.exists(options.protein_sequences):
        print "protein sequences file does not exist"
        sys.exit(0)

    output_files = {}
    input_files = {}
    if  options.gbk_file:
        output_files['gbk'] = options.gbk_file

    nucleotide_seq_dict = {}
    processSequenceFile( options.nucleotide_sequences, nucleotide_seq_dict)
    protein_seq_dict = {}
    processSequenceFile(options.protein_sequences, protein_seq_dict)
    processGffFile(options.gff_file, output_files, nucleotide_seq_dict, protein_seq_dict, input_files)

    return (0,'')


def MetaPathways_create_genbank(argv, extra_command = None, errorlogger = None, runstatslogger =None):
    if errorlogger:
        errorlogger.write("#STEP\tGENBANK_FILE\n")
    createParser(argv)
    main(argv, extra_command = extra_command, errorlogger = errorlogger, runstatslogger = runstatslogger)
    return (0,'')

if __name__ == '__main__':
    createParser(sys.argv[1:])
    main(sys.argv[1:])