#!/usr/bin/python
"""This script BLASTs the preprocessed input agains the Taxonomy databases"""

try:
    import optparse, sys, re, csv, traceback
    from os import path
    import logging.handlers

    from libs.python_modules.utils.sysutil import pathDelim
    from libs.python_modules.utils.metapathways_utils  import fprintf, printf, eprintf,  exit_process
    from libs.python_modules.utils.sysutil import getstatusoutput
except:
    print """ Could not load some user defined  module functions"""
    print """ Make sure your typed 'source MetaPathwaysrc'"""
    print """ """
    print traceback.print_exc(10)
    sys.exit(3)

PATHDELIM=pathDelim()



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

def getFormat(dbstring):
    dbstring = dbstring.rstrip()
    dbstring = dbstring.lstrip()
    dbstring = dbstring.lower()

    # silva, greengenes spellings should be lowercase always
    format = 0
    if re.search(r'silva',dbstring, re.I):
        format = 1

    if re.search(r'greengenes',dbstring, re.I):
        format = 2
    return format

#sequences with no seperate taxonomic name gets the sequences name
def parse_Format_0(line):
    fields = re.split(' ', line)
    if len( fields ) ==1:
        name = fields[0].replace('>','')
        taxonomy = name
    else:
        return( None, None )
    return( name.strip(), taxonomy.strip() )

#silva 
def parse_Format_1(line):
    fields = re.split(' ', line)
    if len( fields ) >= 2:
        name = fields[0].replace('>','')
        fields = fields[1:]
        taxonomy = " "
        taxonomy = taxonomy.join(fields)
    else:
        return( None, None )

    return( name.strip(), taxonomy.strip() )

def parse_Format_2(line):
    fields = re.split(' ', line)
    names = []
    if len( fields ) >= 2:
        name = fields[0].replace('>','')
        for field in fields[1:]:
            field =field.strip()
            if re.match("[a-z]_[a-zA-Z0-9-_]*;",field):
                names.append(field)
        taxonomy = ""
        taxonomy = taxonomy.join(names)
    else:
        return( None, None )

    return( name.strip(), taxonomy.strip() )


def getName_and_Taxonomy(line, format=0):

    if format==0:
        name, taxonomy = parse_Format_0(line)
    elif format==1:
        name, taxonomy = parse_Format_1(line)
    #   print name + "  " + taxonomy
    elif format==2:
        (name, taxonomy) = parse_Format_2(line)
    else:
        return( None, None )

    return(name, taxonomy )

START_PATTERN = re.compile(r'^>')

def read_select_fasta_sequences(candidates, records, input_file_name):
    input_file = open(input_file_name, 'r')
    line = input_file.readline()
    while line:
        if START_PATTERN.search(line):
            name=line.strip()
            name = re.sub(r'>', '',name)
            if name and name in candidates:
                records[name]=""
        else:
            sequence = line.strip()
            if sequence and  name in candidates:
                records[name] += sequence
        line = input_file.readline()

    input_file.close()

def write_selected_sequences(selected_sequences, output_file_name):
    output_file = open(output_file_name, 'w')
    for read in selected_sequences:
        fprintf(output_file, ">%s\n", read)
        fprintf(output_file,"%s\n", selected_sequences[read])

    output_file.close()



def append_taxonomic_information(databaseSequences, table, params):
    try:
        tax_seqs = open(databaseSequences, 'r')
    except IOError:
        print "Cannot read file " + databaseSequences + " !"
        sys.exit(3)


    format = getFormat(databaseSequences)

    taxmapLines = tax_seqs.readlines()
    tax_seqs.close()

    taxMapping={}
    for line in taxmapLines:
        if not re.match('>', line):
            continue
        line = line.strip()
        name, taxonomy = getName_and_Taxonomy(line, format)
        if name:
            taxMapping[name] = taxonomy

    for key in  table:
        key = str(key)
        if int(table[key][5] - table[key][4] ) > params['length'] and table[key][0] > params['similarity'] and table[key][1] < params['evalue'] and table[key][2] > params['bitscore']:
            if  table[key][3] in taxMapping:
                table[key].append(taxMapping[table[key][3]])
            else:
                table[key].append('-')
        else:
            table[key].append('-')



def process_blastout_file(blast_file, database, table, errorlogger = None):
    try:
        blastfile = open(blast_file, 'r')
    except IOError:
        eprintf("ERROR : Cannot write read file " + blast_file + " !" )
        if errorlogger!=None:
            errorlogger.write("STATS_rRNA\tERROR\tCannot write read blast output file " + blast_file + " for database " + database )
        exit_process()

    blastLines = blastfile.readlines()
    blastfile.close()

    for line in blastLines:
        line = line.strip()
        fields = re.split('\t', line)
        if len(fields) < 12:
            continue
        fields[0] =  str(fields[0].strip())
        fields[1] =  str(fields[1].strip())
        fields[2] =  float(fields[2].strip())
        fields[6] =  int(fields[6].strip())
        fields[7] =  int(fields[7].strip())
        fields[10] = float(fields[10].strip())
        fields[11] = float(fields[11].strip())
        table[str(fields[0].strip())] = [fields[2], fields [10], fields[11], fields[1], fields[6], fields[7]]

    #     print table
    #     sys.exit(0)

    #print blast_file + ' ' + tax_maps + ' ' + database

usage = sys.argv[0] + """ -i x.blastout [y.blastout] -d  xdatabase [ydatabase]  -m xtax_maps [ ytax_maps ] -o outputfile -b n -e 0.ddd -s N """
parser = None
def createParser():
    global parser

    epilog = """The input nucleotide sequences are BLASTed against the selected rRNA databases such as  SSU or LSU Silva sequence databases and SSU Greengene database. The hits with high bit scores are flagged as rRNA and the resulting taxonomy from the databases are assigned.  The results from this step are put in the results/rRNA folder, with one tsv file for each rRNA database."""
    epilog = re.sub(r'\s+',' ', epilog)

    parser = optparse.OptionParser(usage=usage, epilog = epilog)

    # Input options

    input_group = optparse.OptionGroup(parser, 'input options')

    input_group.add_option('-i', '--blastout', dest='blast_files',
                           metavar='BLAST_FILE', action='append', default=[],
                           help='BLAST  output file  from blasting of nucleotide sequeces against Taxonomic databases')

    input_group.add_option('-d', '--databases', dest='tax_databases',
                           metavar='TAX_DATABASE', action='append', default=[],
                           help='Taxonomic databases')


    input_group.add_option('-o', '--output', dest='output',
                           metavar='OUTPUTE', help='Taxonomic databases')


    input_group.add_option('-f', '--fasta', dest='fasta',
                           metavar='NUC_SEQUENCES', help='The nucleotide sequences')



    parser.add_option_group(input_group)

    # filtering options
    filtering_options = optparse.OptionGroup(parser, 'filteroptions')

    filtering_options.add_option('-b', '--bit_score', dest='bitscore', metavar='BITSCORE',
                                 default=0, help="bit score cutoff")

    filtering_options.add_option('-e', '--e-value', metavar= 'evalue',
                                 dest='evalue', default=1e-6,
                                 help='e-value cut off for the blast hits')

    filtering_options.add_option('-s', '--similarity', metavar='similarity',
                                 dest='similarity', default = 40,
                                 help='% similarity cut off for the blast hits')

    filtering_options.add_option('-l', '--length', metavar='length',
                                 dest='length', default = 180,
                                 help='length cut off for the blast hits')

    parser.add_option_group(filtering_options)

def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser

    options, args = parser.parse_args(argv)

    if not len(options.blast_files):
        parser.error('At least one taxonomic BLAST output is required')

    if runBlastCommandrRNA(runcommand = runcommand) !=0:
        if errorlogger:
            errorlogger.write("ERROR: Failed to BLAST the sequences against database %s : "  %(options.tax_databases[0]) )
            errorlogger.write("     : " + runcommand)
        exit_process("ERROR: Failed to BLAST the sequences against database %s : "  %(options.tax_databases[0]) + \
                     "     : " + runcommand)

    if not ( len(options.tax_databases) == len( options.blast_files) ):
        parser.error('Number of taxonomic databases and BLAST outputs should be the same')

    if not options.output:
        parser.error('Output file must be specified')
    # Incredible sanity check

    if not files_exist(options.blast_files):
        sys.exit(0)

    if not files_exist( options.tax_databases):
        sys.exit(0)

    params = {'length': int(options.length), 'similarity': float(options.similarity), 'evalue':float(options.evalue), 'bitscore':float(options.bitscore) }
    #print params['bitscore']
    table={}
    for x in range(0, len(options.blast_files)):
        table[options.tax_databases[x]]={}
        process_blastout_file(options.blast_files[x], options.tax_databases[x],table[options.tax_databases[x]], errorlogger = errorlogger)


    priority = 7000
    reads = {}
    for x in range(0, len(options.blast_files)):
        append_taxonomic_information(options.tax_databases[x], table[options.tax_databases[x]],  params)
        for key in table[options.tax_databases[x]]:
            if len(table[options.tax_databases[x]][key][6]) > 1:
                reads[key] = True

        dbname  =  re.sub(r'^.*' + PATHDELIM, '', options.tax_databases[x])
        runstatslogger.write("%s\tTaxonomic hits in %s\t%s\n" %(str(priority),  dbname,  str(len(reads))))
        priority += 1
    outputfile = open(options.output, 'w')
    fprintf(outputfile, "#Similarity cutoff :\t" +  str(params['similarity']) +'\n')
    fprintf(outputfile, "#Length cutoff :\t" +  str(params['length']) +'\n')
    fprintf(outputfile, "#Evalue cutoff :\t" +  str(params['evalue']) +'\n')
    fprintf(outputfile, "#Bit score cutoff :\t" +  str(params['bitscore']) +'\n')
    fprintf(outputfile, "#Number of rRNA sequences detected:\t" +  str(len(reads)) +'\n\n')


    for x in range(0, len(options.tax_databases)):
        #  printf('\t%s\t\t\t', re.sub(r'^.*/','', options.tax_databases[x]))
        fprintf(outputfile, '\t%s\t\t\t', re.sub(r'^.*' + PATHDELIM, '', options.tax_databases[x]))
    #printf('\n')
    fprintf(outputfile,'\n')


    #printf('%s', 'read')
    for x in range(0, len(options.blast_files)):
        fprintf(outputfile, '%s\t%s\t%s\t%s\t%s\t%s\t%s', 'sequence', 'start', 'end', 'similarity', 'evalue', 'bitscore', 'taxonomy')
    fprintf(outputfile,'\n')

    for read in reads:
        #printf('%s', read)
        fprintf(outputfile,'%s', read)
        for x in range(0, len(options.blast_files)):
            if read in table[options.tax_databases[x]]:
                fprintf(outputfile, '\t%s\t%s\t%s\t%s\t%s\t%s', str(table[options.tax_databases[x]][read][4]), str(table[options.tax_databases[x]][read][5]), str(table[options.tax_databases[x]][read][0]),str(table[options.tax_databases[x]][read][1]),str(table[options.tax_databases[x]][read][2]), str(table[options.tax_databases[x]][read][6]))
            else:
                fprintf(outputfile, '\t-\t-\t-\t-\t-\t-')
        fprintf(outputfile,'\n')
    outputfile.close()

    # collect the exact reads 
    database_hits = {}
    for read in reads:
        for x in range(0, len(options.blast_files)):
            if read in table[options.tax_databases[x]]:
                database_hits[read] = [ table[options.tax_databases[x]][read][4], table[options.tax_databases[x]][read][5]]

    # pick the hits, trim them according to the match and write them
    if options.fasta:
        selected_sequences={}
        read_select_fasta_sequences(database_hits, selected_sequences, options.fasta)
        for read in database_hits:
            selected_sequences[read] = selected_sequences[read][database_hits[read][0]:database_hits[read][1]]
        write_selected_sequences(selected_sequences, options.output +'.fasta')

def runBlastCommandrRNA(runcommand = None):
    if runcommand == None:
        return False
    result = getstatusoutput(runcommand)

    return result[0]



def MetaPathways_rRNA_stats_calculator(argv, extra_command = None, errorlogger = None, runstatslogger =None):
    if errorlogger != None:
        errorlogger.write("#STEP\tSTATS_rRNA\n")
    createParser()
    main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (0,'')

if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

