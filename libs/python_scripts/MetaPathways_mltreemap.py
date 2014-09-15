#!/usr/bin/python
##############################
##
## MLTreeMap TSN v. 0.0
##
##############################

try:
    import sys
    import os
    from os import path, _exit
    import shutil
    import re
    import glob
    import subprocess
    from optparse import OptionParser
    import time
except:
    print """ Could not load some user defined  module functions"""
    print """ """
    print traceback.print_exc(10)
    sys.exit(3)


def os_type():
    """Return the operating system of the user."""
    x = sys.platform
    if x:

        hits = re.search(r'darwin', x, re.I)
        if hits :
          return 'mac'
     
        hits = re.search(r'win', x, re.I)
        if hits :
          return 'win'

        hits = re.search(r'linux', x, re.I)
        if hits:
          return 'linux'


def pathDelim():
     """Return the path deliminator based on the operating system of the user."""
     ostype = os_type()
     if ostype == 'win':
         return "\\"
 
     if ostype in ['linux', 'mac']:
         return "/"

PATHDELIM =  str(pathDelim())

class Autovivify(dict):
    """In cases of Autovivify objects, enable the referencing of variables (and sub-variables) without explicitly declaring those variables beforehand."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


# TK class BookKeeping: we will implement a singleton class for bookkeeping
     

parser =  None
def createParser(): 
    global parser
    epilog = "MLTreeMap step searches for a set of 40 universally present COGs as marker genes and taxonomically identifies some sequences in the Tree of Life."
    epilog = re.sub(r'\s+', ' ', epilog)
    """Returns the parser to interpret user options."""

    parser = OptionParser(description=epilog)
    parser.add_option('-i', '--input', dest='input', help='your sequence input file')
    parser.add_option('-b', '--bootstraps', dest='bootstraps', default=0, type=int, help='the number of Bootstrap replicates')
    parser.add_option('-c', '--cluster', dest='cluster', default=0, choices=[0,'s'], help='use a computer cluster? (0 = no cluster; s = sun grid)')
    parser.add_option('-f', '--phylogeny', dest='phylogeny', default='v', choices=['v','p'], help='RAxML algorithm (v = Maximum Likelihood; p = Maximum Parsimony)')
    parser.add_option('-g', '--gblocks', dest='gblocks', default=50, type=int, help='minimal sequence length after Gblocks')
    parser.add_option('-l', '--filelength', dest='filelength', default=2000, type=int, help='long input files will be splint into files containing L sequences each')
    parser.add_option('-m', '--memory', dest='memory',  default=0, type=int, help='minimum memory on a sungrid_cluster in GB')
    parser.add_option('-o', '--output', dest='output', default='output/', help='output directory')
    parser.add_option('-s', '--bitscore', dest='bitscore',  default=60, type=int, help='minimum bitscore for the blast hits')
    parser.add_option('-t', '--reftree', dest='reftree', default='p', choices=['p','g','i'], help='phylogenetic reference tree (p = MLTreeMap reference tree; g = GEBA reference tree; i = fungi tree)')
    parser.add_option('-r', '--reftype', dest='reftype', default='n', choices=['a','n'], help='the type of input sequences (a = Amino Acid; n = Nucleotide)')
    parser.add_option('-e', '--executables', dest='executables', default='sub_binaries', help='locations of executables (e.g. blastx, Gblocks, etc.)')
    parser.add_option('-x', '--mltreemap', dest='mltreemap', default = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) , help='location of MLTreeMap resources (default: directory of mltreemap-improved.py)')
    parser.add_option('-T', '--num_threads', dest='num_threads', default = None , help='specifies the number of CPU threads to use in raxml and blast (default: 1)')
    parser.add_option('-d', '--delete', dest='delete',  default = None, help='the sections of files to be deleted, as separated by colons (1 = Sequence Files; 2 = BLAST Results; 3 = Genewise Results; 4 = hmmalign and Gblocks Results; 5 = Unparsed RAxML Results')
    parser.add_option( '--overwrite', dest='overwrite',  action='store_true',  default = False,  help='overwrites previously processed output folders')
    parser.add_option( '-v', '--verbose', dest='verbose', action='store_true',  default = False,  help='maintains the various_outputs and final_RaXML_outputs directories containing intermediate files')
    parser.add_option( '--showmessages', dest='showmessages', action='store_true',  default = False,  help='shows the detailed  messages')
    parser.add_option( '--mltreemap_data', dest='mltreemap_data',   default = '',  help='The folder where the mltreemap data')
    return parser

def checkParserArguments(options):
    """Returns 'args', a summary of MLTreeMap settings."""

    # Ensure files contain more than 0 sequences

    if options.filelength <= 0:
        sys.exit('Input files require a positive number of sequences!')

    # Set the reference data file prefix and the reference tree name
    if options.reftree == 'g':
        setattr(options, 'reference_data_prefix', 'geba_')
        setattr(options,'reference_tree', 'geba.tree')
    elif options.reftree == 'i':
        setattr(options,'reference_data_prefix','fungi_')
        setattr(options, 'reference_tree', 'fungitr_tree.txt')
    else:
        setattr(options, 'reference_data_prefix', '')
        setattr(options, 'reference_tree','MLTreeMap_reference.tree')

    return options


def removePreviousOutput(args):
    """
    Prompts the user to determine how to deal with a pre-existing output directory.
    Returns an updated version of 'args', a summary of MLTreeMap settings.
    """

    # Add (or replace a trailing (back)slash with) the PATHDELIM to the end of the output directory
    while re.search(r'/\Z', args.output) or re.search(r'\\\Z', args.output):
        args.output = args.output[:-1]
    args.output += PATHDELIM

    # delete previous output folders by force
    if args.overwrite:
       if path.exists(args.output):
          shutil.rmtree(args.output)

    # Prompt the user to deal with the pre-existing output directory
    while os.path.isdir(args.output):
        print('WARNING: Your output directory "' + args.output + '" already exists!')
        print('Overwrite [1], quit [2], or change directory [3]?')
        answer = raw_input()
        answer = int(answer)

        while not answer == 1 and not answer == 2 and not answer == 3:
            answer = raw_input('Invalid input. Please choose 1, 2, or 3.\n')
            answer = int(answer)
        if answer == 1:
            print('Do you really want to overwrite the old output directory?')
            print('All data in it will be lost!')
            answer2 = raw_input('Yes [y] or no [n]?\n')
            while not answer2 == 'y' and not answer2 == 'n':
                answer2 = raw_input('Invalid input. Please choose y or n.\n')
            if answer2 == 'y':
                shutil.rmtree(args.output)
            else:
                sys.exit('Exit MLTreeMap\n')
        elif answer == 2:
            sys.exit('Exit MLTreeMap\n')
        else:
            args.output = raw_input('Please enter the path to the new directory.\n')
    
    # Create the output directories
    args.output_dir_var = args.output + 'various_outputs' + PATHDELIM
    args.output_dir_raxml = args.output + 'final_RAxML_outputs' + PATHDELIM
    args.output_dir_final = args.output + 'final_outputs' + PATHDELIM
    os.makedirs(args.output)
    os.mkdir(args.output_dir_var)
    os.mkdir(args.output_dir_raxml)
    os.mkdir(args.output_dir_final)

    return args


def createCogList(args):
    """Returns an Autovivification of the COGs in the MLTreeMap COG list file, and a list of output text precursors based on the analysis type."""
    
    cog_list = Autovivify()
    text_of_analysis_type = Autovivify()
    alignment_set = args.reftree
    kind_of_cog = ''

    # For each line in the COG list file...    
    cogInputList = open( args.mltreemap_data + PATHDELIM + \
                         'data' + PATHDELIM + 'tree_data' + PATHDELIM + 'cog_list.txt', 'r')
    cogList = [ x.strip() for x in cogInputList.readlines() ] 
    for cogInput in cogList:

        # Get the kind of COG if cogInput is a header line        
        if (re.match(r'\A#(.*)', cogInput)):
            kind_of_cog = re.match(r'\A#(.*)', cogInput).group(1)
            continue

        # Add data to COG list based on the kind of COG it is
        if (kind_of_cog == 'phylogenetic_cogs'):
            cog_list[kind_of_cog][cogInput] = alignment_set
            cog_list['all_cogs'][cogInput] = alignment_set
            text_inset = ''
            if (alignment_set == 'g'):
                text_inset = ' based on the GEBA reference'
            if (alignment_set == 'i'):
                text_inset = ' focusing only on fungi'
            text_of_analysis_type[alignment_set] = 'Phylogenetic analysis' + text_inset + ':'
        elif (kind_of_cog == 'phylogenetic_rRNA_cogs'):
            cog, denominator, text = cogInput.split('\t')
            cog_list[kind_of_cog][cog] = denominator
            cog_list['all_cogs'][cog] = denominator
            text_of_analysis_type[denominator] = 'Phylogenetic analysis, ' + text + ':'
        elif (kind_of_cog == 'functional_cogs'):
            cog, denominator, text = cogInput.split('\t')
            cog_list[kind_of_cog][cog] = denominator
            cog_list['all_cogs'][cog] = denominator
            text_of_analysis_type[denominator] = 'Functional analysis, ' + text + ':'

    # Close the COG list file
    cogInputList.close()

    return (cog_list, text_of_analysis_type)

def calculate_overlap(info):
    """Returns the overlap length of the base and the check sequences."""
    overlap = 0
    base_start = info['base']['start']
    base_end = info['base']['end']
    check_start = info['check']['start']
    check_end = info['check']['end']

    # Calculate the overlap based on the relative positioning of the base and check sequences
    if base_start <= check_start and check_start <= base_end and base_end <= check_end:
        # Base     ----
        # Check    -------
        overlap = base_end - check_start
    elif base_start <= check_start and check_end <= base_end:
        # Base     --------
        # Check    --
        overlap = check_end - check_start
    elif check_start <= base_start and base_start <= check_end and check_end <= base_end:
        # Base     -----
        # Check    -----
        overlap = check_end - base_start
    elif check_start <= base_start and base_end <= check_end:
        # Base     --
        # Check    --------
        overlap = base_end - base_start

    return overlap 


def splitFastaInput(args):
    """
    Splits the input file into multiple files, each containing a maximum 
    number of sequences as specified by the user.
    Ensures each sequence and sequence name is valid.
    Returns a list of the files produced from the input file.
    """

    # Confirm input file is a fasta file
    input = open(args.input, 'r')
    if (not input.read(1) == '>'):
        sys.exit('ERROR: Your file does not appear to be a proper FASTA file!\n')

    # Unread the '>' to prevent problems later
    input.seek(-1,1)

    # Determine the output file names and open the output files
    if (re.match(r'\A.*\/(.*)', args.input)):
        inputFileName = re.match(r'\A.*\/(.*)', args.input).group(1)
    else:
        inputFileName = args.input
    outputSplit = open(args.output_dir_var + inputFileName + '_0.txt', 'w')
    outputFormatted = open(args.output_dir_var + inputFileName + '_formatted.txt', 'w')
    args.formatted_input_file = args.output_dir_var + inputFileName + '_formatted.txt'
    countFiles = 0
    countSequences = 0

    # Iterate through the input file...
    countTotal = 0
    countNucleotides = 0
    countXN = 0
    countUndef = 0
    splitFiles = []

    for line in input:
        # If the line is a sequence name...
        if (re.search('\A>', line)):
            countSequences += 1

            # Replace all non a-z, A-Z, 0-9, or . characters with a _
            # Then replace the initial _ with a >
            line = re.sub(r'[^a-zA-Z0-9.\r\n]', '_', line)
            line = re.sub(r'\A_', '>', line)
            
            # Because RAxML can only work with file names having length <= 125,
            # Ensure that the sequence name length is <= 100
            if (line.__len__() > 100):
                line = line[0:100]

            # Split the file if countSequences > the number of sequences per file specified by the user
            if (countSequences >= args.filelength):
               countSequences = 0
               splitFiles.append(args.output_dir_var + inputFileName + '_%d.txt' %(countFiles))
               countFiles += 1
               outputSplit.close()
               outputSplit = open(args.output_dir_var + inputFileName + '_%d.txt' %(countFiles), 'w')
        # Else, if the line is a sequence...
        else:

            # Remove all non-characters from the sequence
            re.sub(r'[^a-zA-Z]','', line)

            # Count the number of {atcg} and {xn} in all the sequences
            characters = []
            characters = list(line)
            
            for character in characters:
                countTotal += 1
                # If input is nucleotides, count nucleotides
                if args.reftype == 'n':
                    if (re.match(r'[acgtACGT]', character)):
                        countNucleotides += 1
                    elif (re.match(r'[xnXN]', character)):
                        countXN += 1
                    else:
                        countUndef += 1
                # Else, if input is amino acids, count amino acids
                elif args.reftype == 'a':
                    if re.match(r'[abcdefghiklmnpqrstuvwyzABCDEFGHIKLMNPQRSTUVWYZ*]', character):
                        countNucleotides += 1
                    elif re.match(r'[xnXN]', character):
                        countXN += 1
                    else:
                        countUndef += 1

        # Write the lines to the appropriate files
        outputSplit.write(line)
        outputFormatted.write(line)

    # Close the files
    input.close()
    outputSplit.close()
    outputFormatted.close()

    # If there's only one input file, add it to the list of split input files
    if not splitFiles:
        splitFiles.append(args.output_dir_var + inputFileName + '_%d.txt' %(countFiles))

    # Exit the program if character count is 0
    if (countTotal == 0):
        sys.exit('ERROR: Your input file appears to be corrupted. No sequences were found!\n')

    # Exit the program if all sequences are composed only of X or N
    elif (countXN == countTotal):
        sys.exit('ERROR: Your sequence(s) contain only X or N!\n')

    # Exit the program if less than half of the characters are nucleotides
    elif (float(countNucleotides / float(countTotal)) < 0.5):
        if args.reftype == 'n':
            sys.exit('ERROR: Your sequence(s) most likely contain no DNA!\n')
        elif args.reftype == 'a':
            sys.exit('ERROR: Your sequence(s) most likely contain no AA!\n')

    return splitFiles


def runBlast(args, splitFiles):
    """Runs the BLAST algorithm on each of the split input files."""

    print ''
    print 'Run BLAST'

    # For each file containing a maximum of the specified number of sequences...
    alignment_data_dir = args.mltreemap_data + PATHDELIM + \
                         'data' + PATHDELIM + \
                         args.reference_data_prefix + 'alignment_data' + PATHDELIM + \
                         '*.fa'
    db_nt = '-db "'
    db_aa = '-db "'

    for file in glob.glob(alignment_data_dir):
        if re.match(r'.*rRNA\.fa\Z', file):
            db_nt += file + ' '
        else:
            db_aa += file + ' '

    db_nt += '"'
    db_aa += '"'

    for splitFile in sorted(splitFiles):

        # Ensure splitFile is a .txt file; save file name if so, die otherwise
        blastInputFileName = ''
        if (not re.match(r'\A.+/(.+)\.txt\Z', splitFile)):
            sys.exit('ERROR: Something is wrong with the directory of the BLAST input file!\n')
        else:
            blastInputFileName = re.match(r'\A.+/(.+)\.txt\Z', splitFile).group(1)

        # Run the appropriate BLAST command(s) based on the input sequence type
        if args.reftype == 'n':
            command = args.executables + PATHDELIM + 'blastx ' + \
                      '-query ' + splitFile + ' ' + db_aa + ' ' + \
                      '-evalue 0.01 -max_target_seqs 20000 ' + \
                      '-dbsize 1000000 -outfmt 6 '
            if args.num_threads:
                if (int(args.num_threads) >= 1) and (int(args.num_threads) < int(available_cpu_count())):
                    command += '-num_threads ' + str(int(args.num_threads)) + ' '
                else:
                    command += '-num_threads ' + str(1) + ' '
            command += '>> ' + args.output_dir_var + blastInputFileName + '.BLAST_results_raw.txt'
            os.system(command)
            command = args.executables + PATHDELIM + 'blastn ' + \
                      '-query ' + splitFile + ' ' + db_nt + ' ' + \
                      '-evalue 0.01 -max_target_seqs 20000 ' + \
                      '-dbsize 1000000 -outfmt 6 '
            if args.num_threads:
                if (int(args.num_threads) >= 1) and (int(args.num_threads) < int(available_cpu_count())):
                    command += '-num_threads ' + str(int(args.num_threads)) + ' '
                else:
                    command += '-num_threads ' + str(1) + ' '
            command += '>> ' + args.output_dir_var + blastInputFileName + '.rRNA_BLAST_results_raw.txt'
            os.system(command)
        elif args.reftype == 'a':
            command = args.executables + PATHDELIM + 'blastp ' + \
                      '-query ' + splitFile + ' ' + db_aa + ' ' + \
                      '-evalue 0.01 -max_target_seqs 20000 ' + \
                      '-dbsize 1000000 -outfmt 6 '
            if args.num_threads:
                if (int(args.num_threads) >= 1) and (int(args.num_threads) < int(available_cpu_count())):
                    command += '-num_threads ' + str(int(args.num_threads)) + ' '
                else:
                    command += '-num_threads ' + str(1) + ' '
            command += '>> ' + args.output_dir_var + blastInputFileName + '.BLAST_results_raw.txt'
            os.system(command)

        # Remove the BLAST input file
        if path.exists(splitFile):
           os.remove(splitFile)

    # TK? Remove empty BLAST result raw files; store non-empty files in a list

 
def readBlastResults(args):
    """
    Deletes empty BLAST results files.
    Returns a list of non-empty BLAST results files.
    """

    rawBlastResultFiles = []

    for file in glob.glob(args.output_dir_var + '*BLAST_results_raw.txt'):
        file.rstrip('\r\n')
        if path.getsize(file) <= 0:
            os.remove(file)
        else:
            rawBlastResultFiles.append(file)
    
    return rawBlastResultFiles


def parseBlastResults(args, rawBlastResultFiles, cog_list):
    """Returns an Autovivification of purified (eg. non-redundant) BLAST hits."""

    counter = 0
    purifiedBlastHits = Autovivify()
    for file in sorted(rawBlastResultFiles):
        try:     
            blastResults = open(file, 'r')
        except IOError:
            print "ERROR: Cannot open BLAST outputfile " + file
            continue

        contigs = Autovivify()
        identifier = 0
        for line in blastResults:
            # Clear the variables referencing the contig, COG, query start, query end, reference start, reference end, and bitscore
            # Interpret the BLAST hit, and assign the details accordingly
            tempContig, tempDetailedCOG, _, _, _, _, tempQStart, tempQEnd, tempRStart, tempREnd, _, tempBitScore = line.split('\t')
            tempREnd = int(tempREnd)
            tempRStart = int(tempRStart)
            tempQEnd = int(tempQEnd)
            tempQStart = int(tempQStart)
            tempBitScore = float(tempBitScore)

            # Skip to next BLAST hit if bit score is less than user-defined minimum
            if (tempBitScore <= args.bitscore):
                continue

            # Determine the direction of the hit relative to the reference
            direction = 'forward'
            if tempRStart > tempREnd:
                temp = tempRStart
                tempRStart = tempREnd
                tempREnd = temp
                direction = 'reverse'
            if tempQStart > tempQEnd:
                temp = tempQStart
                tempQStart = tempQEnd
                tempQEnd = temp
                if (direction == 'reverse'):
                    sys.exit('ERROR: Parsing error with the BLAST results. Please notify the authors about ' + tempContig + ' at ' +\
                              tempDetailedCOG + 'q('+tempQEnd+'..'+tempQStart+'),r('+tempREnd+'..'+tempRStart+')')
                direction = 'reverse'

            # Trim COG name to last 7 characters of detailed COG name
            # TK - This will be important to note in the user's manual, especially if we enable people to add their own COGs later
            if re.match(r'.*(.{7})\Z', tempDetailedCOG):
                tempCOG = re.match(r'.*(.{7})\Z', tempDetailedCOG).group(1)
            else:
                sys.exit('ERROR: Could not detect the COG of sequence ' + tempDetailedCOG)

            # Save contig details to the list
            contigs[tempContig][identifier]['bitscore'] = tempBitScore
            contigs[tempContig][identifier]['cog'] = tempCOG
            contigs[tempContig][identifier]['seq_start'] = tempQStart
            contigs[tempContig][identifier]['seq_end'] = tempQEnd
            contigs[tempContig][identifier]['direction'] = direction
            contigs[tempContig][identifier]['validity'] = True
            identifier += 1

        # Close the file
        blastResults.close()

        # Purify the BLAST hits

        # For each contig sorted by their stringwise comparison...
        for contig in sorted(contigs.keys()):
            identifier = 0

            # For each blast result for that contig...
            for base_blast_result_raw_identifier in sorted(contigs[contig].keys()):
                base_bitscore = contigs[contig][base_blast_result_raw_identifier]['bitscore']
                base_cog = contigs[contig][base_blast_result_raw_identifier]['cog']
                base_start = contigs[contig][base_blast_result_raw_identifier]['seq_start']
                base_end = contigs[contig][base_blast_result_raw_identifier]['seq_end']
                direction = contigs[contig][base_blast_result_raw_identifier]['direction']
                base_length = base_end - base_start

                # Skip if base_bitscore is less than user specified minimum bitscore
                if (base_bitscore < args.bitscore):
                    continue

                # Set validity to 0 if COG is not in list of MLTreeMap COGs
                if not base_cog in cog_list['all_cogs']:
                    contigs[contig][base_blast_result_raw_identifier]['validity'] = False

                # Compare the BLAST hit (base) against all others
                # There may be several opinions about how to do this. This way is based on the original MLTreeMap
                # ----A----  --C--
                #        ---B---
                # A kills B, B kills C. (Another approach would be to let C live, but the original MLTreeMap authors don't expect C to be useful)

                for check_blast_result_raw_identifier in sorted(contigs[contig]):
                    check_bitscore = contigs[contig][check_blast_result_raw_identifier]['bitscore']
                    check_cog = contigs[contig][check_blast_result_raw_identifier]['cog']
                    check_start = contigs[contig][check_blast_result_raw_identifier]['seq_start']
                    check_end = contigs[contig][check_blast_result_raw_identifier]['seq_end']
                    check_length = check_end - check_start

                    # Don't compare base hit against itself; skip to next iteration
                    if base_blast_result_raw_identifier == check_blast_result_raw_identifier:
                        continue

                    # Compare the base and check BLAST hits
                    info = Autovivify()
                    info['base']['start'] = base_start
                    info['base']['end'] = base_end
                    info['check']['start'] = check_start
                    info['check']['end'] = check_end
                    overlap = calculate_overlap(info)
                    counter +=1

                    # Check for validity for hits with overlap
                    if overlap > 0:
                        if overlap  > 0.5*base_length and base_bitscore < check_bitscore:
                            contigs[contig][base_blast_result_raw_identifier]['validity'] = False
                        elif overlap > 0.5*check_length and check_bitscore < base_bitscore:
                            contigs[contig][check_blast_result_raw_identifier]['validity'] = False
                        elif base_start == check_start and base_end == check_end:
                            # If both are the same, keep only the one with the smaller identifier
                           if check_blast_result_raw_identifier > base_blast_result_raw_identifier:
                                contigs[contig][check_blast_result_raw_identifier]['validity'] = False
                           else:
                                contigs[contig][base_blast_result_raw_identifier]['validity'] = False

                # Save purified hits for valid base hits
                if contigs[contig][base_blast_result_raw_identifier]['validity']:
                     purifiedBlastHits[contig][identifier]['bitscore'] = base_bitscore
                     purifiedBlastHits[contig][identifier]['cog'] = base_cog
                     purifiedBlastHits[contig][identifier]['start'] = base_start
                     purifiedBlastHits[contig][identifier]['end'] = base_end
                     purifiedBlastHits[contig][identifier]['direction'] = direction
                     purifiedBlastHits[contig][identifier]['is_already_placed'] = False
                     identifier += 1

    # Print the BLAST results for each contig
    for contig in sorted(purifiedBlastHits.keys()):
        outfile = args.output_dir_var + contig + '_blast_result_purified.txt'
        out = open(outfile, 'w')
        sorting_hash = {}

        # Identify the first instance of each bitscore
        for identifier in sorted(purifiedBlastHits[contig].keys()):
            if not purifiedBlastHits[contig][identifier]['bitscore'] in sorting_hash:
               sorting_hash[purifiedBlastHits[contig][identifier]['bitscore']] = {}
            sorting_hash[purifiedBlastHits[contig][identifier]['bitscore']][identifier] = 1

        # Print the (potentially reduced set of) BLAST results ordered by decreasing bitscore
        for bitscore in sorted(sorting_hash.keys(), reverse=True):

            for identifier in sorted(sorting_hash[bitscore]):
                out.write(contig + '\t' + str(purifiedBlastHits[contig][identifier]['start']) + '\t' +\
                str(purifiedBlastHits[contig][identifier]['end']) + '\t' +\
                str(purifiedBlastHits[contig][identifier]['direction']) + '\t' +\
                purifiedBlastHits[contig][identifier]['cog'] + '\t' + str( bitscore) + '\n')

        out.close()

    return purifiedBlastHits


def blastpParser(args, shortened_sequence_files, blast_hits_purified):
    """
    For each contig, produces a file similar to the Genewise output file 
     (this is in cases where Genewise is unnecessary because it is already an AA sequence.
    Returns an Autovivification of the output file for each contig.
    """

    blastpSummaryFiles = Autovivify()

    for contig in sorted(blast_hits_purified.keys()): 
        output_file = args.output_dir_var + contig + '_blast_result_summary.txt'
        try:
            output = open(output_file, 'w')
        except IOError:
            sys.exit('ERROR: Unable to open ' + output_file + '!\n')
        blastpSummaryFiles[contig][output_file] = 1
        shortened_sequence_file = args.output_dir_var + contig + '_sequence_shortened.txt'
        try:
            sequence_file = open(shortened_sequence_file, 'r')
        except IOError:
            sys.exit('ERROR: Could not open ' + shortened_sequence_file + '!\n')
        flagSeq = 0
        sequence = ''

        # Get the sequence from the shortened sequence file
        for line in sequence_file:
            if re.search('\A>', line):
                if flagSeq == 1:
                    sys.exit('ERROR: Unexpected multiple shortened sequences found!\n')
                flagSeq = 1
                continue
            else:
                line.strip()
                sequence += line

        # Write the output file to imitate the Genewise results
        for count in sorted(blast_hits_purified[contig].keys()):
            output.write(str(blast_hits_purified[contig][count]['cog']) + '\t')
            output.write(str(blast_hits_purified[contig][count]['start']) + '\t')
            output.write(str(blast_hits_purified[contig][count]['end']) + '\t')
            output.write(str(blast_hits_purified[contig][count]['direction']) + '\t')
            output.write(str(sequence) + '\n')

        output.close()

    return blastpSummaryFiles


def produceGenewiseFiles(args, blast_hits_purified):
    """
    Takes an Autovivification of purified BLAST hits and uses these 
    to produce the input files needed for Genewise.
    Returns an Autovivification mapping the contig to its sequence's 
    start and end positions for Genewise.
    Returns a list of files to be run through Genewise.
    """

    flanking_length = 1000 # Recommended: 1000
    prae_contig_coordinates = Autovivify()
    contig_coordinates = Autovivify()
    shortened_sequence_files = {}

    for contig in sorted(blast_hits_purified.keys()):
        for base_identifier in sorted(blast_hits_purified[contig].keys()):
            # Skip rRNA hits for now (we work with them later)
            if re.search("rRNA", blast_hits_purified[contig][base_identifier]['cog']):
                continue

            # Skip hits which have already been placed; otherwise, mark them as placed
            if blast_hits_purified[contig][base_identifier]['is_already_placed']:
                continue

            blast_hits_purified[contig][base_identifier]['is_already_placed'] = True
            base_start = blast_hits_purified[contig][base_identifier]['start'] - flanking_length
            base_end = blast_hits_purified[contig][base_identifier]['end'] + flanking_length
            nr_of_blast_hits = len(blast_hits_purified[contig].keys())
            check_identifier =0
            while check_identifier < nr_of_blast_hits:
                # Skip rRNA hits for now (we work with them later)
                if re.search(r'rRNA', blast_hits_purified[contig][check_identifier]['cog']):
                    check_identifier +=1
                    continue

                # Skip hits which have already been placed; otherwise, mark them as placed
                if blast_hits_purified[contig][check_identifier]['is_already_placed']:
                    check_identifier +=1
                    continue

                check_start = blast_hits_purified[contig][check_identifier]['start'] - flanking_length
                check_end = blast_hits_purified[contig][check_identifier]['end'] + flanking_length

                # Check for overlap
                if base_start <= check_start and check_start <= base_end and base_end <= check_end:
                    # Base  --------
                    # Check     --------
                    base_end = check_end
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                elif base_start <= check_start and check_end <= base_end:
                    # Base  --------
                    # Check   ----
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                elif check_start <= base_start and base_start <= check_end and check_end <= base_end:
                    # Base      --------
                    # Check --------
                    base_start = check_start
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                elif check_start <= base_start and base_end <= check_end:
                    # Base    ----
                    # Check --------
                    base_start = check_start
                    base_end = check_end
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                check_identifier += 1

            prae_contig_coordinates[contig][base_start][base_end] = 1

    # Produce the input files for Genewise
    input = open(args.formatted_input_file, 'r')
    contig_name = ''
    sequence = ''
    line = 'x'

    while line:
        line= input.readline()
        line =  line.strip()
        line = re.sub(r'\s', '_', line)
        searchmatch =re.search(r'\A>(.+)', line)

        if searchmatch or not line:
            if not line:
               sequence += line
            if contig_name in prae_contig_coordinates:
                sequence_length = len(sequence)
                shortened_sequence="" 

                # Start searching for the information to shorten the file.
                for start_B in sorted(prae_contig_coordinates[contig_name].keys()) :
                    for end_B in sorted(prae_contig_coordinates[contig_name][start_B].keys()) :

                        # Ok, now we have all information about the hit. Correct start and end if needed: 
                        if start_B < 0:
                           start_B = 0 

                        if end_B >= sequence_length:
                           end_B = sequence_length -1 
      
                        # Note: Genewise (GW) positions start with 1, Blast (B) positions with 0 -> thus differenciate between start_B and start_GW
                        shortened_start_GW = len(shortened_sequence) + 1 
                        count = -1
                        for nucleotide in sequence: 
                            count += 1     
                            if not (count >= start_B and count <= end_B):
                               continue
                            shortened_sequence += nucleotide

                        shortened_end_GW = len(shortened_sequence)
                        addition_factor = (start_B + 1) - shortened_start_GW #$start_B + 1 == $start_GW
                        contig_coordinates[contig_name][shortened_start_GW][shortened_end_GW] = addition_factor

                try:
                    with open(args.output_dir_var + contig_name + "_sequence.txt", 'w') as f:
                       fprintf(f, "%s\n", ">"+ contig_name + "\n" + sequence)
                    f.close()
                except:
                    print  "ERROR: Can't create " + args.output_dir_var + contig_name + "_sequence.txt!" 

                try:   
                   with open(args.output_dir_var + contig_name + "_sequence_shortened.txt", 'w') as f:
                      fprintf(f, "%s\n",">" + contig_name + "\n" + shortened_sequence)
                   f.close()
                   shortened_sequence_files[args.output_dir_var +  contig_name + "_sequence_shortened.txt"]=contig_name
                except:
                   print "ERROR: Can't create " + args.output_dir_var +  contig_name +"_sequence_shortened.txt!" 

            if searchmatch:
               contig_name = searchmatch.group(1)
               sequence = ""

        else:
            sequence += line

    input.close()
    return contig_coordinates, shortened_sequence_files


def fprintf(file, fmt, *args): 
   """A helper function used to print to a specified file."""

   file.write(fmt % args)


def startGenewise(args, shortened_sequence_files, blast_hits_purified):
    """
    Runs Genewise on the provided list of sequence files.
    (The provided Autovivification of purified BLAST hits is used for file naming purposes).
    Returns an Autovivification mapping the Genewise output files to each contig.
    """

    print 'Run Genewise'
    genewise_outputfiles = Autovivify()

    # For each file which has been shortened by produceGenewiseFiles...

    for shortened_sequence_file in sorted(shortened_sequence_files.keys()):
        contig = shortened_sequence_files[shortened_sequence_file]

        # For each identifier associated with this contig in the output of parseBlastResults
        for identifier in sorted(blast_hits_purified[contig].keys()):
            cog = blast_hits_purified[contig][identifier]['cog']

            # Prepare the output file name, and store it
            genewise_outputfile = args.output_dir_var + contig + '_' + cog + '_genewise.txt'
            genewise_outputfiles[contig][genewise_outputfile] = 1

            # Prepare the Genewise command and run it
            mltreemap_dir = args.mltreemap_data + PATHDELIM + 'data' + PATHDELIM
            genewise_support = mltreemap_dir + PATHDELIM + 'genewise_support_files' + PATHDELIM
            hmm_dir = mltreemap_dir + "hmm_data" + PATHDELIM
            genewise_command = args.executables + PATHDELIM + 'genewise ' + \
                               hmm_dir + cog + '.hmm ' + \
                               shortened_sequence_file + ' -init local -quiet -gene ' + \
                               genewise_support + 'human.gf -matrix ' + \
                               genewise_support + 'blosum62.bla -codon ' + \
                               genewise_support + 'codon.table -hmmer -subs' + \
                               ' 0.01 -indel 0.01 -gap 11 -ext 1 -both -pep -sum > ' + genewise_outputfile
            os.system(genewise_command)

    # Return the list of output files for each contig
    return genewise_outputfiles

def parse_genewise_results(args, genewise_outputfiles, contig_coordinates):
    """
    Uses the provided Autovivification of Genewise output files and the provided
    Autovivification mapping the contig to its Genewise sequence's start and end
    points to produce files summarizing the purified Genewise results.
    Returns an Autovivification mapping the summary files to each contig.
    """

    genewise_summary_files = Autovivify()

    # For each contig analyzed by Genewise...
    for contig in sorted(genewise_outputfiles.keys()):
        genewise_results_raw = Autovivify()
        genewise_results = Autovivify()
        at_least_one_hit = 0
        count = 0

        # Parse each output file of that contig
        for genewise_outputfile in sorted(genewise_outputfiles[contig].keys()):
            try:     
                input = open(genewise_outputfile, 'r')
            except IOError:
                print "ERROR: Cannot open Genewise outputfile " + genewise_outputfile
                continue

            header_count = 0
            sequence_count = -1

            for line in input:
                line.strip()

                # If the line starts with a digit, parse it
                if re.match(r'\A\d', line):

                    # Split the results based on one or more spaces between the desired data
                    bitscore, query, _, _, _, start, end, _, _ = re.split(' +', line)
                    bitscore = float(bitscore)
                    start = int(start)
                    end = int(end)

                    # If there is at least one query, take note for future use
                    if query is not None:
                        at_least_one_hit = 1

                    # Determine the direction of the predicted amino acid sequence
                    direction = 'forward'
                    if start > end:
                        temp = start
                        start = end
                        end = start
                        direction = 'reverse'

                    # Correct the positions
                    # Genewise is run on a shortened sequence, so the true positions must be calculated
                    for coords_start in sorted(contig_coordinates[contig].keys()):
                        if start >= coords_start:
                            for coords_end in sorted(contig_coordinates[contig][coords_start].keys()):
                                if end <= coords_end:
                                    addition_factor = contig_coordinates[contig][coords_start][coords_end]
                                    start += addition_factor
                                    end += addition_factor
                                    break

                    genewise_results_raw[contig][genewise_outputfile][header_count]['start'] = start
                    genewise_results_raw[contig][genewise_outputfile][header_count]['end'] = end
                    genewise_results_raw[contig][genewise_outputfile][header_count]['cog'] = query
                    genewise_results_raw[contig][genewise_outputfile][header_count]['bitscore'] = bitscore
                    genewise_results_raw[contig][genewise_outputfile][header_count]['direction'] = direction
                    header_count += 1

                # Otherwise, if the line starts with a '>', prepare to intake the sequence
                elif re.match(r'\A>', line):
                    sequence_count += 1
                    genewise_results_raw[contig][genewise_outputfile][sequence_count]['sequence'] = ''

                # If the line begins with any word character, and contains neither 'Bits' (a title line)
                # nor 'Making' (a Genewise comment about the treatment of introns)
                elif re.match(r'\A\w', line) and not re.match(r'\ABits', line) and not re.match(r'\AMaking', line):
                    genewise_results_raw[contig][genewise_outputfile][sequence_count]['sequence'] += line.strip()

            input.close()

        # Skip to next contig if there isn't at least 1 hit
        if at_least_one_hit != 1:
            continue

        # Purify the parsed results
        # For each genewise_outputfile for the contig...
        for base_genewise_outputfile in sorted(genewise_results_raw[contig].keys()):

            # For each count of the genewise_outputfile...
            for base_count in sorted(genewise_results_raw[contig][base_genewise_outputfile].keys()):
                base_start = genewise_results_raw[contig][base_genewise_outputfile][base_count]['start']
                base_end = genewise_results_raw[contig][base_genewise_outputfile][base_count]['end']
                base_cog = genewise_results_raw[contig][base_genewise_outputfile][base_count]['cog']
                base_bitscore = genewise_results_raw[contig][base_genewise_outputfile][base_count]['bitscore']
                base_direction = genewise_results_raw[contig][base_genewise_outputfile][base_count]['direction']
                base_sequence = genewise_results_raw[contig][base_genewise_outputfile][base_count]['sequence']

                # Ensure that the base_cog, base_start, and base_end are defined
                if base_cog is None or base_start is None or base_end is None:
                    error_string = 'ERROR: The file "' + base_genewise_outputfile + '" cannot be parsed!\n' +\
                                   'Please contact the authors about it. As a quick solution to the problem, ' +\
                                   'try to remove the sequence which produced this hit from your input file.\n'
                    sys.exit(error_string)

                base_length = base_end - base_start
                is_valid = 1

                # Check against all other genewise_outputfiles for that contig
                # For each genewise_outputfile for the contig...
                for check_genewise_outputfile in sorted(genewise_results_raw[contig].keys()):

                    # For each count of the genewise_outputfile...
                    for check_count in sorted(genewise_results_raw[contig][check_genewise_outputfile].keys()):

                        # Skip to next iteration if comparing the base to itself
                        if base_count == check_count:
                            continue
                        check_start = genewise_results_raw[contig][check_genewise_outputfile][check_count]['start']
                        check_end = genewise_results_raw[contig][check_genewise_outputfile][check_count]['end']
                        check_cog = genewise_results_raw[contig][check_genewise_outputfile][check_count]['cog']

                        # Ensure that the check_cog, check_start, and check_end are defined
                        if check_cog is None or check_start is None or check_end is None:
                            error_string = 'ERROR: The file "' + check_genewise_outputfile + '" cannot be parsed!\n' +\
                                           'Please contact the authors about it. As a quick solution to the problem, ' +\
                                           'try to remove the sequence which produced this hit from your input file.\n'
                            sys.exit(error_string)

                        check_length = check_end - check_start
                        info = Autovivify()
                        info['base']['start'] = base_start
                        info['base']['end'] = base_end
                        info['check']['start'] = check_start
                        info['check']['end'] = check_end
                        overlap = calculate_overlap(info)

                        # Purify the results
                        # If the size of the overlap is more than half the size of the hit...
                        if float(overlap) / float(base_length) > 0.5:

                            # And if the hit and check are the same COG...
                            if base_cog == check_cog:

                                # Keep only the longer hit, since the major difference between the hits is the length 
                                if base_length < check_length:
                                    is_valid = 0

                            # The COGs are different, so only skip the hit if it is less than half the length of the check
                            elif base_length < check_length / 2:
                                is_valid = 0

                        # But if the overlap is not more than half the size of the hit, and the hit remains valid...
                        if is_valid and base_cog == check_cog:

                            # Invalidate the hit if it is a side hit of the same COG
                            if base_length < check_length * 0.7:
                                is_valid = 0

                # If the hit is valid, save it
                if is_valid == 1:
                    genewise_results[contig][count]['start'] = base_start
                    genewise_results[contig][count]['end'] = base_end
                    genewise_results[contig][count]['cog'] = base_cog
                    genewise_results[contig][count]['direction'] = base_direction
                    genewise_results[contig][count]['sequence'] = base_sequence
                    count += 1

        # Skip to next hit if there are no valid hits
        if count <= 0:
            continue

        # Write the summary file
        genewise_summary_file = args.output_dir_var +  contig + '_genewise_result_summary.txt'
        try: 
            output = open(genewise_summary_file, 'w')
        except IOError:
            print 'ERROR: Cannot open Genewise summary file ' + genewise_summary_file + ' for writing'
            sys.exit(0)

        genewise_summary_files[contig][genewise_summary_file] = 1
        for count in sorted(genewise_results[contig].keys()):
            output.write(genewise_results[contig][count]['cog'] + '\t' +\
                         str(genewise_results[contig][count]['start']) + '\t' +\
                         str(genewise_results[contig][count]['end']) + '\t' +\
                         genewise_results[contig][count]['direction'] + '\t' +\
                         genewise_results[contig][count]['sequence'] + '\n')

        output.close()

    return genewise_summary_files

def get_rRNA_hit_sequences(args, blast_hits_purified, cog_list, genewise_summary_files):
    """
    rRNA does not get translated into protein. Regardless, we want to take the
    rRNA and summarize it in a way that is parallel to the Genewise summary files.
    This function does that using the provided Autovivification of purified BLAST
    hits, list of COGs, and Autovivification of Genewise summary files.
    Returns an Autovivification summarizing the coordinates for each rRNA hit.
    Returns a list of the rRNA summary files.
    """

# TK: ...the list of rRNA hit files is empty.

    contig_rRNA_coordinates = Autovivify()
    rRNA_hit_files = {}
    
    for contig in sorted(blast_hits_purified.keys()) :
        # note: We skipped the Genewise step (we are dealing with rRNA) but we bring the rRNA files in the
        # same structure as the Genewise summary files and bring them back into the ordinary pipeline.
        for identifier in sorted(blast_hits_purified[contig].keys()):
            if not re.search("rRNA", blast_hits_purified[contig][identifier]['cog']):
                continue

            start = blast_hits_purified[contig][identifier]["start"]
            end = blast_hits_purified[contig][identifier]["end"]
            cog = blast_hits_purified[contig][identifier]["cog"]
            direction = blast_hits_purified[contig][identifier]["direction"]
            contig_rRNA_coordinates[contig][identifier]["start"] = start
            contig_rRNA_coordinates[contig][identifier]["end"] = end
            contig_rRNA_coordinates[contig][identifier]["cog"] = cog
            contig_rRNA_coordinates[contig][identifier]["direction"] = direction
            outfile_name = args.output_dir_var +  contig + '_rRNA_result_summary.txt'   
            contig_rRNA_coordinates[contig][identifier]["outfile"] = outfile_name
            genewise_summary_files[contig][outfile_name] = 1     
            try:
               outfile = open(outfile_name, 'w')
               outfile.close() 
            except IOError, e:
               print "ERROR: Can't create " + outfile_name + '!\n' 
               sys.exit(0)

    try:
       input = open(args.input, 'r')
    except IOError, e:
       print "ERROR: Can't create " + args.input +'!\n' 
       sys.exit(0)
    contig_name = ''
    sequence = ''

    line = 'x'
    while line:
        line= input.readline()
        line =  line.strip()
        line = re.sub(r'\s', '_', line)
        searchmatch =re.search(r'\A>(.+)', line)

        if searchmatch or not line:
            if not line:
               sequence += line

            if contig_name in contig_rRNA_coordinates:
                sequence_length = len(sequence)
                shortened_sequence="" 
                #start searching for the information to shorten the file.
                for identifier in sorted(contig_rRNA_coordinates[contig_name].keys()) :
                    start = contig_rRNA_coordinates[contig_name][identifier]["start"]
                    end = contig_rRNA_coordinates[contig_name][identifier]["end"]
                    cog = contig_rRNA_coordinates[contig_name][identifier]["cog"]
                    direction = contig_rRNA_coordinates[contig_name][identifier]["direction"]
                    outfile = contig_rRNA_coordinates[contig_name][identifier]['outfile']
                    denominator = cog_list['all_cogs'][cog]
                    count = -1
                    shortened_sequence = ""
                    for nucleotide in sequence: 
                        count+=1
                        if not (count >= start and count <= end):
                           continue
                        shortened_sequence += nucleotide

                    if direction == 'reverse':
                        #ok, our hit has been on the opposite strand of the reference.
                        #to get a proper alignment, we thus have to produce a negative strand version of the input
                        nucleotides2 = ''.join(reversed(shortened_sequence))
                        shortened_sequence = ""
                        nucleotides2 = nucleotides2.lower()
                        for nucleotide in  nucleotides2:
                            if nucleotide == 't':
                                nucleotide = 'a'
                            elif nucleotide == 'a' :
                                nucleotide = 't'
                            elif nucleotide == 'c':
                                nucleotide = 'g'
                            elif nucleotide == 'g':
                                nucleotide = 'c'

                            shortened_sequence += nucleotide
                    try:
                       out = open(outfile, 'a')
                       fprintf(out,'%s\t%s\t%s\t%s\t%s\n', cog, start, end,'n/a', shortened_sequence)
                       out.close()
                    except IOError, e:
                       print "ERROR: Can't create " + outfile + '!\n' 
                       sys.exit(0)

                try:
                   output_file = open(args.output_dir_var + contig_name + '_sequence.txt', 'w')
                   fprintf(output_file, '>%s\n%s',contig_name, sequence)
                   output_file.close()
                except IOError, e:
                    print "ERROR: Can't create " + args.output_dir_var + contig_name + '_sequence.txt!\n'
                    sys.exit(0)

            if searchmatch:
               contig_name = searchmatch.group(1)
               sequence = ""
        else:
           sequence += line
    input.close()
    return contig_rRNA_coordinates, rRNA_hit_files


def prepare_and_run_hmmalign(args, genewise_summary_files, cog_list):
    """
    Runs hmmalign using the provided COG list and summary of Genewise files.
    Returns an Autovivification of the resulting files from hmmalign.
    """

    reference_data_prefix = args.reference_data_prefix
    hmmalign_singlehit_files = Autovivify(); 
    print 'Run hmmalign'

    # Run hmmalign on each Genewise summary file
    for contig in sorted(genewise_summary_files.keys()):

        for genewise_summary_file in sorted(genewise_summary_files[contig].keys()):
            try:
               input = open(genewise_summary_file, 'r')
            except IOError:  
               print  "ERROR: Can't open " + genewise_summary_file + "!\n"
               sys.exit(0)

            line = input.readline()
            line = line.strip()

            while line:
                cog, start, end, _, sequence = line.split('\t')
                denominator = cog_list["all_cogs"][cog]
                f_contig = denominator + "_" + contig
                genewise_singlehit_file = args.output_dir_var +PATHDELIM +f_contig+'_'+cog+"_"+str(start)+"_"+str(end)
                hmmalign_singlehit_files[f_contig][genewise_singlehit_file + ".mfa"] = True 
                genewise_singlehit_file_fa = genewise_singlehit_file + ".fa" 
                try:
                   outfile = open(genewise_singlehit_file_fa, 'w')
                   fprintf(outfile, '>query\n%s', sequence)
                   outfile.close()
                except IOError:
                   print  'Can\'t create ' + genewise_singlehit_file_fa +'\n'
                   sys.exit(0)                
                mltreemap_resources = args.mltreemap_data + PATHDELIM + 'data' + PATHDELIM
                hmmalign_command = [ args.executables + PATHDELIM + 'hmmalign', '-m', '--mapali',\
                                     mltreemap_resources + reference_data_prefix + 'alignment_data' + PATHDELIM +  cog + '.fa',\
                                     '--outformat', 'Clustal',\
                                     mltreemap_resources + reference_data_prefix + 'hmm_data' + PATHDELIM + cog + '.hmm',\
                                     genewise_singlehit_file_fa, '>', genewise_singlehit_file + '.mfa' ] 
                os.system(' '.join(hmmalign_command))
                line = input.readline()
                line = line.strip()

            input.close()

    return hmmalign_singlehit_files
                   

def get_non_wag_cogs(args):
    """
    Returns an Autovivification listing the COGs which don't follow 
    the WAG evolutionary model.
    """

    non_wag_cog_list = Autovivify()
    try:
        non_wag_cogs_file = args.mltreemap_data + PATHDELIM + \
                            'data' + PATHDELIM + 'tree_data' + PATHDELIM +'non_wag_cogs.txt'
        cogin = open(non_wag_cogs_file, 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + non_wag_cogs_file + '!\n')

    for line in cogin:
        line = line.strip()
        if re.search(r'\A#(.+)', line):
            denominator = re.search(r'\A#(.+)', line).group(1)
        else:
            cog, model = line.split('\t')
            non_wag_cog_list[denominator][cog] = model

    cogin.close()
    return non_wag_cog_list


def concatenate_hmmalign_singlehits_files(args, hmmalign_singlehit_files, non_wag_cog_list):
    """
    Concatenates the hmmalign files using the provided 
    Autovivifications of hmmalign files and non-WAG COGs.
    Returns a list of the files containing the concatenated hmmalign results.
    Returns a list of the model used for each file.
    Returns a list of the number of sequences found in each file.
    """

    # For each type of gene...
    concatenated_mfa_files = {}
    models_to_be_used = {}
    nrs_of_sequences = {}

    for f_contig in sorted(hmmalign_singlehit_files.keys()):
        # Determine what type of gene is currently represented, or die an error
        sequences = Autovivify()
        model_to_be_used = ""
        query_sequence = ""
        denominator = ""
        if re.search(r'\A(.)', f_contig):
            denominator = re.match(r'\A(.)', f_contig).group(1)
        else:
            sys.exit('ERROR: The analysis type could not be parsed from ' + f_contig + '!\n')

        # For each file...
        for hmmalign_singlehit_file in sorted(hmmalign_singlehit_files[f_contig].keys()):
            # Open the file
            try:
                input = open(hmmalign_singlehit_file, 'r')
            except IOError:
                sys.exit('Can\'t open ' + hmmalign_singlehit_file + '!\n')
            reached_data_part = False
            # Determine the best AA model
            if re.search(r'\A.+_(.{7})_\d+_\d+\.mfa\Z', hmmalign_singlehit_file):
                cog = re.search(r'\A.+_(.{7})_\d+_\d+\.mfa\Z', hmmalign_singlehit_file).group(1)
            else:
                sys.exit('ERROR: The COG could not be parsed from ' + hmmalign_singlehit_file + '!\n')
            if non_wag_cog_list[denominator][cog] and model_to_be_used != 'PROTGAMMAWAG':
                model_to_be_used = non_wag_cog_list[denominator][cog]
            else:
                model_to_be_used = 'PROTGAMMAWAG'

            # Get sequence from file
            for _line in input:
                line = _line.strip()
                if re.search(r'query', line):
                    reached_data_part = True
                if not reached_data_part:
                    continue
                searchResult =  re.search(r'\A(.+) (\S+)\Z', line)
                if searchResult:
                    name_long = searchResult.group(1)
                    sequence_part = searchResult.group(2)
                    sequence_name = ''
                    if re.search(r'query', name_long):
                        query_sequence += sequence_part
                    elif re.search(r'(\d+)_', name_long):
                        sequence_name = re.search(r'(\d+)_', name_long).group(1)
                        if sequences[sequence_name]:
                            sequences[sequence_name] += sequence_part
                        else:
                            sequences[sequence_name] = sequence_part

            input.close()

        models_to_be_used[f_contig] = model_to_be_used
        concatenated_mfa_files[f_contig] = args.output_dir_var +  f_contig + '.mfa'
        # Write to the output file
        try:
            output = open(args.output_dir_var + f_contig + '.mfa', 'w')
        except IOError:
            sys.exit('ERROR: Can\'t create ' + args.output_dir_var +  f_contig + '.mfa\n')
        output.write('>query\n' + query_sequence + '\n')
        nrs_of_sequences[f_contig] = 1

        for sequence_name in sorted(sequences.keys()):
            nrs_of_sequences[f_contig] += 1
            sequence = sequences[sequence_name]
            output.write('>' + sequence_name + '\n' + sequence + '\n')
        output.close()

    return (concatenated_mfa_files, nrs_of_sequences, models_to_be_used)


def start_gblocks(args, concatenated_mfa_files, nrs_of_sequences):
    """
    Runs Gblocks using the provided lists of the concatenated 
    hmmalign files, and the number of sequences in each file.
    Returns a list of files resulting from Gblocks.
    """

    gblocks_files = {}
    print 'Run Gblocks'
    
    for f_contig  in sorted(concatenated_mfa_files.keys()) :
        concatenated_mfa_file = concatenated_mfa_files[f_contig]
        nr_of_sequences = nrs_of_sequences[f_contig]
        min_flank_pos = int(nr_of_sequences * 0.55)
        gblocks_file = concatenated_mfa_file+ "-gb"
        gblocks_files[f_contig] = gblocks_file;
        gblocks_command = [  args.executables + PATHDELIM + "Gblocks" ]
        gblocks_command.append(concatenated_mfa_file)
        gblocks_command += ['-t=p', '-s=y', '-u=n', '-p=t', '-b3=15',\
                                 '-b4=3', '-b5=h', '-b2='+str(min_flank_pos),\
                                 '>', '/dev/null']
        os.system(' '.join(gblocks_command))
    
    return gblocks_files


def produce_phy_file(args, gblocks_files, nrs_of_sequences):
    """
    Produces phy files from the provided list of Gblocks result files, 
    and the number of sequences in each file.
    Returns an Autovivification containing the names of the produced phy files.
    """

    phy_files = Autovivify()
    sequence_lengths = Autovivify()

    # Open each Gblocks result file
    for f_contig in sorted(gblocks_files.keys()):
        sequences_for_phy = Autovivify()
        do_not_continue = 0
        sequences_raw = Autovivify()
        gblocks_file = gblocks_files[f_contig]
        try:
            input = open(gblocks_file, 'r')
        except IOError:
            sys.exit('ERROR: Can\'t open ' + gblocks_file + '!\n')

        for line in input:
            line = line.strip()
            seq_name_search =  re.search(r'\A>(.+)', line)
            if seq_name_search:
                seq_name = seq_name_search.group(1)
                # Flag the user if the reference alignment contains the number -666, which is needed later in the code
                if seq_name == -666:
                    sys.exit('ERROR: Your reference alignment contains element with the number -666. ' +\
                             'Please change it, because this number is needed for internal purposes.\n')
                if seq_name == 'query':
                   seq_name = -666
            else:
                line = re.sub(r' ', '', line)
                if seq_name in sequences_raw:
                    sequences_raw[seq_name] += line
                else:
                    sequences_raw[seq_name] = line

        input.close()

        # Ensure the sequences contain only valid characters for RAxML
        for seq_name in sorted(sequences_raw.keys()):
            if do_not_continue == 1:
                continue
            sequence = sequences_raw[seq_name]
            count = 0
            sequence_lengths[f_contig] = len(sequence)
            sequence = re.sub(r'\.', 'X', sequence)
            sequence = re.sub(r'\*', 'X', sequence)
            sequence = re.sub('-', 'X', sequence)
            if re.search(r'\AX+\Z', sequence):
                sequence = re.sub('X', 'V', sequence, 1)
            if seq_name == -666:
                seq_dummy = re.sub('X', '', sequence)
                if len(seq_dummy) < args.gblocks:
                    do_not_continue = 1
                    exit_file_name = args.output_dir_var + f_contig + '_exit_after_Gblocks.txt'
                    try:
                        output = open(exit_file_name, 'w')
                    except IOError:
                        sys.exit('ERROR: Can\'t open ' + exit_file_name + '!\n')
                    output.write('final alignment after gblocks is too short (<' + str(args.gblocks) + 'AAs) ' +\
                                 '-  insufficient number of marker gene residues in query sequence.\n')
                    output.close()
                    continue
            sub_sequences = re.findall(r'.{1,50}', sequence)

            for sub_sequence in sub_sequences:
                sub_sequence = re.sub('U', 'T', sub_sequence) # TK: This for debug; got error from RAxML when encountering Uracil
                sequences_for_phy[f_contig][count][int(seq_name)] = sub_sequence
                count += 1

        if do_not_continue == 1:
            continue

        # Write the sequences to the phy file
        phy_file_name = args.output_dir_var + f_contig + '.phy'
        phy_files[f_contig] = phy_file_name
        try:
            output = open(phy_file_name, 'w')
        except IOError:
            sys.exit('ERROR: Can\'t open ' + phy_file_name + '!\n')
        nr_of_sequences = nrs_of_sequences[f_contig]
        output.write(' ' + str(nr_of_sequences) + '  ' + str(sequence_lengths[f_contig]) + '\n')

        for count in sorted(sequences_for_phy[f_contig].keys()):
            for seq_name in sorted(sequences_for_phy[f_contig][count].keys()):
                sequence_part = sequences_for_phy[f_contig][count][seq_name]
                if count== 0:
                    print_seqname = seq_name
                    if seq_name == -666:
                        print_seqname = 'query'
                    output.write(str(print_seqname))
                    length = len(str(print_seqname))
                    c = length
                    while c < 10:
                        output.write(' ')
                        c += 1
                output.write(sequence_part + '\n')

            output.write('\n')
        output.close()

    return phy_files


#TK undef_hashes on genewise_summary_files, hmmalign_singlehit_files, concatenated_mfa_files, nrs_of_sequences, gblocks_files


def start_RAxML(args, phy_files, cog_list, models_to_be_used):
    """
    Run RAxML using the provided Autovivifications of phy files and COGs, 
    as well as the list of models used for each COG.
    Returns an Autovivification listing the output files of RAxML.
    Returns an Autovivification containing the reference tree 
    file associated with each functional or rRNA COG.
    """

    expected_raxml_outfiles = Autovivify()
    raxml_outfiles = Autovivify()
    print 'Run RAxML'
    # Notify the user, if they've indicated otherwise, that bootstraps cannot be used with the Maximum Parsimony settings of RAxML.
    # TK: Should this be moved to the beginning of the program when the user first specifies their options?
    if args.bootstraps > 1 and args.phylogeny == 'p':
        print 'ATTENTION: You intended to do ' + str(args.bootstraps) + ' bootstrap replications. Unfortunately, bootstrapping is ' +\
              'disabled in the parsimony mode of MLTreeMap. The pipeline will continue without bootstrapping.\n'
        args.bootstraps = 1
    bootstrap_replicates = args.bootstraps
    args2 = Autovivify()
    mltree_resources = args.mltreemap_data + PATHDELIM + 'data' + PATHDELIM

    for f_contig in sorted(phy_files.keys()):
        # Establish the reference tree file to be used for this contig
        reference_tree_file = mltree_resources + 'tree_data' + PATHDELIM + args.reference_tree
        phy_file = phy_files[f_contig]
        if re.search(r'\A(.)', f_contig):
            denominator = re.search(r'\A(.)', f_contig).group(0)
        if not denominator == 'p' and not denominator == 'g' and not denominator == 'i':

            for cog in sorted(cog_list['all_cogs'].keys()):
                if not cog_list['all_cogs'][cog] == denominator:
                    continue
                reference_tree_file = mltree_resources + 'tree_data' + PATHDELIM + cog + '_tree.txt'
                break

        # Determine the output file names, and remove any pre-existing output files
        args2['reference_tree_file_of_denominator'][denominator] = reference_tree_file
        raxml_files = [args.output_dir_var + 'RAxML_info.' + f_contig,\
                       args.output_dir_var + 'RAxML_labelledTree.' + f_contig,\
                       args.output_dir_var + 'RAxML_classification.' + f_contig]

        for raxml_file in raxml_files:
            try:
                shutil.rmtree(raxml_file) 
            except OSError:
                pass

        raxml_option = args.phylogeny
        model_to_be_used = models_to_be_used[f_contig]
        if model_to_be_used is None:
            sys.exit('ERROR: No best AA model could be detected for the ML step!\n')
        # Set up the command to run RAxML
        raxml_command = [ args.executables + PATHDELIM + 'raxmlHPC', '-m', model_to_be_used]
        if bootstrap_replicates > 1:
            raxml_command += [ '-x', '12345', '-#', bootstrap_replicates]
        # Run RAxML using multiple threads, if CPUs available
        if args.num_threads:
            if ( int(args.num_threads) >= 1 ) and ( int(args.num_threads) <= available_cpu_count() ): 
                raxml_command += ['-T', str(int(args.num_threads))]
        raxml_command += [ '-s', phy_file, '-t', reference_tree_file, '-f', str(raxml_option), '-n', f_contig,\
                           '-w', str(args.output_dir_var), '>',\
                           str(args.output_dir_var) + str(f_contig) + '_RAxML.txt']
        os.system(' '.join(raxml_command))

    # Rename the RAxML output files
    for f_contig in sorted(phy_files.keys()):
        denominator = ''
        if re.match(r'\A(.)', f_contig):
            denominator = re.match(r'\A(.)', f_contig).group(1)
        move_command = [ 'mv', str(args.output_dir_var) + 'RAxML_info.' + str(f_contig), \
                         str(args.output_dir_var) + str(f_contig) + '.RAxML_info.txt']
        if os.path.exists(str(args.output_dir_var) + 'RAxML_info.' + str(f_contig)):
            os.system(' '.join(move_command))
        if raxml_option == 'v':
            raxml_outfiles[denominator][f_contig]['classification'] = str(args.output_dir_var) + str(f_contig) + '.RAxML_classification.txt'
            raxml_outfiles[denominator][f_contig]['labelled_tree'] = str(args.output_dir_var) + str(f_contig) + '.originalRAxML_labelledTree.txt'
            move_command1 = [ 'mv', str(args.output_dir_var) + 'RAxML_classification.' + str(f_contig),\
                              str(raxml_outfiles[denominator][f_contig]['classification'])]
            move_command2 = [ 'mv', str(args.output_dir_var) + 'RAxML_originalLabelledTree.' + str(f_contig),\
                              str(raxml_outfiles[denominator][f_contig]['labelled_tree'])]
            remove_command = [ 'rm', str(args.output_dir_var) + 'RAxML_labelledTree.' + str(f_contig)]
            if os.path.exists(str(args.output_dir_var) + 'RAxML_classification.' + str(f_contig)):
                os.system(' '.join(move_command1))
            if os.path.exists(str(args.output_dir_var) + 'RAxML_originalLabelledTree.' + str(f_contig)):
                os.system(' '.join(move_command2))
            if os.path.exists(str(args.output_dir_var) + 'RAxML_labelledTree.' + str(f_contig)):
                os.system(' '.join(remove_command))
            else:
                print "Some files were not successfully created for", str(f_contig) #CML
        elif raxml_option == 'p':
            raxml_outfiles[denominator][f_contig] = str(args.output_dir_var) + str(f_contig) + '.RAxML_parsimonyTree.txt'
            move_command1 = [ 'mv', str(args.output_dir_var) + 'RAxML_parsimonyTree.' + str(f_contig),\
                              str(raxml_outfiles[denominator][f_contig])]
            os.system(' '.join(move_command1))
        else:
            sys.exit('ERROR: The chosen RAxML mode is invalid. This should have been noticed earlier by MLTreeMap.' +\
                     'Please notify the authors\n')

    return raxml_outfiles, args2


def parse_RAxML_output(args, args2, tree_numbers_translation, raxml_outfiles, text_of_analysis_type):
    """
    Parse the RAxML output files.
    Returns an Autovivification of the final RAxML output files.
    """

    raxml_option = args.phylogeny
    print 'Finishing'
    output_directory_final_RAxML = args.output_dir_raxml
    final_RAxML_output_files = Autovivify()

    for denominator in sorted(raxml_outfiles.keys()):
        description_text = '# ' + str(text_of_analysis_type[denominator]) + '\n'
        reference_tree_file = args2['reference_tree_file_of_denominator'][denominator]
        terminal_children_strings_of_reference = read_and_understand_the_reference_tree(reference_tree_file)
        content_of_previous_labelled_tree_file = ''
        rooted_labelled_trees = ''
        insertion_point_node_hash = ''
        final_assignment_target_strings = Autovivify()

        for f_contig in sorted(raxml_outfiles[denominator].keys()):
            denominator = ''
            if re.search(r'\A(.)', f_contig):
                denominator = re.search(r'\A(.)', f_contig).group(1)
            content_of_labelled_tree_file = ''
            assignments = Autovivify()
            nr_of_assignments = 0

            if raxml_option == 'v':
                classification_file = raxml_outfiles[denominator][f_contig]['classification']
                labelled_tree_file = raxml_outfiles[denominator][f_contig]['labelled_tree']
                try:
                    input = open(labelled_tree_file, 'r')
                except IOError:
                    sys.exit('ERROR: Can\'t open ' + str(labelled_tree_file) + '!\n')

                for line in input:
                    line = line.strip()
                    content_of_labelled_tree_file += str(line)

                input.close()
                if not content_of_labelled_tree_file == content_of_previous_labelled_tree_file:
                    rooted_labelled_trees, insertion_point_node_hash = read_understand_and_reroot_the_labelled_tree(labelled_tree_file)
                    final_assingment_target_strings = Autovivify()
                new_assignments = Autovivify()
                at_least_one_new_assignment = 0
                try:
                    input = open(classification_file, 'r')
                except IOError:
                    sys.exit('ERROR: Can\'t open ' + str(classification_file) + '!\n')

                for line in input:
                    line = line.strip()
                    query, insertion_point_l, weight = line.split(' ')
                    assignment = ''
                    if re.search(r'I(\d+)', insertion_point_l):
                        assignment = re.search(r'I(\d+)', insertion_point_l).group(1)
                        nr_of_assignments += 1
                    assignments[assignment] = weight
                    if not assignment in final_assignment_target_strings.keys():
                        new_assignments[assignment] = 1
                        at_least_one_new_assignment = 1

                input.close()
                if at_least_one_new_assignment > 0:
                    prae_assignment_target_strings = identify_the_correct_terminal_children_of_each_assignment(\
                                                         terminal_children_strings_of_reference, rooted_labelled_trees,\
                                                         insertion_point_node_hash, new_assignments)

                    for assignment in sorted(prae_assignment_target_strings.keys()):
                        assignment_target_string = prae_assignment_target_strings[assignment]
                        final_assignment_target_strings[assignment] = assignment_target_string

            elif raxml_option == 'p':
                mp_tree_file = raxml_outfiles[denominator][f_contig]
                assignment = 'mp_root'
                assignments[assignment] = 1
                nr_of_assignments = 1
                prae_assignment_target_strings = get_correct_mp_assignment(terminal_children_strings_of_reference,\
                                                     mp_tree_file, assignments)
                assignment_target_string = prae_assignment_target_strings[assignment]
                final_assignment_target_strings[assignment] = assignment_target_string

            final_RAxML_filename = str(args.output_dir_raxml) + str(f_contig) + '_RAxML_parsed.txt'
            final_RAxML_output_files[denominator][final_RAxML_filename] = 1
            
            try:
                output = open(final_RAxML_filename, 'w')
            except IOError:
                sys.exit('ERROR: Can\'t create ' + str(final_RAxML_filename) + '!\n')
            output.write(str(description_text) + '\n')
            
            for assignment in sorted(assignments.keys()):
                assignment_target_string = final_assignment_target_strings[assignment]
                weight = assignments[assignment]
                relative_weight = int(((int(weight) / int(nr_of_assignments)) * 100) + 0.5)
                assignment_terminal_targets = assignment_target_string.split(' ')
                nr_of_terminal_targets = len(assignment_terminal_targets) - 1
                output.write('Placement weight ' + str(relative_weight) + '%: Assignment of query to ')
                if not nr_of_terminal_targets == 1:
                    output.write('the lowest common ancestor of ')
                count = 1

                while count <= nr_of_terminal_targets:
                    assignment_terminal_target = assignment_terminal_targets[count - 1]
                    is_last_element = 0
                    if count == nr_of_terminal_targets - 1:
                        is_last_element = 1
                    name_of_terminal_target = ''
                    name_of_terminal_target = tree_numbers_translation[denominator][assignment_terminal_target]
                    try:
                        name_of_terminal_target
                    except NameError:
                        sys.exit('ERROR: ' + str(assignment_terminal_target) + ' could not be located in the tree with the denominator ' +\
                                 str(denominator) + '!\n')
                    output.write(str(name_of_terminal_target) + ' (' + str(assignment_terminal_target) + ')')
                    if count < nr_of_terminal_targets - 1:
                        output.write(', ')
                    if count == nr_of_terminal_targets - 1:
                        output.write(' and ')
                    if count == nr_of_terminal_targets:
                        output.write('.')
                    count += 1

            output.close()
            content_of_previous_labelled_tree_file = content_of_labelled_tree_file

    return final_RAxML_output_files


def read_and_understand_the_reference_tree(reference_tree_file):
    reference_tree_elements = read_the_reference_tree(reference_tree_file)
    reference_tree_info = create_tree_info_hash()
    reference_tree_info = get_node_subtrees(reference_tree_elements, reference_tree_info)
    reference_tree_info = assign_parents_and_children(reference_tree_info)
    terminal_children_strings_of_reference = build_terminal_children_strings_of_reference_nodes(reference_tree_info)
    return terminal_children_strings_of_reference


def read_understand_and_reroot_the_labelled_tree(labelled_tree_file):
    labelled_tree_elements, insertion_point_node_hash = read_the_raxml_out_tree(labelled_tree_file)
    labelled_tree_info = create_tree_info_hash()
    labelled_tree_info = get_node_subtrees(labelled_tree_elements, labelled_tree_info)
    labelled_tree_info = assign_parents_and_children(labelled_tree_info)
    labelled_tree_info = build_tree_info_quartets(labelled_tree_info)
    rooted_labelled_trees = build_newly_rooted_trees(labelled_tree_info)
    return rooted_labelled_trees, insertion_point_node_hash


def identify_the_correct_terminal_children_of_each_assignment(terminal_children_strings_of_reference, rooted_labelled_trees, insertion_point_node_hash, assignments):
    terminal_children_strings_of_assignments = build_terminal_children_strings_of_assignments(rooted_labelled_trees, insertion_point_node_hash, assignments)
    real_terminal_children_strings_of_assignments = compare_terminal_children_strings(terminal_children_strings_of_assignments, terminal_children_strings_of_reference)
    return real_terminal_children_strings_of_assignments


def get_correct_mp_assignment(terminal_children_strings_of_reference, mp_tree_file, assignments):
    potential_terminal_children_strings = read_the_raxml_mp_out_tree(mp_tree_file, assignments)
    real_terminal_children_strings_of_assignments = compare_terminal_children_strings(potential_terminal_children_strings, terminal_children_strings_of_reference)
    return real_terminal_children_strings_of_assignments


def read_the_reference_tree(reference_tree_file):
    try:
        input = open(reference_tree_file, 'r')
    except IOError:
        sys.exit('ERROR: Could not open ' + reference_tree_file + '!\n')
    tree_string = ''

    for line in input:
        line = line.strip()
        tree_string += line

    input.close()

    tree_string = re.sub('\(', 'L', tree_string)
    tree_string = re.sub('\)', 'R', tree_string)
    tree_string = re.sub(r':\d+\.\d+', '', tree_string)
    count = -2

    while re.search('R', tree_string):
        tree_string = re.sub('R', 'Q' + str(count), tree_string, 1)
        count += -1

    tree_string = re.sub(r'Q-\d+;', 'Q;', tree_string)
    tree_string = re.sub('L', '(', tree_string)
    tree_string = re.sub('Q', ')', tree_string)
    reference_tree_elements = split_tree_string(tree_string)
    return reference_tree_elements


def read_the_raxml_out_tree(labelled_tree_file):
    insertion_point_node_hash = Autovivify()
    try:
        input = open(labelled_tree_file, 'r')
    except IOError:
        sys.exit('ERROR: Could not open ' + labelled_tree_file + '!\n')
    tree_string = ''

    for line in input:
        line = line.strip()
        tree_string += line

    input.close()
    tree_symbols_raw_1 = list(tree_string)
    bracket_diff = 0
    tree_string_neu = '('
    comma_count = 0

    for tree_symbol_raw_1 in tree_symbols_raw_1:
        if comma_count < 2:
            if tree_symbol_raw_1 == '(':
                bracket_diff += 1
            if tree_symbol_raw_1 == ')':
                bracket_diff += -1
            if tree_symbol_raw_1 == ',' and bracket_diff == 1:
                comma_count += 1
            if comma_count == 2:
                tree_string_neu += '):1.0[I666999666]'
        tree_string_neu += tree_symbol_raw_1

    tree_string = tree_string_neu
    tree_string = re.sub('\(', 'L', tree_string)
    tree_string = re.sub('\)', 'R', tree_string)
    tree_string = re.sub('\[', 'Q', tree_string)
    tree_string = re.sub(':1\.0', '', tree_string)
    
    while re.search(r'((\D(\d+))QI(\d+)])', tree_string):
        to_be_replaced = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(1)
        replacement = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(2)
        terminal_leaf = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(3)
        insertion_point = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(4)
        if terminal_leaf <= 0:
            sys.exit('ERROR: Your tree has terminal leaves with numbers <= 0. Please change them to positive values!\n')
        insertion_point_node_hash[insertion_point] = terminal_leaf
        tree_string = re.sub(to_be_replaced, replacement, tree_string)
    count = -2

    while re.search(r'QI(\d+)]', tree_string):
        insertion_point_node_hash[re.search(r'QI(\d+)]', tree_string).group(1)] = count
        tree_string = re.sub(r'QI(\d+)]', str(count), tree_string, 1)
        count += -1
    
    tree_string = re.sub('L', '(', tree_string)
    tree_string = re.sub('R', ')', tree_string)
    tree_string = re.sub('Q', '[', tree_string)
    tree_elements = split_tree_string(tree_string)
    return tree_elements, insertion_point_node_hash


def read_the_raxml_mp_out_tree(mp_tree_file, assignments):
    potential_terminal_children_strings = Autovivify()
    assignment = ''

    for assig in sorted(assignments.keys()):
        assignment = assig
        break

    try:
        input = open(mp_tree_file, 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + str(mp_tree_file) + '\n')
    tree_string = ''

    for line in input:
        line = line.strip()
        tree_string += line

    input.close()
    tree_string = re.sub('\(', 'L', tree_string)
    tree_string = re.sub('\)', 'R', tree_string)
    if not re.search(r',queryR;\Z', tree_string):
        sys.exit('ERROR: The query is not at the root of ' + str(mp_tree_file) + '!\n')
    else:
        tree_string = re.sub(r',queryR;\Z', 'R;', tree_string)
    tree_string = re.sub(r':\d+\.\d+', '', tree_string)
    count = -2

    while re.search('R', tree_string):
        tree_string = re.sub('R', 'Q' + str(count), tree_string, 1)
        count += -1

    tree_string = re.sub(r'Q-\d+;', 'Q;', tree_string)
    tree_string = re.sub('L', '(', tree_string)
    tree_string = re.sub('Q', ')', tree_string)
    tree_symbols = list(tree_string)
    bracket_diff = 0
    comma_count = 0
    substrings = ['', ',']

    for tree_symbol in tree_symbols:
        if comma_count < 1:
            if tree_symbol == '(':
                bracket_diff += 1
            if tree_symbol == ')':
                bracket_diff += -1
            if tree_symbol == ',' and bracket_diff == 1:
                comma_count += 1
            substrings[0] += tree_symbol
        else:
            substrings[1] += tree_symbol

    for substring in substrings:
        terminal_children = Autovivify()

        for eachGroup in re.findall(r'(\D)(\d+)', str(substring)):
            if eachGroup[0] == '-':
                continue
            terminal_children[eachGroup[1]] = 1

        potential_terminal_children_string = ''

        for potential_terminal_child in sorted(terminal_children.keys(), key=int):
            potential_terminal_children_string += str(potential_terminal_child) + ' '

        potential_terminal_children_strings[assignment][potential_terminal_children_string] = 1

    return potential_terminal_children_strings


def split_tree_string(tree_string):
    tree_symbols_raw = list(str(tree_string))
    count = -1
    previous_symbol = ''
    tree_elements = Autovivify()
    
    for tree_symbol_raw in tree_symbols_raw:
        if re.search(r'\d', tree_symbol_raw) and (re.search(r'\d', previous_symbol) or previous_symbol == '-'):
            tree_elements[count] += tree_symbol_raw
        else:
            count += 1
            tree_elements[count] = tree_symbol_raw
        previous_symbol = tree_symbol_raw
    return tree_elements


def create_tree_info_hash():
    tree_info = Autovivify()
    return tree_info


def get_node_subtrees(tree_elements, tree_info):
    bracket_l_count = 0
    bracket_r_count = 0
    parents_of_node = Autovivify()
    tree_element_nr = -1
    
    for tree_element in tree_elements.values():
        tree_element_nr += 1
        
        if str(tree_element) == '(':
            bracket_l_count = 1
            bracket_r_count = 0
            tree_sub_element_nr = tree_element_nr
            subtree_string = '('
            
            while True:
                tree_sub_element_nr += 1
                tree_sub_element = tree_elements[tree_sub_element_nr]
                if str(tree_sub_element) == '(':
                    bracket_l_count += 1
                if str(tree_sub_element) == ')':
                    bracket_r_count += 1
                if bracket_l_count == bracket_r_count:
                    nodename = tree_elements[tree_sub_element_nr + 1]
                    if str(nodename) == ';':
                        nodename = -1
                    subtree_string += ')' + str(nodename)
                    tree_info['subtree_of_node'][nodename] = subtree_string
                    break
                else:
                    subtree_string += str(tree_sub_element)
    
    for tree_element in tree_elements.values():
        if not re.search(r'\d+', str(tree_element)):
            continue
        if tree_element in tree_info['subtree_of_node'].keys():
            continue
        tree_info['subtree_of_node'][tree_element] = tree_element
    
    return tree_info


def assign_parents_and_children(tree_info):
    
    for node in sorted(tree_info['subtree_of_node'].keys()):
        if node == -1:
            continue
        subtree = str(tree_info['subtree_of_node'][node])
        parent = None
        
        for potential_parent in sorted(tree_info['subtree_of_node'].keys()):
            if node == potential_parent:
                continue
            potential_parent_subtree = str(tree_info['subtree_of_node'][potential_parent])
            subtree = re.sub('\(', 'L', subtree)
            subtree = re.sub('\)', '#', subtree)
            potential_parent_subtree = re.sub('\(', 'L', potential_parent_subtree)
            potential_parent_subtree = re.sub('\)', '#', potential_parent_subtree)
            potential_parent = str(potential_parent)
            if re.search(r'\AL'+re.escape(subtree)+r',.+#'+re.escape(potential_parent)+r'\Z', potential_parent_subtree) or \
               re.search(r'\AL.+,'+re.escape(subtree)+r'#'+re.escape(potential_parent)+r'\Z', potential_parent_subtree):
                parent = potential_parent
                break
        tree_info['parent_of_node'][node] = parent
        tree_info['children_of_node'][parent][node] = 1

    return tree_info


def build_tree_info_quartets(tree_info):
    for node in sorted(tree_info['parent_of_node'].keys(), key=int):
        parent = tree_info['parent_of_node'][node]
        if int(parent) == -1:
            for roots_child in sorted(tree_info['children_of_node']['-1'].keys(), key=int):
                if roots_child == node:
                    continue
                parent = roots_child

        tree_info['quartets'][node][parent] = 1
        if node in tree_info['children_of_node']:
            for child in sorted(tree_info['children_of_node'][node].keys(), key=int):
                tree_info['quartets'][node][child] = 1

    return tree_info


def build_newly_rooted_trees(tree_info):
    tree_number = 0
    list_of_already_used_attachments = Autovivify()
    rooted_trees = Autovivify()
    
    for node in sorted(tree_info['quartets'].keys(), key=int):
        if node in list_of_already_used_attachments:
            continue
        for attachment in sorted(tree_info['quartets'][node].keys(), key=int):
            list_of_already_used_attachments[attachment] = 1
            tree_string = ''
            root = -1
            node_infos = Autovivify()
            node_infos['previous_node'] = ''
            node_infos['node'] = ';'
            node_infos['open_attachments'][node] = 1
            node_infos['open_attachments'][attachment] = 1
            new_tree = recursive_tree_builder(tree_info, node_infos, tree_string)
            rooted_trees[tree_number] = new_tree
            tree_number += 1
    
    return rooted_trees


def recursive_tree_builder(tree_info, node_infos, tree_string):
    node = node_infos['node']
    count = 0

    for attachment in sorted(node_infos['open_attachments'].keys(), key=int):
        count += 1
        if count == 1:
            tree_string += '('
        node_infos2 = Autovivify()
        node_infos2['previous_node'] = node
        node_infos2['node'] = attachment
        count2 = 0

        for attachment_of_used_attachment in sorted(tree_info['quartets'][attachment].keys()):
            if attachment_of_used_attachment in node_infos['open_attachments']:
                continue
            if attachment_of_used_attachment == node:
                continue
            count2 += 1
            node_infos2['open_attachments'][attachment_of_used_attachment] = 1

        if count2 > 0:
            tree_string = recursive_tree_builder(tree_info, node_infos2, tree_string)
        else:
            tree_string += str(attachment)
        if count == 1:
            tree_string += ','
        if count == 2:
            tree_string += ')' + str(node)

    return tree_string


def build_terminal_children_strings_of_assignments(rooted_trees, insertion_point_node_hash, assignments):
    terminal_children_strings_of_assignments = Autovivify()

    for assignment in sorted(assignments.keys()):
        internal_node_of_assignment = insertion_point_node_hash[assignment]

        for rooted_tree in rooted_trees.keys():
            rooted_tree_elements = split_tree_string(rooted_trees[rooted_tree])
            rooted_tree_info = create_tree_info_hash()
            rooted_tree_info = get_node_subtrees(rooted_tree_elements, rooted_tree_info)
            assignment_subtree = str(rooted_tree_info['subtree_of_node'][str(internal_node_of_assignment)])
            terminal_children = Autovivify()
            if re.search(r'\A(\d+)\Z', assignment_subtree):
                terminal_children[re.search(r'\A(\d+)\Z', assignment_subtree).group(1)] = 1
            else:

                for each_hit in re.findall(r'(\D)(\d+)', assignment_subtree):
                    if each_hit[0] == '-':
                        continue
                    terminal_children[each_hit[1]] = 1

            terminal_children_string_of_assignment = ''

            for terminal_child_of_assignment in sorted(terminal_children.keys(), key=int):
                terminal_children_string_of_assignment += str(terminal_child_of_assignment) + ' '

            terminal_children_strings_of_assignments[assignment][terminal_children_string_of_assignment] = 1

    return terminal_children_strings_of_assignments


def build_terminal_children_strings_of_reference_nodes(reference_tree_info):
    terminal_children_strings_of_reference = Autovivify()

    for node in sorted(reference_tree_info['subtree_of_node'].keys()):
        reference_subtree = reference_tree_info['subtree_of_node'][node]
        terminal_children = Autovivify()
        if re.search(r'\A(\d+)\Z', str(reference_subtree)):
            terminal_children[re.search(r'\A(\d+)\Z', str(reference_subtree)).group(1)] = 1
        else:

            for each_hit in re.findall(r'(.)(\d+)', str(reference_subtree)):
                if each_hit[0] == '-':
                    continue
                terminal_children[each_hit[1]] = 1

        terminal_children_string_of_reference = ''

        for terminal_child_of_reference in sorted(terminal_children.keys(), key=int):
            terminal_children_string_of_reference += str(terminal_child_of_reference) + ' '

        terminal_children_strings_of_reference[terminal_children_string_of_reference] = 1

    return terminal_children_strings_of_reference


def compare_terminal_children_strings(terminal_children_strings_of_assignments, terminal_children_strings_of_reference):
    real_terminal_children_strings_of_assignments = Autovivify()
    there_was_a_hit = 0

    for assignment in sorted(terminal_children_strings_of_assignments.keys()):
        real_terminal_children_string = ''

        for terminal_children_string_of_assignment in sorted(terminal_children_strings_of_assignments[assignment].keys()):
            if terminal_children_string_of_assignment in terminal_children_strings_of_reference:
                real_terminal_children_string = terminal_children_string_of_assignment
                real_terminal_children_strings_of_assignments[assignment] = real_terminal_children_string
                there_was_a_hit = 1
                break

        if str(real_terminal_children_string) == '' and not str(assignment) == 'mp_root':
            sys.exit('ERROR: The RAxML output tree could not be rooted correctly!!!\n')

    if there_was_a_hit <= 0:
        sys.exit('ERROR: The RAxML output tree could not be rooted correctly!!!\n')
    return real_terminal_children_strings_of_assignments


def concatenate_RAxML_output_files(args, final_RAxML_output_files, text_of_analysis_type):
    output_directory_final = args.output_dir_final
    
    for denominator in sorted(final_RAxML_output_files.keys()):
        nr_of_files = 0
        assignments = Autovivify()
        description_text = '# ' + str(text_of_analysis_type[denominator]) + '\n'
        final_output_file_name = str(output_directory_final) + str(denominator) + '_concatenated_RAxML_outputs.txt'
        
        for final_RAxML_output_file in sorted(final_RAxML_output_files[denominator].keys()):
            nr_of_files += 1
            try:
                input = open(final_RAxML_output_file, 'r')
            except IOError:
                sys.exit('ERROR: Can\'t open ' + str(final_RAxML_output_file) + '!\n')
            
            for line in input:
                line = line.strip()
                if re.search(r'Placement weight (\d+)%: (.+)\Z', line):
                    weight = re.search(r'Placement weight (\d+)%: (.+)\Z', line).group(1)
                    assignment = re.search(r'Placement weight (\d+)%: (.+)\Z', line).group(2)
                    if assignment in assignments.keys():
                        assignments[assignment] += int(weight)
                    else:
                        assignments[assignment] = int(weight)
                else:
                    continue

            input.close()

        assignments_with_relative_weights = Autovivify()

        for assignment in sorted(assignments.keys(), reverse=True):
            weight = assignments[assignment]
            relative_weight = (int(((float(weight)/float(nr_of_files))*10000.0)+0.5))/10000.0
            assignments_with_relative_weights[relative_weight][assignment] = 1

        try:
            output = open(final_output_file_name, 'w')
        except IOError:
            sys.exit('ERROR: Can\'t create ' + str(final_output_file_name) + '!\n')

        if args.showmessages:
          print str(denominator) + '_ results concatenated:'

        output.write(str(description_text) + '\n')
        sum_of_relative_weights = 0

        for relative_weight in sorted (assignments_with_relative_weights.keys(), reverse=True):

            for assignment in sorted(assignments_with_relative_weights[relative_weight].keys(), reverse=True):
                sum_of_relative_weights += relative_weight
                if args.showmessages:
                   print 'Placement weight ' + str(relative_weight) + '%: ' + str(assignment)
                output.write('Placement weight ' + str(relative_weight) + '%: ' + str(assignment) + '\n')

        output.close()
        if args.showmessages:
           print str(denominator) + '_ sum of placement weights (should be 100): ' + str(sum_of_relative_weights)


def read_species_translation_files(args, cog_list):
    tree_numbers_translation = Autovivify()
    translation_files = Autovivify()
    phylogenetic_denominator = args.reftree
    if phylogenetic_denominator == 'g':
        translation_files[phylogenetic_denominator] = args.mltreemap_data + PATHDELIM + \
                                                      'data' + PATHDELIM + 'tree_data' + \
                                                      PATHDELIM + 'tax_ids_geba_tree.txt'
    elif phylogenetic_denominator == 'i':
        translation_files[phylogenetic_denominator] = args.mltreemap_data + PATHDELIM + \
                                                      'data' + PATHDELIM + 'tree_data' + \
                                                      PATHDELIM + 'tax_ids_fungitr.txt'
    else:
        translation_files[phylogenetic_denominator] = args.mltreemap_data + PATHDELIM + \
                                                      'data' + PATHDELIM + 'tree_data' + \
                                                      PATHDELIM + 'tax_ids_nr.txt'

    for functional_cog in sorted(cog_list['functional_cogs'].keys()):
        denominator = cog_list['functional_cogs'][functional_cog]
        filename = 'tax_ids_' + str(functional_cog) + '.txt'
        translation_files[denominator] = args.mltreemap_data + PATHDELIM + \
                                         'data' + PATHDELIM + 'tree_data' + \
                                         PATHDELIM + filename

    for phylogenetic_rRNA_cog in sorted(cog_list['phylogenetic_rRNA_cogs'].keys()):
        denominator = cog_list['phylogenetic_rRNA_cogs'][phylogenetic_rRNA_cog]
        filename = 'tax_ids_' + str(phylogenetic_rRNA_cog) + '.txt'
        translation_files[denominator] = args.mltreemap_data + PATHDELIM + \
                                         'data' + PATHDELIM + 'tree_data' + \
                                         PATHDELIM + filename

    for denominator in sorted(translation_files.keys()):
        filename = translation_files[denominator]
        try:
            input = open(filename, 'r')
        except IOError:
            sys.exit('ERROR: Can\'t open ' + str(filename) + '!\n')

        if args.showmessages:
          print 'opened ' + str(filename) # TK

        for line in input:
            line = line.strip()
            try:
                number, translation = line.split('\t')
            except ValueError:
                sys.exit('ValueError: .split(\'\\t\') on ' + str(line))
            tree_numbers_translation[denominator][number] = translation

        input.close()

    return tree_numbers_translation


def available_cpu_count():
    """ Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # cpuset
    # cpuset may restrict the number of *available* processors
    try:
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
                      open('/proc/self/status').read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    # http://code.google.com/p/psutil/
    try:
        import psutil
        return psutil.NUM_CPUS
    except (ImportError, AttributeError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])

        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0:
            return res
    except ImportError:
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                  stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)

        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')

        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        res = 0
        for pd in pseudoDevices:
            if re.match(r'^cpuid@[0-9]+$', pd):
                res += 1

        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0:
            return res
    except OSError:
        pass

    raise Exception('Can not determine number of CPUs on this system')


def deleteFiles(args):
    print 'Deleting the files other than \'final_outputs\' folder'
    sectionsToBeDeleted = []
    if args.delete:
        sectionsToBeDeleted = args.delete.split(':')

    filesToBeDeleted = []

    for section in sectionsToBeDeleted:
        if section == '1':
            filesToBeDeleted += glob.glob(args.output_dir_var + '*.fa')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*sequence.txt')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*sequence_shortened.txt')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*.fasta_formatted.txt')
        if section == '2':
            filesToBeDeleted += glob.glob(args.output_dir_var + '*BLAST_results*')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*blast_result_purified.txt')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*rRNA_result_summary.txt')
        if section == '3':
            filesToBeDeleted += glob.glob(args.output_dir_var + '*genewise.txt')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*genewise_result_summary.txt')
        if section == '4':
            filesToBeDeleted += glob.glob(args.output_dir_var + '*.mfa')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*.mfa-gb')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*.mfa-gb.txt')
        if section == '5':
            filesToBeDeleted += glob.glob(args.output_dir_var + '*_exit_after_Gblocks.txt')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*_RAxML.txt')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*RAxML_classification.txt')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*RAxML_info.txt')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*RAxML_labelledTree.txt')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*.phy')
            filesToBeDeleted += glob.glob(args.output_dir_var + '*.phy.reduced')

    for file in filesToBeDeleted:
        if path.exists(file):
            os.remove(file)
    
    if not args.verbose:
        dirsToBeDeleted = [args.output_dir_var, args.output_dir_raxml]
        for dir in dirsToBeDeleted:
            if path.exists(dir):
                shutil.rmtree(dir)
            else:
                pass


def main(argv, errorlogger = None, runstatslogger = None):
    global parser

    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = createParser()

    options, args = parser.parse_args(argv)

    options = checkParserArguments(options)

    removePreviousOutput(options)

    cog_list, text_of_analysis_type = createCogList(options)
    non_wag_cog_list = get_non_wag_cogs(options)
    splitFiles = splitFastaInput(options)

    # STAGE 2: Run BLAST to determine which COGs are present in the input sequence(s)
    runBlast(options, splitFiles)
    blastResults =  readBlastResults(options)
    blast_hits_purified = parseBlastResults(options, blastResults, cog_list)

    # STAGE 3: Run Genewise (or not) to produce amino acid sequences based on the COGs found in the input sequence(s)
    contig_coordinates, shortened_sequence_files = produceGenewiseFiles(options, blast_hits_purified)
    if options.reftype == 'n':
        genewise_outputfiles = startGenewise(options, shortened_sequence_files, blast_hits_purified)
        genewise_summary_files = parse_genewise_results(options, genewise_outputfiles, contig_coordinates)
        contig_rRNA_coordinates, rRNA_hit_files =  get_rRNA_hit_sequences(options, blast_hits_purified, cog_list, genewise_summary_files)
    elif options.reftype == 'a':
        genewise_summary_files = blastpParser(options, shortened_sequence_files, blast_hits_purified)

    # STAGE 4: Run hmmalign and Gblocks to produce the MSAs required to perform the subsequent ML/MP estimations
    hmmalign_singlehit_files  = prepare_and_run_hmmalign(options, genewise_summary_files, cog_list);
    concatenated_mfa_files, nrs_of_sequences, models_to_be_used = concatenate_hmmalign_singlehits_files(options, hmmalign_singlehit_files, non_wag_cog_list)
    gblocks_files  = start_gblocks(options, concatenated_mfa_files, nrs_of_sequences)

    # STAGE 5: Run RAxML to compute the ML/MP estimations
    phy_files = produce_phy_file(options, gblocks_files, nrs_of_sequences)
    raxml_outfiles, options2 = start_RAxML(options, phy_files, cog_list, models_to_be_used)
    tree_numbers_translation = read_species_translation_files(options, cog_list)
    final_RAxML_output_files = parse_RAxML_output(options, options2, tree_numbers_translation, raxml_outfiles, text_of_analysis_type)
    concatenate_RAxML_output_files(options, final_RAxML_output_files, text_of_analysis_type)

    # STAGE 6: Delete files as determined by the user
    deleteFiles(options)


def MetaPathways_mltreemap(argv, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\tMLTREEMAP_CALCULATION\n")
    createParser()
    main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger)
    return (0,'')

if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

