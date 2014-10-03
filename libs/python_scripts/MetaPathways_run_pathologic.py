#!/usr/bin/python

"""This script run the pathologic """

try:
   import optparse, sys, re, csv, traceback
   from optparse import OptionGroup
   import pickle
   import math
   from libs.python_modules.taxonomy.LCAComputation import *
   import operator

   from os import path, _exit, remove, rename
   import logging.handlers
   from glob import glob
   from libs.python_modules.utils.sysutil import pathDelim
   from libs.python_modules.utils.metapathways_utils  import fprintf, printf, eprintf,  exit_process
   from libs.python_modules.utils.sysutil import getstatusoutput

   from libs.python_modules.utils.pathwaytoolsutils import *

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

def files_exist( files , errorlogger = None):
    status = True    
    for file in files:
       if not path.exists(file):
          if errorlogger:
             errorlogger.write( 'ERROR\tCould not find ptools input  file : ' +  file )
          status = False
    return not status



usage = sys.argv[0] + """ -i input_folder -p pgdb_dir --ptoolsExec pathwaytools_executable """
parser = None
def createParser():
    global parser

    epilog = """The pathway prediction algorithm, Pathologic in the Pathway Tools software, is run with the folder ptools as the input. The result of this step is an ePGDB (environmental pathway genome database).
The resulting ePGDB is in the ~/ptools-local/pgdbs/user folder. They can be viewed using the Pathway Tools software."""

    epilog = re.sub(r'\s+', ' ', epilog)

    parser = optparse.OptionParser(usage=usage, epilog = epilog)
    standard_options_group = OptionGroup(parser, "Standard Ptools group" )
    # Input options
    standard_options_group.add_option('-i', '--input', dest='inputfolder', default=None,
                           help='runs pathologic on the input folder')

    standard_options_group.add_option('-s', '--sample', dest='sample_name', default=None,
                           help='sample name')

    standard_options_group.add_option('-r', '--reactions', dest='reactions_list', default=None,
                           help='creates the metacyc reaction lists extracted from PGDB')

    standard_options_group.add_option('-p', '--pgdb', dest='pgdbdir', default=None,
                           help='folder of the PGDB')

    standard_options_group.add_option('--ptoolsExec', dest='ptoolsExec', default=None,
                           help='PathoLogic Executable')

    standard_options_group.add_option('--no-taxonomic-pruning', dest='no_taxonomic_pruning', default=True,
                           help='Option to stop taxonomic pruning')

    standard_options_group.add_option('--no-web-cel-overview', dest='no_web_cel_overview', default=True,
                           help='Option to turn off cellular overview')

    standard_options_group.add_option("-o", "--output-pwy-table", dest="table_out",
        help='the output table for the pathways [REQUIRED]')

    # WTD options [OPTIONAL]
    wtd_options_group = OptionGroup(parser, "Weighted Taxonomic Distance group")
    wtd_options_group.add_option("-d", "--wtd", dest="wtd", action="store_true",
        help='flag to add the WTD to each pathway')

    wtd_options_group.add_option("-a", "--annotation-table", dest="annotation_table",
        help='ORF annotation table for WTD')

    wtd_options_group.add_option("-n", "--ncbi-tree", dest="ncbi_tree",
        help='add the ncbi taxonomy map')

    wtd_options_group.add_option("--lca-min-score", dest="lca_min_score",  type='float', default=20,
        help='minimum BLAST/LAST score to consider as for LCA rule')

    wtd_options_group.add_option("--lca-top-percent", dest="lca_top_percent",  type='float', default=90,
        help='set of considered matches are within this percent of the highest score hit')

    wtd_options_group.add_option("--lca-min-support", dest="lca_min_support",  type='int', default=1,
        help='minimum number of reads that must be assigned to a taxon for ' + \
             'that taxon to be present otherwise move up the tree until there ' +
             'is a taxon that meets the requirement')

    wtd_options_group.add_option("--ncbi-megan-map", dest="ncbi_megan_map", help="MEGANs preferred mapping NCBI IDs" )


import os, signal
TIME = 10

def __StopPathwayTools():
    processPATT = re.compile(r'pathway-tools-runtime')
    for line in os.popen("ps xa"):
        fields = line.split()
        pid = fields[0]
        process = fields[4]
        result = processPATT.search(process)
        if result :
            os.kill(int(pid), signal.SIGHUP)


def StopPathwayTools():
  try:
     __StopPathwayTools()
     time.sleep(TIME)
     __StopPathwayTools()
     time.sleep(TIME)

     if path.exists("/tmp/ptools-socket"): 
        remove("/tmp/ptools-socket")
  except:
    pass


def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser

    options, args = parser.parse_args(argv)
    if options.inputfolder ==None:
       parser.error('ERROR\tInput folder for Pathologic not found')
    else:
      # required files to be able to build ePGDB
      files = [ 
                options.inputfolder + PATHDELIM + '0.pf',
                # options.inputfolder + PATHDELIM + '0.fasta',
                options.inputfolder + PATHDELIM + 'genetic-elements.dat',  
                options.inputfolder + PATHDELIM + 'organism-params.dat'
              ]

      if files_exist( files , errorlogger = errorlogger):
        exit_process("ERROR\tCannot find all inputs for Pathologic in folder %s : "  %(options.inputfolder) )

    # is there a pathwaytools executable installed
    if not path.exists(options.ptoolsExec):
       eprintf("ERROR\tPathwayTools executable %s not found!\n", options.ptoolsExec)
       if errorlogger:
          errorlogger.printf("ERROR\tPathwayTools executable %s not found!\n",  options.ptoolsExec)
       exit_process("ERROR\tPathwayTools executable %s not found!\n" %(options.ptoolsExec))


    # command to build the ePGDB
    command = "%s -patho %s"  %(options.ptoolsExec, options.inputfolder)
    if options.no_taxonomic_pruning:
       command += " -no-taxonomic-pruning "

    if options.no_web_cel_overview:
       command += " -no-web-cel-overview"

    command += " -api"

    status =0
    fix_pgdb_input_files(options.pgdbdir, pgdbs = [])

    if not path.exists(options.pgdbdir):
      status  = runPathologicCommand(runcommand = command) 
      fix_pgdb_input_files(options.pgdbdir, pgdbs = [])
    if status!=0:
       eprintf("ERROR\tFailed to run Pathologic on input %s : \n" %(options.inputfolder))
       eprintf("INFO\tKill any other PathwayTools instance running on the machine and try again\n")
       if errorlogger:
          errorlogger.write("ERROR\tFailed to run Pathologic on input %s : " %(options.inputfolder))
          errorlogger.write("INFO\tKill any other PathwayTools instance running on the machine and try again")
          errorlogger.write("     : " + command)
       exit_process("ERROR\tFailed to run Pathologic on input %s : "  %(options.inputfolder) )


    if not path.exists(options.reactions_list):
       try:
           pythonCyc = startPathwayTools(options.sample_name.lower(), options.ptoolsExec, True)
           pythonCyc.setDebug() # disable pathway debug statements
           printf("INFO\tExtracting the reaction list from ePGDB " + options.sample_name + "\n")
           resultLines = pythonCyc.getReactionListLines()
           #pythonCyc.stopPathwayTools()
           reaction_list_file = open(options.reactions_list + ".tmp", 'w')
           for line in resultLines:
              fprintf(reaction_list_file,"%s\n",line.strip())
           reaction_list_file.close()
           rename(options.reactions_list + ".tmp", options.reactions_list)

           StopPathwayTools()

       except:
           print traceback.print_exc(10)
           eprintf("ERROR\tFailed to run extract pathways for %s : \n" %(options.sample_name))
           eprintf("INFO\tKill any other PathwayTools instance running on the machine and try again")
           if errorlogger:
               errorlogger.write("ERROR\tFailed to run extract pathways for %s : " %(options.sample_name))
               errorlogger.write("INFO\tKill any other PathwayTools instance running on the machine and try again\n")
           StopPathwayTools()

    if not path.exists(options.table_out):
        ExtractPathway_WTD(options)
   

def startPathwayTools(organism, ptoolsExec, debug):
    StopPathwayTools()
    pythonCyc = PythonCyc()
    pythonCyc.setDebug(debug = debug)
    pythonCyc.setOrganism(organism)
    pythonCyc.setPToolsExec(ptoolsExec)
    pythonCyc.startPathwayTools()

    return pythonCyc


def  ExtractPathway_WTD(options):
    # Extract pathways and WTD
   # place to store list of expected taxonomic range(s)
    printf('INFO\tEntering the WTD calculations!\n')
    serialized_metacyc_taxa_ranges = "/tmp/metacyc_pwy_taxa_range.pk"
    serialized_metacyc_taxa_ranges_tmp = "/tmp/metacyc_pwy_taxa_range.pk.tmp"
    try:
        if options.wtd and not path.isfile(serialized_metacyc_taxa_ranges):
            # get MetaCyc's expected taxonomic range(s) and serialize for later use in /tmp
            # try:
            printf('INFO\tGetting MetaCyc Expected Taxonomic Range(s)\n')
            pythonCyc = startPathwayTools('meta', options.ptoolsExec, True)

            pwys = pythonCyc.getAllPathways()

            pwy_taxa_range = {} # hash from pwy to expected taxonomic range(s)
            pwy_taxa_range_pk = open(serialized_metacyc_taxa_ranges_tmp ,"w")

            # get expected taxonomic ranges for each pathway
            for pwy in pwys:
                # printf(" " + pwy)
                my_expected_taxonomic_range = pythonCyc.getExpectedTaxonomicRange(pwy)
                pwy_taxa_range[pwy] = my_expected_taxonomic_range
            # printf(" " + pwy)

            # write the pathway
            pickle.dump(pwy_taxa_range, pwy_taxa_range_pk)
            pwy_taxa_range_pk.close()
            StopPathwayTools()
            rename(serialized_metacyc_taxa_ranges_tmp, serialized_metacyc_taxa_ranges) 
        else:
            # read expected taxonomic range from serialized file
            exepected_taxa_in = open(serialized_metacyc_taxa_ranges ,"r")
            pwy_taxa_range = pickle.load(exepected_taxa_in)

        # create mapping of preferred NCBI to MEGAN taxonomy
        megan_map = {}
        if options.ncbi_megan_map:
            with open(options.ncbi_megan_map) as megan_map_file:
                for line in megan_map_file:
                    fields = line.split("\t")
                    fields = map(str.strip, fields)
                    megan_map[ fields[0] ] = fields[1]

        # get ORF to taxa map from annotation_table
        printf("INFO\tGetting ORF to Taxa Map from AnnotationTable\n")
        orf_lca = {}
        with open(options.annotation_table) as f:
            for line in f:
                fields = line.split("\t")
                orf_lca[fields[0].strip()] = fields[8].strip()

        # get pathway ORFs and Rxns
        pwy_to_orfs = {}
        pwy_to_long = {}
        pwy_to_rxns = {}
        try:
            pythonCyc = startPathwayTools(options.sample_name.lower(), options.ptoolsExec, True)
            pwys = pythonCyc.getAllPathways()

            for pwy in pwys:
                # printf(" " + pwy)
                genes = pythonCyc.getPathwayORFs(pwy)
                rxns = pythonCyc.getPathwayReactionInfo(pwy)
                pwy_to_orfs[pwy] = genes
                pwy_to_long[pwy] = cleanup(pythonCyc.get_slot_value(pwy, "common-name"))
                pwy_to_rxns[pwy] = rxns
            # printf("\n")
            StopPathwayTools()

        except:
            print """
            Problem connecting to Pathway Tools. Check the /tmp/ptools-socket file.
            """
    except:
        print """
        Problem calculating WTD via Pathway Tools. Check the /tmp/ptools-socket file.
        """

    # get LCA per pathway
    pwy_lca = {}
    # load NCBI taxonomy map
    printf("INFO\tLoading NCBI Taxonomy Map\n")
    lca = LCAComputation([ options.ncbi_tree ], )

    for pwy in pwy_to_orfs:
        orfs = pwy_to_orfs[pwy]
        taxa_ids = []
        for orf in orfs:
            if orf in orf_lca:
                # could strip out id here
                res = re.search("(.+?)\(([0-9]+?)\)",  orf_lca[orf] )
                if res:
                    taxa_annotation = res.group(1)
                    id = res.group(2)
                else:
                    id = lca.get_a_Valid_ID([ orf_lca[orf] ])
                taxa_ids.append(id)
        pwy_lca_id = lca.get_lca(taxa_ids, True)
        # print "In run_pathologic"
        # print pwy_lca_id
        # print pwy_lca_id
        lca.clear_cells(taxa_ids)

        pwy_lca[pwy] = [pwy_lca_id, lca.translateIdToName(pwy_lca_id)]

    # calculate weighted taxonomic distance
    pwy_to_wtd = {}
    for pwy in pwy_lca:

        C = [] # list of distances
        C_taxa = [] # list of parallel observed-expected taxa pairs
        C_pos = [] # list of non-negative distances
        C_pos_taxa = [] # list of parallel observed-expected taxa pairs
        C_neg = [] # list of negative distances
        C_neg_taxa = [] # list of parallel observed-expected taxa pairs

        if len(pwy_taxa_range[pwy]) > 0:
            for expected in pwy_taxa_range[pwy]:
                dist = lca.wtd(expected[0], pwy_lca[pwy][0])
                if dist or dist == 0:
                    # valid distance
                    # add distance respective lists
                    C.append(dist) # add distance
                    C_taxa.append([ expected[0], pwy_lca[pwy][0] ])
                    if dist >= 0:
                        C_pos.append(dist)  # add to non-negative list
                        C_pos_taxa.append([ expected[0], pwy_lca[pwy][0] ])
                    else:
                        C_neg.append(dist)  # add to negative list
                        C_neg_taxa.append([ expected[0], pwy_lca[pwy][0] ])
                else:
                    print "Not a valid distance"
                    continue
        else:
            # no expected taxonomy, set to root
            min_taxa = "1"
            dist = lca.wtd(min_taxa, pwy_lca[pwy][0])
            # add distance respective lists
            C.append(dist) # add distance
            C_taxa.append([ min_taxa, pwy_lca[pwy][0] ])
            if dist >= 0:
                C_pos.append(dist)  # add to non-negative list
                C_pos_taxa.append([ min_taxa, pwy_lca[pwy][0] ])
            else:
                C_neg.append(dist)  # add to negative list
                C_neg_taxa.append([ min_taxa, pwy_lca[pwy][0] ])

        # find index with max distance (closest to expected taxonomy)
        max_index, max_dist = max(enumerate(C), key=operator.itemgetter(1))
        max_taxa = C_taxa[max_index]

        # remap to preferred names
        observed = get_preferred_taxa_name(max_taxa[1], megan_map, lca.id_to_name)
        expected = get_preferred_taxa_name(max_taxa[0], megan_map, lca.id_to_name)

        pwy_to_wtd[pwy] = [ max_dist, observed, expected ]

    # write out pathway table
    table_out_tmp  = options.table_out + ".tmp"
    try:
        out = open(table_out_tmp, "w")
    except:
        print "Had problems opening file: " + options.table_out

    # write appropreate header
    if options.wtd:
        header = "SAMPLE\tPWY_NAME\tPWY_COMMON_NAME\tNUM_REACTIONS\tNUM_COVERED_REACTIONS\tORF_COUNT\tWTD\tOBSERVED\tEXPECTED\tORFS\n"
    else:
        header = "SAMPLE\tPWY_NAME\tPWY_COMMON_NAME\tNUM_REACTIONS\tNUM_COVERED_REACTIONS\tORF_COUNT\tORFS\n"
    out.write(header)

    sample = options.sample_name # sample name
    for pwy in pwy_to_orfs:
        # generate output line
        line = []
        line.append(sample) # sample name
        line.append(pwy) # pathway name
        line.append(pwy_to_long[pwy]) # pathway longname
        line.append(pwy_to_rxns[pwy][0]) # pathway num reactions
        line.append(pwy_to_rxns[pwy][1]) # pathway covered reactions
        line.append(len(pwy_to_orfs[pwy])) # num orfs
        if options.wtd:
            line.append(pwy_to_wtd[pwy][0]) # wtd
            line.append(pwy_to_wtd[pwy][1]) # wtd observed taxa
            line.append(pwy_to_wtd[pwy][2]) # wtd expected taxa
        line.append("[" + ",".join(pwy_to_orfs[pwy]) + "]") # list of ORFs

        line = map(str, line) # cast all to string

        out.write("\t".join(line) + "\n") # write out line
    try:
        out.close() # close file
        rename(table_out_tmp, options.table_out)
    except:
        print "Had problems closing file: " + options.table_out



def runPathologicCommand(runcommand = None):
    if runcommand == None:
      return False
    result = getstatusoutput(runcommand)
    return result[0]


# this is the portion of the code that fixes the name

def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

     return attributes


# this is the function that fixes the name
def  fix_pgdb_input_files(pgdb_folder, pgdbs = []):
     pgdb_list = glob(pgdb_folder + '/*/input/organism.dat')     

     for pgdb_organism_file in pgdb_list:
        process_organism_file(pgdb_organism_file)


def fixLine(line, id):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[0]+'\t' + id
     

def getID(line):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[1]
     
def process_organism_file(filel):
     patternsToFix = [ re.compile(r'NAME\tunclassified sequences'), re.compile(r'ABBREV-NAME\tu. sequences') ]
     patternID =  re.compile(r'^ID\t.*')
     try:
         orgfile = open(filel,'r')
     except IOError:
         print "ERROR : Cannot open organism file" + str(filel)
         return 

     lines = orgfile.readlines()
     newlines = []

     needsFixing = False

     id = None
     for line in lines:
         line = line.strip()
         if len(line)==0:
            continue
         flag = False

         result = patternID.search(line)
         if result:   
             id = getID(line)
          
         for patternToFix in patternsToFix:
             result = patternToFix.search(line)
             if result and id:
                 newline = fixLine(line, id)
                 newlines.append(newline)
                 flag= True
                 needsFixing = True

         if flag==False:
            newlines.append(line)

     orgfile.close()
     if needsFixing:
       write_new_file(newlines, filel)


def write_new_file(lines, output_file):
    
    print "Fixing file " + output_file 
    try:
       outputfile = open(output_file,'w')
       pass
    except IOError:
         print "ERROR :Cannot open output file "  + output_file
   
    for line in lines:
       fprintf(outputfile, "%s\n", line)

    outputfile.close()


def cleanup(string):
    """
    Cleans up pathway long-names for presentation.
    :param string:
    :return:
    """
    string = re.sub("|", "", string) # vertical bar
    string = re.sub("&", "", string) # ampersand
    string = re.sub(";", "", string) # semicolon
    string = re.sub("<[^<]+?>", '', string) # HTML tags
    string = re.sub("\'", "", string) # remove quotes

    return string

def get_preferred_taxa_name(taxa_id, megan_map, id_to_name):
    """
    Helper function to format NCBI IDs into preferred names. First checks for MEGAN name,
    if not found moves to current taxonomy in loaded NCBI taxonomy tree, failing that
    gives the taxonomy of 'Unknown', but still provides the id, e.g., 'Unknown (12345)'.
    :param taxa_id: numeric taxa id to translate
    :param megan_map: preferred megan mapping hash
    :param id_to_name: local ncbi tree hash
    :return: "perferred name (id)"
    """
    taxa_id = str(taxa_id)
    if taxa_id in megan_map:
        taxa = megan_map[ taxa_id ] + " (" + taxa_id + ")"
    elif taxa_id in id_to_name:
        taxa = id_to_name[ taxa_id ] + " (" + taxa_id + ")"
    else:
        taxa = "Unknown" + " (" + taxa_id + ")"

    return taxa

def MetaPathways_run_pathologic(argv, extra_command = None, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\tBUILD_PGDB\n")
    createParser()
    main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (0,'')

if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

