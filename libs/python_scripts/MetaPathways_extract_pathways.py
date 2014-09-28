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
    import os
    from sys import path
    import re
    import operator
    import math
    import pickle
    from threading import Thread
    from time import sleep
    from optparse import OptionParser, OptionGroup
    from libs.python_modules.taxonomy.LCAComputation import *

    # from libs.python_modules.utils.metapaths_utils  import parse_command_line_parameters, fprintf, printf
    from libs.python_modules.utils.pathwaytoolsutils import PythonCyc
    from libs.python_modules.utils.sysutil import getstatusoutput
except:
    print """ Could not load some user defined  module functions"""
    print """ Make sure your typed \"source MetaPathwaysrc\""""
    print """ """
    sys.exit(3)


usage= sys.argv[0] + """-o/--output table.txt -p/--pgdb pgdbname """

parser = None

def createParser():
    global parser
    parser = OptionParser(usage)

    # standard options [REQUIRED]
    standard_options_group = OptionGroup(parser, "Input/Output/Ptools group" )
    standard_options_group.add_option("-o", "--output-pwy-table", dest="table_out",
        help='the output table for the pathways [REQUIRED]')
    standard_options_group.add_option("-p", "--pgdb", dest="pgdb_name",
        help='the pgdb name [REQUIRED]')
    standard_options_group.add_option("-t", "--ptools", dest="pathway_tools",
        help='the pathway tool executable [REQUIRED]')

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
    wtd_options_group.add_option("--ncbi-megan-map", dest="ncbi_megan_map", help="MEGANs prefered mapping NCBI IDs" )

def check_arguments(opts, args):
    # standard options
    if opts.table_out == None:
        print "Table file name should be provided"
        return False
    if opts.pgdb_name == None:
        print "PGDB name should be provided"
        return False
    if opts.pathway_tools == None:
        print "The pathway tools executable name should be provided"
        return False

    # wtd options
    if hasattr(opts, 'wtd'):
        if opts.annotation_table == None:
            print "Need to specify annotation-table and ncbi-taxonomy-map for WTD"
            return False
        if opts.ncbi_tree == None:
            print "Need to specify annotation-table and ncbi-taxonomy-map for WTD"
            return False

    return True


def start_pathway_tools_api_mode(pathway_tools_exe):
    command = pathway_tools_exe + " -api"
    result = getstatusoutput(command)

def cleanup(string):
    """
    Cleans up pathway long-names for presentation.
    :param string:
    :return:
    """
    string = re.sub("|", "", string) # vertical bar
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

# the main function
def main(argv):
    global parser
    (opts, args) = parser.parse_args()

    if not check_arguments(opts, args):
        print usage
        sys.exit(0)

    # place to store list of expected taxonomic range(s)
    serialized_metacyc_taxa_ranges = "/tmp/metacyc_pwy_taxa_range.pk"

    if opts.wtd and not os.path.isfile(serialized_metacyc_taxa_ranges):
        # get MetaCyc's expected taxonomic range(s) and serialize for later use in /tmp
        try:
            print 'Getting MetaCyc Expected Taxonomic Range(s)'

            # connect to Pathway Tools
            cyc = PythonCyc()
            cyc.setOrganism('meta')
            cyc.setPToolsExec(opts.pathway_tools)
            cyc.startPathwayTools()

            pwys = cyc.getAllPathways()

            pwy_taxa_range = {} # hash from pwy to expected taxonomic range(s)
            pwy_taxa_range_pk = open(serialized_metacyc_taxa_ranges ,"w")

            # get expected taxonomic ranges for each pathway
            for pwy in pwys:
                my_expected_taxonomic_range = cyc.getExpectedTaxonomicRange(pwy)
                pwy_taxa_range[pwy] = my_expected_taxonomic_range

            # write the pathway
            pickle.dump(pwy_taxa_range, pwy_taxa_range_pk)
            pwy_taxa_range_pk.close()

            # close Pathway Tools
            cyc.stopPathwayTools()
        except:
            print """
            Problem connecting to Pathway Tools. Check the /tmp/ptools-socket file.
            """
    else:
        # read expected taxonomic range from serialized file
        exepected_taxa_in = open(serialized_metacyc_taxa_ranges ,"r")
        pwy_taxa_range = pickle.load(exepected_taxa_in)

    # create mapping of preferred NCBI to MEGAN taxonomy
    megan_map = {}
    if opts.ncbi_megan_map:
        with open(opts.ncbi_megan_map) as megan_map_file:
            for line in megan_map_file:
                fields = line.split("\t")
                fields = map(str.strip, fields)
                megan_map[ fields[0] ] = fields[1]

    # get ORF to taxa map from annotation_table
    print "Getting ORF to Taxa Map from AnnotationTable"
    orf_lca = {}
    with open(opts.annotation_table) as f:
        for line in f:
            fields = line.split("\t")
            orf_lca[fields[0].strip()] = fields[8].strip()

    # get pathway ORFs and Rxns
    pwy_to_orfs = {}
    pwy_to_long = {}
    pwy_to_rxns = {}
    try:
        cyc = PythonCyc()
        cyc.setOrganism(opts.pgdb_name)
        cyc.setPToolsExec(opts.pathway_tools)
        cyc.startPathwayTools()
        pwys = cyc.getAllPathways()
        for pwy in pwys:
            genes = cyc.getPathwayORFs(pwy)
            rxns = cyc.getPathwayReactionInfo(pwy)
            pwy_to_orfs[pwy] = genes
            pwy_to_long[pwy] = cleanup(cyc.get_slot_value(pwy, "common-name"))
            pwy_to_rxns[pwy] = rxns

        cyc.stopPathwayTools()
    except:
        print """
        Problem connecting to Pathway Tools. Check the /tmp/ptools-socket file.
        """

    # get LCA per pathway
    pwy_lca = {}
    # load NCBI taxonomy map
    print "Loading NCBI Taxonomy Map"
    lca = LCAComputation([ opts.ncbi_tree ])
    lca.setParameters(opts.lca_min_score, opts.lca_top_percent, opts.lca_min_support)

    for pwy in pwy_to_orfs:
        orfs = pwy_to_orfs[pwy]
        taxa_ids = []
        for orf in orfs:
            if orf in orf_lca:
                id = lca.get_a_Valid_ID([ orf_lca[orf] ])
                taxa_ids.append(id)
        pwy_lca_id = lca.get_lca(taxa_ids, True)
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
        if pwy in pwy_taxa_range:
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
    try:
        out = open(opts.table_out, "w")
    except:
        print "Had problems opening file: " + opts.table_out

    # write appropreate header
    if opts.wtd:
        header = "SAMPLE\tPWY_NAME\tPWY_COMMON_NAME\tNUM_REACTIONS\tNUM_COVERED_REACTIONS\tORF_COUNT\tWTD\tOBSERVED\tEXPECTED\tORFS\n"
    else:
        header = "SAMPLE\tPWY_NAME\tPWY_COMMON_NAME\tNUM_REACTIONS\tNUM_COVERED_REACTIONS\tORF_COUNT\tORFS\n"
    out.write(header)

    sample = opts.pgdb_name # sample name
    for pwy in pwy_to_orfs:
        # generate output line
        line = []
        line.append(sample) # sample name
        line.append(pwy) # pathway name
        line.append(pwy_to_long[pwy]) # pathway longname
        line.append(pwy_to_rxns[pwy][0]) # pathway num reactions
        line.append(pwy_to_rxns[pwy][1]) # pathway covered reactions
        line.append(len(pwy_to_orfs[pwy])) # num orfs
        if opts.wtd:
            if pwy in pwy_to_wtd:
                line.append(pwy_to_wtd[pwy][0]) # wtd
                line.append(pwy_to_wtd[pwy][1]) # wtd observed taxa
                line.append(pwy_to_wtd[pwy][2]) # wtd expected taxa
            else:
                line.append("NA")
                line.append("NA")
                line.append("NA")
        line.append("[" + ",".join(pwy_to_orfs[pwy]) + "]") # list of ORFs

        line = map(str, line) # cast all to string

        out.write("\t".join(line) + "\n") # write out line
    try:
        out.close() # close file
    except:
        print "Had problems closing file: " + opts.table_out

if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

