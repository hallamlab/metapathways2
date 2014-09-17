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
    standard_options_group.add_option("-o", "--output", dest="table_out",
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
    wtd_options_group.add_option("-n", "--ncbi-taxonomy-map", dest="ncbi_taxonomy_map",
        help='add the ncbi taxonomy map')
    wtd_options_group.add_option("--lca-min-score", dest="lca_min_score",  type='float', default=20,
        help='minimum BLAST/LAST score to consider as for LCA rule')
    wtd_options_group.add_option("--lca-top-percent", dest="lca_top_percent",  type='float', default=90,
        help='set of considered matches are within this percent of the highest score hit')
    wtd_options_group.add_option("--lca-min-support", dest="lca_min_support",  type='int', default=1,
        help='minimum number of reads that must be assigned to a taxon for ' + \
             'that taxon to be present otherwise move up the tree until there ' +
             'is a taxon that meets the requirement')

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
        if opts.ncbi_taxonomy_map == None:
            print "Need to specify annotation-table and ncbi-taxonomy-map for WTD"
            return False

    return True


def start_pathway_tools_api_mode(pathway_tools_exe):
    command = pathway_tools_exe + " -api"
    result = getstatusoutput(command)

# the main function
def main(argv):
    global parser
    (opts, args) = parser.parse_args()

    if not check_arguments(opts, args):
        print usage
        sys.exit(0)

    if opts.wtd and not os.path.isfile("/tmp/expected_taxa.txt"):
        # get expected taxonomic range(s)
        print 'Get Expected Taxonomic Range'
        cyc = PythonCyc()
        cyc.setOrganism('meta')
        cyc.setPToolsExec(opts.pathway_tools)
        cyc.startPathwayTools()

        pwys = cyc.getAllPathways()

        pwy_taxa_range = {} # hash from pwy to expected taxonomic range(s)
        for pwy in pwys:
            pwy_taxa_range[pwy] = cyc.getExpectedTaxonomicRange(pwy)

        cyc.stopPathwayTools()

    # get ORF to taxa map from annotation_table
    print "Getting ORF to Taxa Map from AnnotationTable"
    orf_lca = {}
    with open(opts.annotation_table) as f:
        for line in f:
            fields = line.split("\t")
            orf_lca[fields[0].strip()] = fields[8].strip()

    # get pathway ORFs and Rxns
    pwy_to_orfs = {}
    pwy_to_rxns = {}
    cyc = PythonCyc()
    cyc.setOrganism(opts.pgdb_name)
    cyc.setPToolsExec(opts.pathway_tools)
    cyc.startPathwayTools()
    pwys = cyc.getAllPathways()
    for pwy in pwys:
        genes = cyc.getPathwayORFs(pwy)
        rxns = cyc.getPathwayReactionInfo(pwy)
        pwy_to_orfs[pwy] = genes
        pwy_to_rxns[pwy] = rxns

    cyc.stopPathwayTools()

    # get LCA per pathway
    pwy_lca = {}
    # load NCBI taxonomy map
    print "Loading NCBI Taxonomy Map"
    lca = LCAComputation([ opts.ncbi_taxonomy_map ])
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

    for pwy in pwy_lca:
        min_taxa = "1"
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

        max_index, max_dist = max(enumerate(C), key=operator.itemgetter(1))
        max_taxa = C_taxa[max_index]

        # alternative min-magnitude
        if len(C_pos) > 0:
            min_mag_index, min_mag_dist = min(enumerate(C_pos), key=operator.itemgetter(1))
            min_mag_taxa = C_pos_taxa[min_mag_index]
        else:
            min_mag_index, min_mag_dist = max(enumerate(C_neg), key=operator.itemgetter(1))
            min_mag_taxa = C_neg_taxa[min_mag_index]

    # write-out pathway table



if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

