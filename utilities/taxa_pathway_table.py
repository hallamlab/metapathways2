#!/usr/bin/python
# File created on Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar, Niels W. Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     import re
     import sys
     from optparse import OptionParser, OptionGroup
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwayrc\""""
     print """ """
     sys.exit(3)

script_name = "taxa_pathway_table.py"
usage= script_name + """ --taxa taxa_file --path pathway_file -o output.csv"""
parser = OptionParser(usage)

parser.add_option( "--taxa", dest="taxa_file", help='the taxa file relating read-name to taxa name')

parser.add_option( "--path", dest="pathway_file", help='the pathway file listing pathways and associated read names')

parser.add_option( "--path_rxn", dest="path_rxn_file", help='a pathway and reaction file')

parser.add_option( "-o", dest="output_file", help='the output file')

# The C programming language is a way of life
def fprintf(file, fmt, *args):
    file.write(fmt % args)


# check the arguments
def check_arguments(opts, args):
    
    #TODO get this to work
    # if len(args) == 0:
    #     print script_name + ": a script to create pathway-taxa matrix of pathways from a read-taxa" 
    #     print "tab-separated file and the pathway summary file from the MetaPathways pipeline. Output"
    #     print "is a comma-separated matrix."
    #     return False
    
    if opts.taxa_file == None:
         print """Must have the \"taxa\" file"""
         return False
         
    if opts.pathway_file == None:
         print """Must have the \"path\" file"""
         return False
         
    if opts.output_file == None:
         print """Must have an output file"""
         return False
         
    return True

def process_read_taxa_file(taxa_file):
    # try to open taxa file
    try:
        mytaxa_file = open(taxa_file,'r')
    except IOError:
        print "Cannot open " + str(mytaxa_file)
    
    # read in the lines
    lines=mytaxa_file.readlines()
    mytaxa_file.close()
    
    read_to_taxa = {}
    for line in lines:
        words = [x.strip() for x in line.rstrip().split("\t") ]
        if len(words)==2:
            read_to_taxa[words[0]] = words[1]
            
    return read_to_taxa

def process_path_short_long(pathway_file):
    # try to open pathway file
    try:
        mypathway_file = open(pathway_file,'r')
    except IOError:
        print "Cannot open " + str(mypathway_file)
    
    # read in the lines
    lines=mypathway_file.readlines()
    mypathway_file.close()
    
    
    path_short_long = {} # short to long pathway names dictionary
    for line in lines:
        words = [x.strip() for x in line.rstrip().split("\t") ]
        path_short_long[words[0]] = re.sub("<.*?>","",words[1])
        
    return path_short_long

def process_pathways(pathway_file):
    # try to open pathway file
    try:
        mypathway_file = open(pathway_file,'r')
    except IOError:
        print "Cannot open " + str(mypathway_file)
    
    # read in the lines
    lines = mypathway_file.readlines()
    mypathway_file.close()
    
    # create a pathway object that contains the short_name, long_name, 
    # pathway_length, num_orfs, and orf_names
    pathways_object = {}
    for line in lines:
        words = [x.strip() for x in line.rstrip().split("\t") ]
        
        # collect all the information from the line
        short_name = words[0].strip() # MetaCyc short name
        long_name = words[1].strip() # long-name
        pathway_length = words[2].strip() # number of rxn in pathway
        num_orfs = words[3].strip() # number of ORFs assocated with pathway
        orf_names = [] # list of ORF names
        
        # collect all the ORF names on the line
        for orf in words[4:len(words)]:
            orf_names.append(orf.strip())
        
        # fill the object
        pathways_object[short_name] = { "long_name":None, \
                                        "pathway_length":None, \
                                        "num_orfs":None, \
                                        "orf_names":None }
        pathways_object[short_name]["long_name"] = long_name
        pathways_object[short_name]["pathway_length"] = pathway_length
        pathways_object[short_name]["num_orfs"] = num_orfs
        pathways_object[short_name]["orf_names"] = orf_names
        
    return pathways_object


def write_pathways_taxa(reads_to_taxa, pathways_short_to_long, pathways, output_file):
    # try to open output file
    try:
        myoutput_file = open(output_file,'w')
    except IOError:
        print "Cannot open " + str(myoutput_file)
        
    
    taxa_list = []
    for read in reads_to_taxa:
        if reads_to_taxa[read] not in taxa_list:
            taxa_list.append(reads_to_taxa[read])
    
    # count number of taxa for each pathway
    path_tax_table = {} # pathway and count of taxas to return
    for path in pathways:
        # set the taxa count to zero
        taxa_count = {}
        for taxa in taxa_list:
            taxa_count[taxa] = 0
        # for each orf name in the pathway
        for orf in pathways[path]["orf_names"]:
            if orf in reads_to_taxa:
                taxa_count[reads_to_taxa[orf]] += 1 # increment taxa count
        path_tax_table[path] = taxa_count
        
    # write taxa vs. pathway table
    # print header
    line = ""
    for taxa in taxa_list:
        line = line + "," + "\"" + taxa + "\""
    fprintf(myoutput_file, "%s\n", line)
    for path in path_tax_table:
        line = ""
        line = "\"" + path + "\""
        for taxa in taxa_list:
            line = line + "," + str(path_tax_table[path][taxa])
        fprintf(myoutput_file, "%s\n", line)
    
    myoutput_file.close()
    
    # pathway short_to_long file
    try:
        myoutput_file = open(output_file + str("short_to_long"),'w')
    except IOError:
        print "Cannot open " + str(myoutput_file)
    
    line = "short_name" + "\t" + "long_name" + "\t" + "pathway_length"
    fprintf(myoutput_file, "%s\n", line)
    for i in pathways_short_to_long:
        line = "\"" + i + "\"" +"\t"+ "\"" + pathways_short_to_long[i] + "\"" + "\t" + pathways[i]["pathway_length"]
        fprintf(myoutput_file, "%s\n", line)
    
    myoutput_file.close()
    

def process_pathway_reaction(path_rxn_file):
    # try to open pathway file
    try:
        mypath_rxn_file = open(path_rxn_file,'r')
    except IOError:
        print "Cannot open " + str(path_rxn_file)
    
    # read in the lines
    lines = mypath_rxn_file.readlines()
    mypath_rxn_file.close()
    
    # definition of a pathway line
    path_pattern = re.compile("^PATHWAY\:\t([^\t\n\r\f\v]+)\t([^\t\n\r\f\v]+)\t([0-9]+)\t([0-9]+)(.*)\n")
    rxn_pattern = re.compile("^RXN\:\t([^\t\n\r\f\v]+)\t([^\t\n\r\f\v]*)\t([0-9]*)\t{0,1}(.*)\n")
    
    # create a structure
    pathway_rxns = {}
    
    current_pathway = "" # takes the current pathway
    for line in lines:
        pathway_hits = path_pattern.search(line)
        rxn_hits = rxn_pattern.search(line)
        if pathway_hits:
            
            # pull the groups
            short_name = pathway_hits.group(1).strip()
            long_name = pathway_hits.group(2).strip()
            length = pathway_hits.group(3).strip()
            num_reads = pathway_hits.group(4).strip()
            reads = pathway_hits.group(5).strip()
            
            # create the pathway
            pathway_rxns[short_name] = { "long_name":None, "length":None, "num_reads":None, "rxns":None }
            pathway_rxns[short_name]["long_name"] = long_name
            pathway_rxns[short_name]["length"] = length
            pathway_rxns[short_name]["num_reads"] = num_reads
            pathway_rxns[short_name]["rxns"] = {}
            
            # set the current pathway
            current_pathway = short_name
            
        if rxn_hits:
            # parse reaction line
            rxn_short_name = rxn_hits.group(1).strip()
            rxn_long_name = rxn_hits.group(2).strip()
            rxn_num_reads = rxn_hits.group(3).strip()
            rxn_reads = rxn_hits.group(4).strip()
            
            # create reaction object
            pathway_rxns[current_pathway]["rxns"][rxn_short_name] = {"rxn_long_name":None, "rxn_num_reads":None, "rxn_reads":None, "taxa":{} }
            pathway_rxns[current_pathway]["rxns"][rxn_short_name]["rxn_long_name"] = rxn_long_name
            pathway_rxns[current_pathway]["rxns"][rxn_short_name]["rxn_num_reads"] = rxn_num_reads
            pathway_rxns[current_pathway]["rxns"][rxn_short_name]["rxn_reads"] = rxn_reads.split("\t")
            
    
    return pathway_rxns
    
def write_pathway_rxn_taxa_table(read_to_taxa,pathway_rxns,output_file):
    # try to open pathway file
    try:
        mypath_output_file = open(output_file,'w')
    except IOError:
        print "Cannot open " + str(output_file)
        
    
    
    for path in pathway_rxns:
        for rxn in pathway_rxns[path]["rxns"]:
            for read in pathway_rxns[path]["rxns"][rxn]["rxn_reads"]:
                if read in read_to_taxa:
                    cur_taxa = read_to_taxa[read]
                    if cur_taxa not in pathway_rxns[path]["rxns"][rxn]["taxa"]:
                        pathway_rxns[path]["rxns"][rxn]["taxa"] = {cur_taxa:1}
                    else:
                        pathway_rxns[path]["rxns"][rxn]["taxa"][cur_taxa] += 1
    
    line =""
    for path in pathway_rxns:
        fprintf(mypath_output_file, "%s\n", "-------")
        fprintf(mypath_output_file, "%s\n", str(path))
        fprintf(mypath_output_file, "%s\n", "-------")
        pathway_taxa_list = []
        for rxn in pathway_rxns[path]["rxns"]:
            for taxa in pathway_rxns[path]["rxns"][rxn]["taxa"]:
                if str(taxa) not in pathway_taxa_list:
                   pathway_taxa_list.append(str(taxa))
        line = "" # print header
        for i in pathway_taxa_list:
            line += "," + i
        fprintf(mypath_output_file, "%s\n", line)
        line = ""
        for rxn in pathway_rxns[path]["rxns"]:
            line += str(rxn)
            for i  in pathway_taxa_list:
                if i in pathway_rxns[path]["rxns"][rxn]["taxa"]:
                    line += "," + str(pathway_rxns[path]["rxns"][rxn]["taxa"][i])
                else:
                    line += "," + str(0)
            fprintf(mypath_output_file, "%s\n", line)
            line = ""
        
        
    mypath_output_file.close()


# the main function
def main(argv): 
    # grab and check arguments
    (opts, args) = parser.parse_args() 
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)
    
    # create a read to taxa dictionary
    read_to_taxa = process_read_taxa_file(opts.taxa_file)
    
    # create a MetaCyc pathway shortname to longname hash
    pathway_short_to_long =  process_path_short_long(opts.pathway_file)
    
    # collect pathway names, reads, and reactions
    pathways = process_pathways(opts.pathway_file)
    
    if(opts.path_rxn_file != None):
        pathway_rxns = process_pathway_reaction(opts.path_rxn_file)
        write_pathway_rxn_taxa_table(read_to_taxa,pathway_rxns,opts.output_file)
    
    
    # create table of pathways vs taxa in the output file
    # write_pathways_taxa(read_to_taxa, pathway_short_to_long, pathways, opts.output_file)
    
# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

