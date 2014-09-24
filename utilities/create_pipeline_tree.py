#!python
"""
create_pipeline_tree.py

Created by Niels Hanson
Copyright (c) 2013 Steven J. Hallam Laboratory. All rights reserved.
"""

from __future__ import division

__author__ = "Niels W Hanson"
__copyright__ = "Copyright 2014"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Niels W Hanson"
__status__ = "Release"

try:
     import os
     import re
     import argparse
     import glob
     import random
     from os import makedirs, sys, remove
     from sys import path
except:
     print """ Could not load some modules """
     print """ """
     sys.exit(3)

what_i_do = "Creates an NCBI Tree Compatible with MetaPathways."
parser = argparse.ArgumentParser(description=what_i_do)

parser.add_argument('--ncbi_names', dest='ncbi_names', type=str, nargs='?',
          required=True, help='NCBI Taxonomy tree names.dmp file', default=None)
parser.add_argument('--ncbi_nodes', dest='ncbi_nodes', type=str, nargs='?',
          required=True, help='NCBI Taxonomy tree nodes.dmp file', default=None)                
parser.add_argument('-o', dest='output_file', type=str, nargs='?',
          required=True, help='name out output tree file', default=None)
parser.add_argument('--modifications', dest='modifications', type=str, 
          required=False, help='key-value file for lines apply modifications for MDM, prokaryotes, etc.', default=None)



def main(argv):
   args = vars(parser.parse_args())
   
   names = {}
   fh = open(args['ncbi_names'], "r")
   l = fh.readline()
   while l:
      fields = l.split("|")
      fields = map(str.strip, fields, "\t")
      fields = map(str.strip, fields)
      if fields[0] not in names:
          names[fields[0]] = []
      names[fields[0]].append(fields[1])
      l = fh.readline()
   fh.close()
   
   fh = open(args['ncbi_nodes'], "r")
   l = fh.readline()
   output = open(args['output_file'], "w")
   while l:
       fields = l.split("|")
       fields = map(str.strip, fields, "\t")
       child_id = fields[0]
       parent_id = fields[1]
       if child_id in names:
           while names[child_id]:
               output.write("\t".join([names[child_id].pop(), child_id, parent_id]) + "\n")
       else:
           print "Error: missing child"
       l = fh.readline()
   fh.close()
   output.close()
   
   # apply modifications
   if args['modifications']:
       mods = {}
       pattern = re.compile("\[\[(.*?)\]\]") # modification pattern

       with open(args['modifications'], "r") as mod_file:
           for line in mod_file:
               hits = pattern.findall(line)
               if hits:
                   mods[hits[0]] = hits[1]
       
       out = open(args['output_file'], "r")
       out2 = open(args['output_file'] + ".tmp", "w")
       for line in out:
           line = line.strip("\n")
           if line in mods:
               out2.write(mods[line] + "\n")
           else:
               out2.write(line + "\n")
   
       out.close()
   
       # add SAGs to nearest parent
       out2.write("\t".join(["candidate division KSB1 bacterium SCGC AAA255-G23", "99999998", "1046950"]) + "\n")
       out2.write("\t".join(["candidate division OP11 bacterium SCGC AAA011-L6", "99999997", "1046980"]) + "\n")
       out2.close()
       
       # rename file
       os.rename(args['output_file'] + ".tmp", args['output_file'])
       


# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
