#!/usr/bin/python
# File created on Nov 27 Jan 2012
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

script_name = sys.argv[0]
usage= script_name + """ --subsys_file subsys_file --subsys2peg subsys2peg_file -o output"""
parser = OptionParser(usage)

parser.add_option( "--subsys_file", dest="subsys_file",  
                  help='the subsys.txt file, where the subsystems are listed')

parser.add_option( "--subsys2peg", dest="subsys2peg_file",  
                  help='the subsystem to peg file')

parser.add_option( "-o", dest="output_file",  
                  help='the output file')

def fprintf(file, fmt, *args):
    file.write(fmt % args)



def check_arguments(opts, args):
    if opts.subsys_file == None:
         print """Must have the \"subsys.txt\" file"""
         return False

    if opts.subsys2peg_file == None:
         print """Must have the \"subsys2peg\" file"""
         return False

    if opts.output_file == None:
         print """Must have an output file"""
         return False

    return True


def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

     return attributes


def create_dictionary(databasemapfile, annot_map):
       seq_beg_pattern = re.compile(">")

       dbmapfile = open( databasemapfile,'r')
       lines=dbmapfile.readlines()
       dbmapfile.close()
       for line in lines:
          if seq_beg_pattern.search(line):
              words = line.rstrip().split()
              name = words[0].replace('>','',1)
               
              words.pop(0)
              annotation = ' '.join(words)
              annot_map[name]= annotation
           


def process_subsys2peg_file(subsys2peg, peg2subsys,  org_file):
     try:
         orgfile = open(org_file,'r')
     except IOError:
         print "Cannot open " + str(org_file)

     lines = orgfile.readlines()
     for line in lines:
        hits = line.split('\t')
        newhits = []
        for hit in hits:
          if hit.strip():
             newhits.append(hit.strip())

        if len(newhits) > 2:
           if not newhits[1] in  subsys2peg:
              subsys2peg[newhits[1]]= []
           subsys2peg[newhits[1]].append(newhits[2])
           peg2subsys[newhits[2]]= newhits[1]





def process_function_file(catalogue, func_file):
     try:
         funcfile = open(func_file,'r')
     except IOError:
         print "Cannot open " + str(func_file)

     lines=funcfile.readlines()
     funcfile.close()

     for line in lines:
        hits = line.split("\t")
        newhits = []
        if len(hits) > 2:
          for hit in hits:
            if hit.strip():
               newhits.append(hit.strip())

        tempcatalogue = catalogue
        for i in  range(0, len(newhits) ):
           if not newhits[i] in tempcatalogue:
              tempcatalogue[newhits[i]]= {}
           tempcatalogue = tempcatalogue[newhits[i]]


def process_catalogue_file(catalogue, key,  subsys2peg, depth):

     prefixTab = ""
     for i in range(0, depth):
        prefixTab = prefixTab + "\t"

     if catalogue[key].keys():
         for skey in catalogue[key]:
            print prefixTab +  skey
            process_catalogue_file(catalogue[key], skey, subsys2peg, depth + 1)
#     else:
#         if key in subsys2peg:
#            for skey in subsys2peg[key]:
#               print prefixTab + skey
#         else:
#            print prefixTab +  key


def  process_whog_file(whog_file,whog_scheme) :
     square_brace = re.compile("\[([A-Z]*)\]\s*(\S*)\s*(.*)")
     map_pattern = re.compile("(.*):(.*)")

     try:
         whogfile = open(whog_file,'r')
     except IOError:
         print "Cannot open " + str(whog_file)

     lines=whogfile.readlines()
     whogfile.close()
     for line in lines:

        hits = square_brace.search(line)
        isCogLine = False
        if hits != None and len(hits.groups())==3:
          letter= hits.group(1).strip() 
          cogid = hits.group(2).strip() 
          function = hits.group(3).strip() 
          isCogLine = True

        if not isCogLine:  
           hits = map_pattern.search(line)
        else:
           hits = None

        isMapLine = False
        if hits!=None and len(hits.groups())==2:
           triplet= hits.group(1).strip()
           seqid= hits.group(2).strip()
           isMapLine = True


        if cogid!=None:
            if not cogid in whog_scheme:   
              whog_scheme[cogid]  = { "function":None, "seqmap":None, "letter":None }
              whog_scheme[cogid]["function"] = function
              whog_scheme[cogid]["letter"] = letter
              whog_scheme[cogid]["seqmap"] = {}

            if isMapLine:
               seqids = [ x.strip() for x in seqid.split() ]
               for s in seqids:
                   whog_scheme[cogid]["seqmap"][s]=triplet



def write_output(whog_scheme, letter_function_map, trip_organism_map, seqid_ginumber, output_file):
    
    try:
         outputfile = open(output_file,'w')
    except IOError:
         print "Cannot open output file" 

    for cogid in whog_scheme:
       for seqid in whog_scheme[cogid]['seqmap']:
          output = ""       
          letters = list(whog_scheme[cogid]["letter"]) 

          output+= ">" + seqid + " " + cogid 


          output+= " # Protein_GI_number: "
          if seqid in seqid_ginumber:
             output+= str(seqid_ginumber[seqid])
          else:
             output+= "00000000" 

          output+= " # Func_class: " + whog_scheme[cogid]["letter"] + " "

          letters = list(whog_scheme[cogid]["letter"]) 
          for x in letters:
            output+=letter_function_map[x] + " "

          output+=" # Function: "+ whog_scheme[cogid]["function"]
          output+=" # Organism: " + trip_organism_map[whog_scheme[cogid]['seqmap'][seqid]]
          fprintf(outputfile, "%s\n", output);
    outputfile.close()

def  process_ginumber_file(myva_gb_file) :
     try:
        myvagbfile = open(myva_gb_file,'r')
     except IOError:
        print "Cannot open " + str(myva_gb_file)

     lines=myvagbfile.readlines()
     myvagbfile.close()

     seqid_ginumber = {}
     for line in lines:
         words = [x.strip() for x in line.rstrip().split() ]
         if len(words)==2:
            seqid_ginumber[words[0]] = words[1]

     return seqid_ginumber


# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)
    catalogue = {}
    catalogue['SEED'] = {}
    process_function_file(catalogue["SEED"], opts.subsys_file)

    subsys2peg = {}
    peg2subsys = {}
    process_subsys2peg_file(subsys2peg, peg2subsys, opts.subsys2peg_file)
    process_catalogue_file(catalogue,'SEED',  subsys2peg, 0)

#    write_output(whog_scheme, letter_function_map, functionTriplets_organism_map, seqid_ginumber, opts.output_file)

    
    

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

