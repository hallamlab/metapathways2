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

script_name = "format_cog_names.py"
usage= script_name + """ --whog whog_file --org org.txt_file --fun func.txt [ --myva_gb gi_accession_file ] -o output"""
parser = OptionParser(usage)

parser.add_option( "--whog", dest="whog_file",  
                  help='the whog file where the classifications are listed')

parser.add_option( "--org", dest="org_file",  
                  help='the org where the Organism names are mentioned')

parser.add_option( "--fun", dest="func_file",  
                  help='the letter based function listing file ')

parser.add_option( "--myva_gb", dest="myva_gb_file",  
                  help='the mapping to GI number file [OPTIONAL] ')

parser.add_option( "-o", dest="output_file",  
                  help='the output file')

def fprintf(file, fmt, *args):
    file.write(fmt % args)



def check_arguments(opts, args):
    if opts.whog_file == None:
         print """Must have the \"whog\" file"""
         return False

    if opts.org_file == None:
         print """Must have the \"org\" file"""
         return False

    if opts.func_file == None:
         print """Must have the \"func.txt\" file"""
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
           


def process_orginasim_file(org_file):
     organism_name = re.compile("(\S*)\s*\S*\s*\S*\s*(.*)")
     try:
         orgfile = open(org_file,'r')
     except IOError:
         print "Cannot open " + str(org_file)

     functionTriplets_organism_map = {}
     lines=orgfile.readlines()
     orgfile.close()
     for line in lines:
        hits = organism_name.search(line)
        if hits:
#           print hits.group(1) + " " + hits.group(2)
           functionTriplets_organism_map[hits.group(1)] = hits.group(2).strip()

     return functionTriplets_organism_map




def process_function_file(func_file):
     square_brace = re.compile("\[([A-Z])\](.*)")
     try:
         funcfile = open(func_file,'r')
     except IOError:
         print "Cannot open " + str(func_file)

     letter_function_map = {}
     lines=funcfile.readlines()
     funcfile.close()
     for line in lines:
        hits = square_brace.search(line)
        if hits:
           letter_function_map[hits.group(1)] = hits.group(2).strip()

     return letter_function_map


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


    letter_function_map = process_function_file(opts.func_file)

    functionTriplets_organism_map =  process_orginasim_file(opts.org_file)

    whog_scheme = {}
    process_whog_file(opts.whog_file, whog_scheme) 

    seqid_ginumber = {}
    if opts.myva_gb_file: 
       seqid_ginumber= process_ginumber_file(opts.myva_gb_file) 

    write_output(whog_scheme, letter_function_map, functionTriplets_organism_map, seqid_ginumber, opts.output_file)

    
    

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

