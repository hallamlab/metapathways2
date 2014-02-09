#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division
import threading
import time 
from multiprocessing import Process

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
     import traceback
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwayrc\""""
     print """ """
     sys.exit(3)

script_name = "format_cog_names.py"
usage= script_name + """ --product product_file --subsystems subsystems_file -o output"""
parser = OptionParser(usage)

parser.add_option( "--products", dest="products",  
                  help='the product file')

parser.add_option( "--subsystems", dest="subsystems",  
                  help='the subsystem file')

parser.add_option( "-o", dest="output_file",  
                  help='the output file')

parser.add_option( "-N", dest="numthreads",  type = 'int',
                  help='the number of threads')



def fprintf(file, fmt, *args):
    file.write(fmt % args)



def check_arguments(opts, args):
    if opts.products == None:
         print """Must have the \"products\" file"""
         return False

    if opts.subsystems == None:
         print """Must have the \"subsystems\" file"""
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


def create_dict_from_list(list, dictionary):
    for word in list:
      dictionary[word] = True
   

stopwords_list= ["a", "able", "about", "across", "after", "all", "almost", "also", "am",\
                  "among", "an", "and", "any", "are", "as", "at", "be",\
                  "because", "been", "but", "by", "can", "cannot", "could", "dear", "did",\
                  "do", "does", "either", "else", "ever", "every", "for", "from", "get",\
                  "got", "had", "has", "have", "he", "her", "hers", "him", "his", "how", "however",\
                   "i", "if", "in", "into", "is", "it", "its", "just", "least", "let", "like",\
                   "likely", "may", "me", "might", "most", "must", "my", "neither", "no",\
                   "nor", "not", "of", "off", "often", "on", "only" , "or", "other", "our",\
                   "own", "rather", "said", "say", "says", "she", "should", "since", "so", "some",\
                   "than", "that", "the", "their", "them", "then", "there", "these", "they",\
                   "this", "tis", "to", "too", "twas", "us", "wants", "was", "we", "were", "what",\
                   "when", "where", "which", "while", "who", "whom", "why", "will", "with",\
                   "would", "yet", "you", "your", "enzyme", "hypothetical", "protein"]

def remove_repeats(filtered_words):
    word_dict = {}
    newlist = []
    for word in filtered_words:
       if not word in word_dict:
           word_dict[word]=1
           newlist.append(word)
    return  newlist

def remove_stop_words(list):
    newlist = []
    for item in list:
      if not item  in stopwords_list:
         newlist.append(item)
    return newlist

def format_product(product):
    product0 = re.sub(r'\[.+?\]', '', product)
    product1 = re.sub(r'\(', '', product0)
    product2 = re.sub(r'\)', '', product1)
    product3 = re.sub(r':', '', product2)
    product4 = re.sub(r'\/', '', product3)

    product_end = product4
 #######   subproducts = re.split('\/', product0)
    list = create_words_list(product_end) 

    dict =  create_dictionary_from_list(list)
    return dict


def create_dictionary_from_list(list):
    dict = {}
    for item in list:
       dict[item] = 1;
    return dict

   
def create_words_list(string):
     Splits = re.split(' ', string)
     list = []
     for Split in Splits:
        smallsplits = Split.split(',') 
        list +=smallsplits

     list = remove_repeats(list)
     list = remove_stop_words(list)
     
     list_final = [] 
     for item in list:
        if len(item):  
          list_final.append(item)
     return list_final  


def process_products_file(output_file, t,  N):
       outputfile = open( output_file,'w')
       seq_beg_pattern = re.compile(">")
       stopwords = {}
       create_dict_from_list(stopwords_list, stopwords)
       count = 0
       success = 0 

       print 'Thread :' + str(t) + ' subsystems :' + str(len(subsystems)) + ' ' + str(len(lines))

       #if t==0:
       #   print completed

       for line in lines:

          if t==10 and count % 1000 ==0:
             print count

          if count % N != t:
             count += 1
             continue
          
          count += 1
          if seq_beg_pattern.search(line):
              words = line.rstrip().split()
              seqname = re.sub('>', '', words.pop(0))
              product = ' '.join(words)
              productdict = format_product(product) 
              best_match, score = get_best_match(productdict, subsystems)
              if  best_match != None:
                fprintf(outputfile, ">%s\n", seqname  + '\t' + best_match)
                success +=1
              else:
                fprintf(outputfile, ">%s\n", seqname + '\t' + product)


       outputfile.close()

def get_best_match( productdict, subsystems = None, minMatchScore = 80):
     perfectmatch = 85
     maxScore = 0
     maxMatch = None
     maxFail= ( 1 - minMatchScore) + 1
     for key, value in subsystems.iteritems(): 
         score =  computescore(productdict, value, maxFail) 
         if score > maxScore and score > minMatchScore: 
             maxScore = score 
             maxMatch = key
             if maxScore > perfectmatch:
                break
             
     return maxMatch, maxScore

def computescore(productlist, value, maxFail) :
     matchSize = 0
     if len(productlist.keys()) > len(value.keys()):
         length = len(value.keys()) 
         maxFailNum = maxFail*length
         fails = 0

         for word in value:
            if word in productlist:
               matchSize += 1
            else:
               fails +=1
               if fails > maxFailNum:
                  return 0
     else:
         length =  len(productlist.keys())
         maxFailNum = maxFail*length
         fails = 0
         for word in productlist:
            if word in value:
               matchSize += 1
            else:
               fails +=1
               if fails > maxFailNum:
                  return 0

     if length == 0 :
        return 0

     return 100*matchSize/length
     


def  process_subsystems_file(subsystems_file, subsystems) :
     try:
        subsystemsfile = open(subsystems_file,'r')
     except IOError:
        print "Cannot open " + str(subsystems_file)

     sublines=subsystemsfile.readlines()
     subsystemsfile.close()

     for line in sublines:
         words = [x.strip() for x in line.rstrip().split('\t') ]
         maxField = len(words) - 1
         if maxField < 2:
           continue
         subsystems[words[maxField]] = create_dictionary_from_list(create_words_list(words[maxField]))


class ThreadClass(threading.Thread):
     N = None
     output_file = None
     subsystems = None
     #lines = None
     i = None
     #completed = [ 0 for i in range(20) ]
     def __init__(self,  output_file, i, N,  subsystems):
         threading.Thread.__init__(self)
         self.output_file = output_file
         self.i = i
         self.N = N
         self.subsystems = subsystems
       

     def run(self):
       pass  
       # process_products_file(self.output_file, self.i, self.N,  None)
  



# the main function

completed = [ 0 for i in range(0, 20)]
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)

    N = opts.numthreads 

    global subsystems
    global lines

    subsystems = {}
    lines = []
    process_subsystems_file(opts.subsystems, subsystems) 

    Threads = []
    try:
        productsfile = open( opts.products,'r')
        lines =  productsfile.readlines()
        productsfile.close()
    except:
        print 'cannot open ' + opts.products

#    print 'lines : ' + str(len(ThreadClass.lines))
    try:
       for i in range(0, N):
         print i,   opts.products,  opts.output_file + str(i),  N #, len(subsystems), len(ThreadClass.lines) 
         t =  Process(target =  process_products_file, args = (opts.output_file + str(i), i,  N, ) ) 
         Threads.append(t)
    except:
       print "Error: unable to start thread"

    for t in Threads:
       t.start()

    for t in Threads:
       t.join()

#    letter_function_map = process_function_file(opts.func_file)

#    functionTriplets_organism_map =  process_orginasim_file(opts.org_file)

#    whog_scheme = {}
#    process_whog_file(opts.whog_file, whog_scheme) 

#    seqid_ginumber = {}
#    if opts.myva_gb_file: 
#       seqid_ginumber= process_ginumber_file(opts.myva_gb_file) 

#    write_output(whog_scheme, letter_function_map, functionTriplets_organism_map, seqid_ginumber, opts.output_file)

    
    

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])


