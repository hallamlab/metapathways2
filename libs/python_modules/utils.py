#!/usr/bin/env python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


"""Contains general utility code for the metapaths project"""

try:
    from shutil import rmtree
    from StringIO import StringIO
    from os import getenv, makedirs, path, remove
    from operator import itemgetter
    from os.path import abspath, exists, dirname, join, isdir
    from collections import defaultdict
    from optparse import make_option
    import sys
    import os
    import traceback
    from  python_modules.parsers.fastareader  import FastaReader
    from sysutil import pathDelim
    from metapaths_utils import printf, eprintf, fprintf
except:
    print "Cannot load some modules"
    sys.exit(0)
   

PATHDELIM = pathDelim()

def load_job_status_file(filename, A) :
       if path.exists(filename):
           listfile = open(filename, 'r')
           lines = listfile.readlines()
           listfile.close()
           for line in lines:
                fields = [ x.strip() for x in line.strip().split('\t') ]
                if len(fields) == 6:
                   if not fields[0] in  A:
                        A[fields[0]] = {}

                   if not fields[1] in  A[fields[0]]:
                        A[fields[0]][fields[1]] = {}

                   if not fields[2] in  A[fields[0]][fields[1]]:
                        A[fields[0]][fields[1]][fields[2]] = {}

                   if not fields[3] in  A[fields[0]][fields[1]][fields[2]]:
                        A[fields[0]][fields[1]][fields[2]][fields[3]] = {}

                   A[fields[0]][fields[1]][fields[2]][fields[3]][fields[4]]=int(fields[5])


def remove_files(dir, filenames):
   for file in filenames:
      try:
        if path.exists(dir + PATHDELIM + file):
         remove(dir + PATHDELIM + file)
      except IOError:
         print "Cannot remove file  " + dir + PATHDELIM + file + " !"
         sys.exit(0)


# (Re)create the sequence blocks along with the necessary log files 
def create_splits(outputdir, listfilename, input_filename, maxMBytes,   maxSize, splitPrefix = 'split', splitSuffix=''):
     maxBytes = 1024*1024*maxMBytes
     if splitSuffix:
        suffix = '.' + splitSuffix
     else:
        suffix = ''

     try:
        if path.exists( listfilename):
           listfile = open( listfilename, 'r')
           listfilenames = [ x.strip() for x in listfile.readlines() ]
           remove_files(outputdir, listfilenames)
           listfile.close()
     except IOError:
        print "Cannot read file " +  listfilename + " !"
        sys.exit(0)

     try:
        listfile = open(listfilename, 'w')
     except IOError:
        print "Cannot read file " + listfilename + " !"
        sys.exit(0)


     fragments= []
     seq_beg_pattern = re.compile(">")
     splitno = 0
     currblocksize = 0
     currblockbyteSize = 0

     fastareader = FastaReader(input_filename)
     # Read sequences from sorted sequence file and write them to block files

     for name in fastareader:
           fragments.append(fastareader.seqname) 
           fragments.append(fastareader.sequence)

           if currblocksize >= maxSize -1 or currblockbyteSize >= maxBytes:
               splitfile = open(outputdir +  PATHDELIM + splitPrefix + str(splitno) + suffix, 'w')
               fprintf(splitfile, "%s",'\n'.join(fragments))
               fragments=[]
               splitfile.close()
                # Add this block name to the blocklistfile
               fprintf(listfile, "%s\n", splitPrefix + str(splitno) + suffix)
               splitno += 1
               currblocksize = 0
               currblockbyteSize = 0
           else: 
               currblocksize += 1
               currblockbyteSize += len(fastareader.sequence)


     if fragments:
        splitfile = open(outputdir +  PATHDELIM + splitPrefix + str(splitno) + suffix, 'w')
        fprintf(splitfile, "%s",'\n'.join(fragments))
        splitfile.close()
        fragments = []
        fprintf(listfile, "%s\n", splitPrefix + str(splitno) + suffix)
        splitno += 1

     #Add this block name to the blocklistfile
     currblocksize = 0
     currblockbyteSize = 0

     listfile.close()
     return True

def countNoOfSequencesInFile(file):
    fastareader = FastaReader(file)
    count = 0
    for record in fastareader:
       count+=1
    return count

def number_of_lines_in_file(filename):
     
     try:  
         file = open(filename, 'r')
         lines = file.readlines()
         file.close()
         size = len(lines)
     except:   
         return 0
         
     return size

def  read_one_column(listfilename, dictionary, col=0) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields) > col:
          dictionary[fields[col]] = True
    listfile.close()
  except:
    traceback.print_exc(1)

def  enforce_number_of_fields_per_row(listfilename, col):
  needsSanitization = False
  try:
    listfile = open(listfilename, 'r+')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') if len(x.strip())  ]
       if len(fields) !=  col:
         needsSanitization = True

    if needsSanitization: 
       listfile.seek(0) 
       listfile.truncate()
    
       for line in lines:
          fields = [ x.strip() for x in line.strip().split('\t') if len(x.strip()) ]
          if len(fields) == col:
            fprintf(listfile, line)

    listfile.close()
  except:
    traceback.print_exc(1)

  return  needsSanitization


# if the the folder is found all the files
# in the folder and but DO NOT  delete the folder 
def clearFolderIfExists(folderName):
    if path.exists(folderName) :
       files = glob(folderName)
       for  f in files:
         remove(f)

# if the the folder is found all the files
# in the folder and then delete the folder too
def removeFolderIfFound(folderName):
    if path.exists(folderName) :
       files = glob(folderName)
       for  f in files:
         remove(f)
       if path.exists(folderName): 
         shutil.rmtree(origFolderName)


# if folder does not exist then create one
def createFolderIfNotFound( folderName ):
    if not path.exists(folderName) :
        makedirs(folderName)
        return False
    else:
        return True 
 
# does folder does ?
def doesFolderExist( folderName ):
    if not path.exists(folderName) :
        return False
    else:
        return True

# does file exist ?
def doesFileExist( fileName ):
    if not path.exists(fileName) :
        return False
    else:
        return True
#"""This module defines classes for working with GenBank records."""
import re
import sys



class FastaReader():
    """Parses a GenBank record from a string or file."""
    stop = False
    START_PATTERN = re.compile(r'^>')
    name = None
    future_name =None
    sequence=""
    def __init__(self, fasta_filename):
        try:
            self.file = open(fasta_filename, 'r')
        except IOError:
            print "Cannot open fasta file " + fasta_filename

    def __iter__(self):
        return self

 
    def next(self):
        if self.stop:
          raise StopIteration

        try:
           if not self.name: 
               self.name = self.file.readline().strip()
           line = self.file.readline().strip()
        except:
           line = None


        if not line:
           self.stop = True
           raise StopIteration


        fragments = []
        while line and not self.START_PATTERN.search(line):
            fragments.append(line.strip()) 
            line = self.file.readline()

       # print line
        if self.future_name:
            self.name = self.future_name

        if line:
          self.future_name = line.strip()

        self.sequence =''.join(fragments)
        self.seqname = self.name
        
        return self.name

# Read the contents of a file into a dictionary (col begin with 0)
def read_list(listfilename, dictionary, col=0) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields) > col:
        dictionary[fields[0]] = fields[col]
    listfile.close()
  except:
    traceback.print_exception()

def hasInput(expected_input):
    if  path.exists(expected_input):
        return True 
    else:
        return False

