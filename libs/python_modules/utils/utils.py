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
    from glob import glob
    import sys, os, traceback, shutil

    from libs.python_modules.parsers.fastareader  import FastaReader
    from libs.python_modules.utils.sysutil import pathDelim
except:
    print "Cannot load some modules"
    sys.exit(0)
   

def fprintf(file, fmt, *args):
   file.write(fmt % args)
   
   

def printf(fmt, *args):
   sys.stdout.write(fmt % args)
   sys.stdout.flush()
 
def eprintf(fmt, *args):
   sys.stderr.write(fmt % args)
   sys.stderr.flush()

PATHDELIM = pathDelim()




def isFastaFile(filename):
    ''' this function checks if the given file is a fasta file
        by examining the first 100 lines of the file
    '''

    fastaNamePATT = re.compile(r'^>')
    fastaAlphabetPATT = re.compile(r'[a-zA-Z]+')
    isFasta = True
    seenNamePatt = False 

    try:
      c = 0
      with open(filename) as fp:
        for line in fp:
          '''trim the line'''
          line_trimmed = line.strip()
          
          if line_trimmed:
             if fastaNamePATT.search(line_trimmed):
               ''' is a name line '''
               seenNamePatt = True 
             else:
               ''' not a seq name '''
               if fastaAlphabetPATT.search(line_trimmed):
                  ''' it is of the alphabet'''
                  if not seenNamePatt:
                     ''' am i seeing sequence before the name'''
                     isFasta = False 
               else:
                  isFasta = False 
          c+=1
          if c > 500:
             break
      fp.close()
    except:
       eprintf("ERROR:\tCannot open filee " + filename)
       print traceback.print_exc(10)
       return False

    if seenNamePatt==False:
       isFasta = False
    
    return isFasta


def isGenbank(filename):
    ''' this function decides if a file is in genbank format or not
        by reading the first 100 lines and look for the key words that 
        usually appear in the genbank file formsts
    '''
    locusPATT = re.compile(r'^\s*LOCUS')
    versionPATT = re.compile(r'^\s*VERSION')
    featuresPATT = re.compile(r'^\s*FEATURES')
    originPATT = re.compile(r'\s*ORIGIN')
    accessionPATT = re.compile(r'^\s*ACCESSION')
    sourcePATT = re.compile(r'^\s*SOURCE')

    patterns = [locusPATT, versionPATT, featuresPATT, originPATT, accessionPATT, sourcePATT ]
    countPatterns  = [ 0 for i in range(0, len(patterns)) ]   

    try:
      c = 0
      with open(filename) as fp:
        for line in fp:
          '''trim the line'''
          line_trimmed = line.strip()
          if line_trimmed:
            for i in range(0, len(patterns) ):
               if patterns[i].search(line_trimmed.upper()):
                  countPatterns[i] = 1
          c+=1
          if c > 500:
            break

    except:
       eprintf("ERROR:\tCannot open filex " + filename)
       print traceback.print_exc(10)
       return False


    numPattsSeen = 0
    for val in countPatterns:
       numPattsSeen +=  val

    if numPattsSeen >= 3:
      '''if you have seen more than 3 of the above patters
         then we decide that it is a genbank file
      '''
      return True

    return False


def isNucleotide( filename):
    ''' checks if a fasta file is a nucleotide file format'''
    fastaNamePATT = re.compile(r'^>')
    isFasta = True
    nucCount = 0.0
    nonNucCount = 0.0

    try:
      c = 0
      with open(filename) as fp:
        for line in fp:
          '''trim the line'''
          line_trimmed = line.strip()
          if line_trimmed:
            if not fastaNamePATT.search(line_trimmed):
               for a in line_trimmed.upper():
                  if  a in ['A', 'T', 'C', 'G', 'N' ]:
                    nucCount+= 1
                  else:
                    nonNucCount+= 1
          c+=1
          if c > 500:
            break
    except:
       eprintf("ERROR:\tCannot open file " + filename)
       return False

    if nucCount ==0:
      return False

    if float(nucCount)/float(nonNucCount + nucCount) > 0.9 :
      return True

    return False


def check_file_types(filenames):
    filetypes={}

    for filename in filenames:
      if not path.exists(filename):
         filetypes[filename] = ['UNKNOWN', 'UNKNOWN', False]

      if isFastaFile(filename):
         if isNucleotide(filename):
            filetypes[filename] = ['FASTA', 'NUCL', False]
         else: # assume amino
            filetypes[filename] = ['FASTA', 'AMINO', False]
      elif isGenbank(filename):
            filetypes[filename] = ['GENBANK', 'NOT-USED', False]
      else:
         filetypes[filename] = ['UNKNOWN', 'UNKNOWN', False]

    return filetypes


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

def read_list(listfilename, dictionary, col=0) :
  """ Read the contents of a file into a dictionary (col begin with 0) """
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
    """ checks if the expected input, a file or folder is present"""
    if  path.exists(expected_input):
        return True 
    else:
        return False

def sQuote(string):
    """ Puts double quotes around a string"""
    return "\'" + string + "\'"



def shouldRunStep1(run_type, dir , expected_outputs):
    """ decide if a command should be run if it is overlay,
        when the expected outputs are present """
    if  run_type =='overlay'  and  doFilesExist(expected_outputs, dir =  dir):
        return False
    else:
        return True


def shouldRunStep(run_type, expected_output):
    """ decide if a command should be run if it is overlay,
      when results are alread computed decide not to run """
    if  run_type =='overlay'  and  path.exists(expected_output):
        return False
    else:
        return True

def hasResults(expected_output):
    """ has the results to use """
    if  path.exists(expected_output):
        return True
    else:
        return False


def hasResults1(dir , expected_outputs):
    """ has the results to use """
    if  doFilesExist(expected_outputs, dir =  dir):
        return True
    else:
        return False

def shouldRunStepOnDirectory(run_type, dirName):
    """if the directory is empty then there is not precomputed results
        and so you should decide to run the command
    """
    dirName = dirName + PATHDELIM + '*'
    files = glob(dirName)
    if len(files)==0:
      return True
    else:
      return False

def removeDirOnRedo(command_Status, origFolderName):
    """ if the command is "redo" then delete all the files
        in the folder and then delete the folder too """
    if command_Status=='redo' and path.exists(origFolderName) :
       folderName = origFolderName + PATHDELIM + '*'
       files = glob(folderName)
       for  f in files:
         remove(f)
       if path.exists(origFolderName):
         shutil.rmtree(origFolderName)

def removeFileOnRedo(command_Status, fileName):
    """ if the command is "redo" then delete the file """
    if command_Status=='redo' and path.exists(fileName) :
        remove(fileName)
        return True
    else:
        return False


def cleanDirOnRedo(command_Status, folderName):
    """ remove all the files in the directory on Redo """
    if command_Status=='redo':
       cleanDirectory(folderName)


def cleanDirectory( folderName):
    """ remove all the files in the directory """
    folderName = folderName + PATHDELIM + '*'
    files = glob(folderName)
    for  f in files:
       remove(f)

def checkOrCreateFolder( folderName ):
    """ if folder does not exist then create one """
    if not path.exists(folderName) :
        makedirs(folderName)
        return False
    else:
        return True

def doFilesExist( fileNames, dir="" ):
    """ does the file Exist? """
    for fileName in fileNames:
       file = fileName
       if dir!='':
         file = dir + PATHDELIM + fileName
       if not path.exists(file):
          return False
    return True


def Singleton(class_):
  instances = {}
  def getinstance(*args, **kwargs):
    if class_ not in instances:
        instances[class_] = class_(*args, **kwargs)
    return instances[class_]
  return getinstance


def extractSampleName(sampleName, type = None):
     sample_name  = sampleName 

     if type == 'fasta' or type==None:
         sample_name = re.sub(r'^.*/','',sample_name, re.I)
         sample_name = re.sub(r'^.*\\','',sample_name, re.I)
         sample_name = re.sub(r'\.fasta$','',sample_name, re.I)
         sample_name = re.sub(r'\.fna$','',sample_name, re.I)
         sample_name = re.sub(r'\.faa$','',sample_name, re.I)
         sample_name = re.sub(r'\.fas$','',sample_name, re.I)
         sample_name = re.sub(r'\.fa$','',sample_name, re.I)
     elif type in ['gbk-unannotated', 'gbk-annotated']  or type==None:
         sample_name = re.sub(r'^.*/','',sample_name, re.I)
         sample_name = re.sub(r'^.*\\','',sample_name, re.I)
         sample_name = re.sub(r'\.gbk$','',sample_name, re.I)
     else:
         eprintf("ERROR: Incorrect type %s to function extractSampleName\n", sQuote(type))

     return sample_name 


