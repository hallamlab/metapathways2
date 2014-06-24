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
    from optparse import make_option
except:
    print "Cannot load some modules"
    sys.exit(0)
   

class Tools():
    """Contains the acceptable configurations of the tool 
    If you want to add any new module, you will have to add the new algorithm 
    this configuration file."""

    tools = {}
    specialAlternatives =  { }



    def __init__(self):
        try:
            self.file = self.initializeTools()
        except IOError:
            print "ERROR : Cannot  initialize \"Configuration\" in file " + sys.argv[0]


    def isValid1(self, key1):
        if key1 in self.Tools:
            return True
        return False 

    def getValidValues1(self, key1):
        if self.isValid1(key1):
           return self.Tools.keys()
        return []


    def isValid2(self, key1, key2):
        if not self.isValid1(key1):
            return False

        if not key2 in  self.Tools[key1]:
            return False

        return True 

    def getValidValues2(self, key1, key2):
        if self.isValid2(key1, key2):
           return self.Tools[key1][key2].keys()
        return []

    def getTools(self):
        return self.Tools


    def getExecutables(self, step, params):
        choice = None
        result = {}
        if step in self.specialAlternatives:
           _keys = self.specialAlternatives[step]
           choice = self.getHashValue(params, _keys) 

        if choice == None:
           return self.Tools[step]['exec']
        else:
           #if there are choices
           return self.Tools[step]['exec'][choice.upper()]

        return {}


    def getHashValue(self, dictionary, keys): 
         _dictionary = dictionary
         value = None
         for key in keys: 
            if key in _dictionary:
               _dictionary = _dictionary[key]
               value = _dictionary

         return value


    def initializeTools(self):
        self.specialAlternatives =  { 
                                  'COMPUTE_REFSCORE'  : [ 'annotation', 'algorithm' ],
                                  'BLAST_REFDB'       : [ 'annotation', 'algorithm' ]
                               }

        self.Tools = { 
               'PREPROCESS_FASTA'  : {
                               'exec': { 'PREPROCESS_FASTA': None }
                                     }, 
               'ORF_PREDICTION'  : {
                              'exec':  {'ORF_PREDICTION': None }
                                   },
               'GFF_TO_AMINO':  {
                              'exec': { 'GFF_TO_FASTA': None }
                                }, 
               'FILTERED_FASTA': { 
                              'exec': { 'PREPROCESS_FASTA': None }
                                 },
               'COMPUTE_REFSCORE'  : {
                              'exec': { 
                                         'BLAST': {'COMPUTE_REFSCORE':None, 'BLASTP_EXECUTABLE' : None, 'FORMATDB_EXECUTABLE': None }, 
                                         'LAST':  {'COMPUTE_REFSCORE':None, 'LAST_EXECUTABLE' : None, 'LASTDB_EXECUTABLE':None }  
                                      },
                              },
               'BLAST_REFDB':{  
                              'exec': { 
                                         'BLAST': {'COMPUTE_REFSCORE':None, 'BLASTP_EXECUTABLE' : None, 'FORMATDB_EXECUTABLE': None }, 
                                         'LAST':  {'COMPUTE_REFSCORE':None, 'LAST_EXECUTABLE' : None, 'LASTDB_EXECUTABLE':None }  
                                      }
                               },
               'PARSE_BLAST':{
                              'exec': { 'PARSE_BLAST': None }
                             }, 
               'SCAN_rRNA'  : {
                              'exec': { 'SCAN_rRNA': None, 'BLASTN_EXECUTABLE':None }
                              }, 
               'SCAN_tRNA'  : {
                              'exec': { 'SCAN_tRNA': None }
                              }, 
               'ANNOTATE'  : {
                              'exec': { 'ANNOTATE': None }
                             }, 
               'PATHOLOGIC_INPUT'  : {
                              'exec': { 'GENBANK_FILE': None }
                                     }, 
               'GENBANK_FILE'  : {
                              'exec': { 'GENBANK_FILE': None }
                                 }, 
               'CREATE_SEQUIN_FILE'  : {
                              'exec': { 'CREATE_REPORT_FILES': None }
                                 }, 
               'CREATE_REPORT_FILES'  : {
                              'exec':  { 'CREATE_REPORT_FILES': None }
                                 }, 
               'MLTREEMAP_CALCULATION'  : {
                              'exec':  {'MLTREEMAP_CALCULATION': None }
                                          }, 
               'MLTREEMAP_IMAGEMAKER'  : {
                              'exec':  { 'MLTREEMAP_IMAGEMAKER': None }
                                         }, 
               'PATHOLOGIC': { 
                              'exec': {'PATHOLOGIC_EXECUTABLE': None }
                             }
           }





if __name__=="__main__":
     v =  Tools()
