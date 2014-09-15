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
    import re, sys, os, traceback
except:
    print "Cannot load some modules"
    sys.exit(0)
   

class Parameters():
    """Contains the acceptable values of the configuration 
    If you want to add any new module, you will have to add the new algorithm 
    this configuration file."""

    acceptableValues = {}


    def __init__(self):
        try:
            self.file = self.initializeExpectations()
        except IOError:
            print "ERROR : Cannot  initialize \"Configuration\" in file " + sys.argv[0]


    def isValid1(self, key1):
        if key1 in self.acceptableValues:
            return True
        return False 

    def getValidValues1(self, key1):
        if self.isValid1(key1):
           return self.acceptableValues.keys()
        return []


    def isValid2(self, key1, key2):
        if not self.isValid1(key1):
            return False

        if not key2 in  self.acceptableValues[key1]:
            return False

        return True 

    def getValidValues2(self, key1, key2):
        if self.isValid2(key1, key2):
           return self.acceptableValues[key1][key2].keys()
        return []

    def getAcceptableParameters(self):
        return self.acceptableValues 


    def getRunSteps(self, activeOnly = False):
        steps = []
        for step in self.acceptableValues['metapaths_steps']:
            if activeOnly==False  or self.acceptableValues['metapaths_steps'][step] in [ 'redo', 'yes' ]:
               steps.append(step)
        return steps




    def initializeExpectations(self):
        self.acceptableValues = { 
             'INPUT': { 
                        'format':{ 'fasta':True,
                                   'gbk-annotated': True,
                                   'gbk-unannotated':True, 
                                   'fasta-amino':True, 
                                 } 
                      },

              'quality_control': {
                        'delete_replicates': {
                                               'yes': True,
                                               'no': True 
                                             }
                      },

              'orf_prediction': {
                        'algorithm': {
                                       'prodigal':True
                                     }
                      },

              'annotation' : {
                         'algorithm' : { 
                                        'LAST': True ,
                                        'BLAST': True 
                                     }
                      },
# e.g. blast or last
#annotation:dbs metacyc-v4-2011-07-03,refseq-nr-2014-01-18,COG_2013-12-27,kegg-pep-2011-06-18,seed-2014-01-30

# rRNA annotation parameters
#rRNA:refdbs GREENGENES_gg16S-2012-11-06,LSURef_115_tax_silva,SSURef_NR99_115_tax_silva
# e.g. rRNA:refdbs GREENGENES_gg16S,SSURef_111_NR_tax_silva,LSURef_111_tax_silva

# pathway tools parameters
               'ptools_settings': {
                          'taxonomic_pruning' : {
                                         'yes': True, 
                                         'no': True 
                                         }
                       },
               'rRNA': { 'refdbs': { }
               }, 
# grid settings
               'metapaths_steps': {
                          'PREPROCESS_INPUT'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'ORF_PREDICTION'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'GFF_TO_AMINO'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'FILTER_AMINOS'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'COMPUTE_REFSCORES'  : {'yes':True, 'skip':True, 'redo':True, }, 
                          'FUNC_SEARCH'  : {'yes':True, 'skip':True, 'redo':True, 'grid':True }, 
                          'PARSE_FUNC_SEARCH'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'SCAN_rRNA'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'SCAN_tRNA'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'ANNOTATE_ORFS'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'PATHOLOGIC_INPUT'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'GENBANK_FILE'  : {'yes':True, 'skip':True, 'redo':True, }, 
                          'CREATE_ANNOT_REPORTS'  : {'yes':True, 'skip':True, 'redo':True }, 
                          'MLTREEMAP_CALCULATION'  : {'yes':True, 'skip': True, 'redo':True }, 
                          'BUILD_PGDB'  : {'yes':True, 'skip':True, 'redo':True, },
                          'COMPUTE_RPKM'  : {'yes':True, 'skip':True, 'redo':True, }
                }
           }


if __name__=="__main__":
     v =  Parameters()
     print v.isValid1('metapaths_steps')
     print v.isValid1('mtapaths_steps')
     print v.getValidValues1('metapaths_steps')
     print v.getValidValues2('INPUT', 'format')
