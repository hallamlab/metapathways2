#!/usr/bin/env python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


"""Contains general utility code for the metapaths project"""

try:
    import os, re

    from shutil import rmtree
    from optparse import make_option
    from os import path, _exit, remove

    from libs.python_modules.utils.utils import *
except:
    print "Cannot load some modules"
    sys.exit(0)
   
PATHDELIM = pathDelim()


class Context:
    """ This class holds the context of a stage """
    inputs = {}
    inputs1 = {}
    outputs = {}
    outputs1 = {}
    name = None
    status = None
    commands = []
    message = "Message not set"

    def __init__(self):
        pass

    def isOutputAvailable(self):

        return doFilesExist(self.outputs.values())

    def isInputAvailable(self, errorlogger = None):
        #print self.inputs.values()
        status = True
        for file in self.inputs.values():
            if not doesFileExist(file):
          #      print file
                if errorlogger!=None:
                   errorlogger.printf("#STEP\t%s\n", self.name)
                   errorlogger.printf("ERROR\tMissing input %s\n", file)
                status = False
        return status



    def getMissingList(self, errorlogger = None):
        #print self.inputs.values()
        missingList = []
        status = True
        for file in self.inputs.values():
            if not doesFileExist(file):
                missingList.append(file)
                if errorlogger!=None:
                   errorlogger.printf("ERROR\tMissing input %s\n", file)
                status = False
        return missingList




    def removeOutput(self, errorlogger = None):
        #print self.inputs.values()
        annotationPATT = re.compile(r'annotation_table')
        for item in self.outputs.values():
           if not path.exists(item):
              continue

           if path.isdir(item):
              if annotationPATT.search(item):
                 pass
              else:
                 rmtree(item)
           else:
              remove(item)

