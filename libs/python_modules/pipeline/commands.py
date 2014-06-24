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
   

class Command():

    inputs = []
    outputs = []
    commands = []
    flags = []

    def __init__(self):
        pass


class Commands():
    """Contains the list of commands for different steps of the pipeline """

    commands = {}


    def __init__(self):
        pass


    def addCommand(self, stepName, command):
        self.commands[stepName] = command
        return True 

if __name__=="__main__":
    pass
