#!/usr/bin/python
# File created on Nov 27 Jan 2012

from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     from os import makedirs, sys, remove
     from sys import path
     import re
     from threading import Thread
     from time import sleep
     from optparse import OptionParser, OptionGroup

    # from libs.python_modules.utils.metapaths_utils  import parse_command_line_parameters, fprintf, printf
     from libs.python_modules.utils.sysutil import getstatusoutput
    # from sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)


usage= sys.argv[0] + """-o/--output table.txt -p/--pgdb pgdbname """

parser = None

def createParser():
    global parser
    parser = OptionParser(usage)
    parser.add_option("-o", "--output", dest="table_out", 
                  help='the output table for the pathways [REQUIRED]')

    parser.add_option("-p", "--pgdb", dest="pgdb_name", 
                  help='the pgdb name [REQUIRED]')

    parser.add_option("-t", "--ptools", dest="pathway_tools", 
                  help='the pathway tool executable [REQUIRED]')


def check_arguments(opts, args):
    if not hasattr(opts, 'table_out') == 0:
         print "Table file name  should be provided"  
         return False

    if not hasattr(opts, 'pgdb_name') == 0:
         print "PGDB name  should be provided"  
         return False

    if not hasattr(opts, 'pathway_tools') == 0:
         print "The pathway tools executable name should be provided"  
         return False

    return True


def start_pathway_tools_api_mode(pathway_tools_exe):
     command = pathway_tools_exe + " -api"
     result = getstatusoutput(command)
      
    
      

# the main function
def main(argv): 
    global parser
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)


    thread = Thread(target = start_pathway_tools_api_mode, args = (opts.pathway_tools, ))
    thread.start()
    sleep(10)
    print "going to kill "
    if thread.isAlive():
        try:
            thread._Thread__stop()
        except:
           print(str(thread.getName()) + ' could not be terminated')
    


# the main function
       

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

