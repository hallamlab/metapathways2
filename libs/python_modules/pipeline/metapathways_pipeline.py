#!/usr/bin/env python
# File created on 30 Dec 2009.


from __future__ import division
try:
      import traceback
      import sys
      from subprocess import Popen, PIPE, STDOUT
      from os import makedirs, listdir, _exit
      from glob import glob
      from optparse import OptionParser
      from os.path import split, splitext, join, dirname, abspath
      from datetime import datetime

      from libs.python_modules.utils.metapathways_utils import printf, eprintf
      from libs.python_modules.utils.sysutil import getstatusoutput
      from libs.python_modules.utils.pathwaytoolsutils import *

      from libs.python_modules.grid.BlastGrid import *
      import libs.python_scripts  as python_scripts

except:
      print """ Could not load some user defined  module functions"""
      print """ Make sure your typed \"source MetaPathwaysrc\""""
      print """ """
      print traceback.print_exc(10)
      sys.exit(3)


__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


"""
This file contains the metapaths workflow functions which string together 
independent scripts. 
"""

def print_commands(commands,
                   status_update_callback,
                   logger):
    """Print list of commands to run """
    #logger.write("Printing commands only.\n\n")
    #for c in commands:
    #    for e in c:
    #        status_update_callback('#%s' % e[0])
    #        print '%s' % e
    #        logger.write('# %s command\n%s\n\n' % e)
            
def execute_pipeline_stage(pipeline_command, extra_command = None,  errorlogger = None, runstatslogger = None):

     argv = [ x.strip() for x in pipeline_command.split() ]

     funcname = re.sub(r'.py$','', argv[0])
     funcname = re.sub(r'^.*/','', funcname)
     args = argv[1:] 
     

     if hasattr(python_scripts, funcname):
        methodtocall = getattr( getattr(python_scripts, funcname), funcname)
        if extra_command == None:
           result = methodtocall(args, errorlogger = errorlogger, runstatslogger = runstatslogger)
        else:
#           print extra_command
           result = methodtocall(args, errorlogger = errorlogger, extra_command = extra_command, runstatslogger = runstatslogger)
     else:
        result = getstatusoutput(pipeline_command)
     return result



def execute_tasks(s, verbose = False, block = 0):
    """Run list of commands, one after another """
    #logger.write("Executing commands.\n\n")
    contextBlocks = s.getContextBlocks()
       
    contextBlock = contextBlocks[block]

    for c in contextBlock:
        #print c.name, c.status, 'status'
        if c.status=='stop':
           print "Stopping!"
           s.stepslogger.write('%s\t%s\n' %(c.name, "STOPPED"))
           return (0,'')

        if verbose:
             eprintf("\n\n\nEXECUTED COMMAND : %s\n", ', '.join(c.commands) )

        eprintf("%s" %(c.message))

        if c.status in ['redo']:
            c.removeOutput(s)
            if c.isInputAvailable( errorlogger = s.errorlogger):
               s.stepslogger.write('%s\t%s\n' %(c.name, "RUNNING"))
               result = execute(s,c)
               if result[0] == 0 :
                  eprintf('..... Redo Success!\n')
                  s.stepslogger.write('%s\t%s\n' %( c.name, "SUCCESS"))
               else:
                  eprintf('..... Failed!\n')
                  s.stepslogger.write('%s\t%s\n' %( c.name, "FAILED"))
            else:
               eprintf('..... Skipping [NO INPUT]!\n')
               s.stepslogger.write('%s\t%s\n' %( c.name, "MISSING_INPUT"))

        elif c.status in ['yes']:
           if not c.isOutputAvailable():
               if c.isInputAvailable(errorlogger = s.errorlogger):
                  s.stepslogger.write('%s\t%s\n' %(c.name, "RUNNING"))
                  result = execute(s,c)
                  if result[0] == 0 :
                     eprintf('..... Success!\n')
                     s.stepslogger.write('%s\t%s\n' %( c.name, "SUCCESS"))
                  else:
                     eprintf('..... Failed!\n')
                     s.stepslogger.write('%s\t%s\n' %( c.name, "FAILED"))
               else:
                  eprintf('..... Skipping [NO INPUT]!\n')
                  s.stepslogger.write('%s\t%s\n' %(  c.name, "SKIPPED"))
           else:
               eprintf('..... Already Computed!\n')
               s.stepslogger.write('%s\t%s\n' %( c.name, "ALREADY_COMPUTED"))

        elif c.status in ['skip']:
           eprintf('..... Skipping!\n')
           s.stepslogger.write('%s\t%s\n' %(  c.name, "SKIPPED"))
        elif c.status=='grid':
           blastgrid(c.commands[0])



def execute(s, c):
       
       if len(c.commands) == 2:
             result = execute_pipeline_stage(c.commands[0], extra_command =  c.commands[1], errorlogger= s.errorlogger, runstatslogger = s.runstatslogger )
       else:
             result = execute_pipeline_stage(c.commands[0], errorlogger= s.errorlogger, runstatslogger = s.runstatslogger)
       return result





def print_to_stdout(s):
    print s
    
def no_status_updates(s):
    pass

def get_params_str(params):
    result = []
    for param_id, param_value in params.items():
        result.append('--%s' % (param_id))
        if param_value != None:
            result.append(param_value)
    return ' '.join(result)


## End  workflow and related functions
