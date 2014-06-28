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
           result = methodtocall(args, errorlogger = errorlogger, extra_command = extra_command, runstatslogger = runstatslogger)
     else:
        result = getstatusoutput(pipeline_command)
     return result



def call_commands_serially(commands, status_update_callback, logger, stepslogger, errorlogger, runstatslogger,  params):
    """Run list of commands, one after another """
    #logger.write("Executing commands.\n\n")
    from sys import exit
    
    make = False
    for c in commands:
        if c[3]=='stop':
           print "Stopping!"
           stepslogger.write('%s\t%s\n' %(c[2], "STOPPED"))
           return (0,'')

        if params['verbose']:
             eprintf("\n\n\nEXECUTED COMMAND : %s\n", c[1])

        eprintf("%s" %(c[0]))
        stepslogger.write('%s\t%s\n' %(c[2], "RUNNING"))
        if c[3]=='missing':
           print "..... Input Missing!"
           stepslogger.write('%s\t%s\n' %(c[2], "INPUT_MISSING"))
           return (0,'')


        if c[3] in ['yes', 'redo' ] and c[4]:
           #result = getstatusoutput(c[1])
           if type(c[1]) is tuple:
              result = execute_pipeline_stage(c[1][0], extra_command =  c[1][1], errorlogger= errorlogger, runstatslogger = runstatslogger )
           else:
              result = execute_pipeline_stage(c[1], errorlogger= errorlogger, runstatslogger = runstatslogger)

           if result[0] == 0 :
             if c[3] in ['redo']:
                eprintf('..... Redo Success!\n')
                stepslogger.write('%s\t%s\n' %(c[2], "SUCCESS"))
             else:
                eprintf('..... Success!\n')
                stepslogger.write('%s\t%s\n' %(c[2], "SUCCESS"))
           else:
             eprintf('..... Failed!\n')
             stepslogger.write('%s\t%s\n' %(c[2], "FAILED"))
             #print c
        elif c[3]=='grid' and c[4]:
           blastgrid(c[1])
           continue  
        elif c[3]=='skip':
           eprintf('..... Skipping!\n')
           stepslogger.write('%s\t%s\n' %(c[2], "SKIPPED"))
           continue
        elif c[3]=='missing':
           eprintf('..... Input Missing!\n')
           stepslogger.write('%s\t%s\n' %(c[2], "INPUT_MISSING"))
           continue
        else:
           eprintf('..... Already Computed!\n')
           stepslogger.write('%s\t%s\n' %(c[2], "ALREADY_COMPUTED"))
           continue

        #status_update_callback('%s\n%s' %(result[0], result[0]))

        logger.write('COMMAND : \n%s \n' %( ' '.join(c[1]))  )
        if result[0] == 0 :
            logger.write('Success!:\n%s\n' %( result[1]))
        else:
            #print result[1]
            try:
               printf('Error! : %s\n' %(result[1]))
               logger.write('Error!:\n%s\n' %( result[1]))
            except:
               pass

            break 

      #  timestamp_pattern='%Y-%m-%d %H:%M:%S'
      #  timestamp = datetime.now().strftime(timestamp_pattern)


def execute_tasks(s, verbose = False):
    """Run list of commands, one after another """
    #logger.write("Executing commands.\n\n")

    for c in s.getContexts():

        #print c.name, c.status, 'status'
        if c.status=='stop':
           print "Stopping!"
           s.stepslogger.write('%s\t%s\n' %(c.name, "STOPPED"))
           return (0,'')

        if verbose:
             eprintf("\n\n\nEXECUTED COMMAND : %s\n", ', '.join(c.commands) )

        eprintf("%s" %(c.message))
        s.stepslogger.write('%s\t%s\n' %(c.name, "RUNNING"))

        if c.status in ['redo']:
            c.removeOutput(s)
            if c.isInputAvailable( errorlogger = s.errorlogger):
               result = execute(s,c)
               if result[0] == 0 :
                  eprintf('..... Redo Success!\n')
                  s.stepslogger.write('%s\t%s\n' %( ','.join(c.name), "SUCCESS"))
               else:
                  eprintf('..... Failed!\n')
                  s.stepslogger.write('%s\t%s\n' %(','.join(c.name), "FAILED"))
            else:
               eprintf('..... Skipping [NO INPUT]!\n')
               s.stepslogger.write('%s\t%s\n' %( c.name, "MISSING_INPUT"))

        elif c.status in ['yes']:
           if not c.isOutputAvailable():
               if c.isInputAvailable(errorlogger = s.errorlogger):
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
