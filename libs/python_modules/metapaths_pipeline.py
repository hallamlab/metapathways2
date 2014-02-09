#!/usr/bin/env python
# File created on 30 Dec 2009.


from __future__ import division
try:
      from subprocess import Popen, PIPE, STDOUT
      from os import makedirs, listdir
      from glob import glob
      from os.path import split, splitext, join, dirname, abspath
      from datetime import datetime
      from metapaths_utils import printf, eprintf
      from sysutil import getstatusoutput
      import sys, traceback
      from BlastGrid import *
      from optparse import OptionParser
      import python_scripts 

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

## Start utilities used by the pipeline functions
def generate_log_fp(output_dir,
                    basefile_name='metapathways_run_log',
                    suffix='txt',
                    timestamp_pattern=''):
    timestamp = datetime.now().strftime(timestamp_pattern)
    filename = '%s.%s' % (basefile_name,suffix)
    return join(output_dir,filename)

def generate_steps_log_fp(output_dir,
                    basefile_name='metapathways_steps_log',
                    suffix='txt'):
    filename = '%s.%s' % (basefile_name,suffix)

    return join(output_dir,filename)

class WorkflowError(Exception):
    pass


def contract_key_value_file(fileName):

     file = open(fileName,'r')
     lines = file.readlines()
     if len(lines) < 20:
        file.close()
        return

     keyValuePairs = {}
     
     for line in lines:
       fields = [ x.strip() for x in line.split('\t') ] 
       if len(fields) == 2:
          keyValuePairs[fields[0]] = fields[1]
     file.close()

     file = open(fileName,'w')
     for key, value in  keyValuePairs.iteritems():
          fprintf(file, "%s\t%s\n",key, value)
     file.close()

     
class WorkflowLogger(object):
    
    def __init__(self,log_fp=None,params=None,metapaths_config=None,open_mode='w'):
        if log_fp:

        #contract the file if we have to
            if open_mode=='c':
                try:
                   contract_key_value_file(log_fp)
                except:
                   pass 
                open_mode='a'
            self._f = open(log_fp,open_mode)
        else:
            self._f = None
        self._filename = log_fp
        #start_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.writemetapathsConfig(metapaths_config)
        self.writeParams(params)

    def get_log_filename(self): 
        return self._filename

    def write(self,s):
        if self._f:
            self._f.write(s)
            # Flush here so users can see what step they're
            # on after each write, since some steps can take
            # a long time, and a relatively small amount of 
            # data is being written to the log files.
            self._f.flush()
        else:
            pass
    
    def writemetapathsConfig(self,metapaths_config):
        if metapaths_config == None:
            #self.write('#No metapaths config provided.\n')
            pass
        else:
            self.write('#metapaths_config values:\n')
            for k,v in metapaths_config.items():
                if v:
                    self.write('%s\t%s\n' % (k,v))
            self.write('\n')
            
    def writeParams(self,params):
        if params == None:
            #self.write('#No params provided.\n')
            pass 
        else:
            self.write('#parameter file values:\n')
            for k,v in params.items():
                for inner_k,inner_v in v.items():
                    val = inner_v or 'True'
                    self.write('%s:%s\t%s\n' % (k,inner_k,val))
            self.write('\n')
    
    def close(self):
        end_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('\nLogging stopped at %s\n' % end_time)
        if self._f:
            self._f.close()
        else:
            pass

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
            
def execute_pipeline_stage(pipeline_command):
     argv = [ x.strip() for x in pipeline_command.split() ]
     funcname = re.sub(r'.py$','', argv[0])
     funcname = re.sub(r'^.*/','', funcname)
     args = argv[1:] 
     
     if hasattr(python_scripts, funcname):
        methodtocall = getattr( getattr(python_scripts, funcname), funcname)
        result = methodtocall(args)
     else:
        result = getstatusoutput(pipeline_command)
     return result



def call_commands_serially(commands, status_update_callback, logger, stepslogger, params ):
    """Run list of commands, one after another """
    #logger.write("Executing commands.\n\n")
    from sys import exit
    
    make = False
    for c in commands:
        if c[3]=='stop':
           print "Stopping!"
           stepslogger.write('%s\t%s\n' %(c[2], "STOPPED"))
           sys.exit(0)
           return (0,'')

        if params['verbose']:
             eprintf("\n\n\nIssuing Command : %s\n", c[1])

        eprintf("%s" %(c[0]))
        stepslogger.write('%s\t%s\n' %(c[2], "RUNNING"))
        if c[3]=='missing':
           print "..... Input Missing!"
           stepslogger.write('%s\t%s\n' %(c[2], "INPUT_MISSING"))
           return (0,'')


        if c[3] in ['yes', 'redo' ] and c[4]:
           #result = getstatusoutput(c[1])
           result = execute_pipeline_stage(c[1])
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

        logger.write('COMMAND : \n%s \n' %( c[1]))
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
