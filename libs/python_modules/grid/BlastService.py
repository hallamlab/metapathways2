#!/usr/bin/python

from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = [""]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
   import subprocess, sys, re, inspect, time, shutil , traceback
   from os import makedirs, sys, listdir, environ, path
   from optparse import OptionParser
   from glob import glob

   from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
   from libs.python_modules.utils.metapathways_utils import printf, eprintf, fprintf, Job, Performance
   from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
   from libs.python_modules.utils.utils import number_of_lines_in_file
except:
   print """ Could not load some user defined  module functions"""
   print """ Make sure your typed \"source MetaPathwaysrc\""""
   print traceback.print_exc(10)
   print """ """
   sys.exit(3)

PATHDELIM = pathDelim()

class BlastService:
     SSH = 'ssh'
     SCP = ''
     PATHDELIM = pathDelim()
     WIN_RSA_KEY_ARGS = ["-i", 'executables' + PATHDELIM + 'win' + PATHDELIM + 'bit64' + PATHDELIM + 'win_rsa.ppk', '-batch']
     

     sub_string = ""
     submission_type = '0'


     def buildSSHLogin(self, connType='ssh'):  
         user, server = self.getUserServer()
         args = [ connType]

         if hasattr(self, 'keyfile'):
             args += ['-i', self.keyfile ]

         if connType=='ssh':
             args += [user+'@'+server]  

         return args


     def  isValid(self, opts):
         
         if not hasattr(opts,'sample_name'):
            return False
     
         if not hasattr(opts,'user') or not hasattr(opts,'server'):
            return False
     
         if opts.sample_name:
            return True 
         return False
     
     
     def getUserServer(self):
         return (self.user, self.server)


     def create_a_process(args):
        p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return p

     def _submit_remote_location(basicargs):

        args = self.buildSSHLogin()  
        args += ['python' ]
        args.extend(basicargs)
    
        print ' '.join(args)

        p = create_a_process(args)

        result = p.communicate()
        return result

     def _interpret_results(self, _result, expect):
         
          result = ' '.join(_result)
          if expect in result:
             return (True, result)
          else:
             return (False, result)
     
     
     def _parse_array_results(self, result, delim = ','):
          array = []
          for  x in result[0].split(delim):
             if x.strip():
               array.append(x.strip())
          return array

    
     def  __remote_check_if_server_is_up(self, user, server):
     #    user, server = self.getUserServer()

         args = self.buildSSHLogin()  

         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS

         args += ['echo','hello']
         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, 'hello')
         return boolean 
     
     def  __remote_Remove_Sample_Folders(self, samples,  home_dir ='~', working_dir = '~'):
         args = self.buildSSHLogin()  
         fullFolderNames = []
         for s in samples:
            fullFolderNames.append('--remove-sample-dirs')
            fullFolderNames.append(working_dir + '/MetaPathways/samples/' + s)
            
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS
         args = args + ['python', 'daemon.py','--home-dir', home_dir ] + fullFolderNames
         #print ' '.join(args)
         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, '<<Success!>>')
         return boolean

     def __local_No_of_Lines(self, filename):
         return number_of_lines_in_file(filename)

     def __remote_No_of_Lines(self, filename, home_dir='~', working_dir = '~' ):
         args = self.buildSSHLogin()  
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS
         args = args + ['python',  'daemon.py','--home-dir', home_dir, '--number-of-lines-in-file', working_dir + '/' + filename]
         #print ' '.join(args)
         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, 'number')
         try: 
            return int(message)
         except:
            return -1

     def __remote_deleteRemoteFile(self, filename, home_dir = '~',  working_dir='~'):
         args = self.buildSSHLogin()  
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS
         args = args + ['python',  'daemon.py','--home-dir', home_dir, '--remove-file', working_dir + '/' + filename]
         #print ' '.join(args)
         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, '<<Success!>>')
         return boolean


     def __remote_copy_file(self, source, target):
          user, server = self.getUserServer()
          args  = self.buildSSHLogin( 'scp')  
          args += [ source , user+'@'+server+':'+ target]

          p = self.create_a_process(args)
          result = p.communicate()
          return result[0]==''

     def __remote_copy_file_back(self, source, target):
          user, server = self.getUserServer()
          args  = self.buildSSHLogin( 'scp')  
          args += [user+'@'+server+':'+ source, target]

          p = self.create_a_process(args)
          result = p.communicate()
          (boolean, message)  = self._interpret_results(result, '')
          return boolean

     def __remote_DownloadFile(self, source, target):
         self.__remote_copy_file_back(source, target) 
         if path.exists(target) :
            return True
         else:
            return False

     def __remote_createFile(self, file_name, home_dir='~', working_dir = '~'):
         self.__remote_copy_file(self.MetaPathwaysDir  + PATHDELIM + self.Files[file_name][0] + PATHDELIM + file_name, working_dir + '/' + self.Files[file_name][1] + '/' + file_name) 
         if self.__remote_doesFileExist(working_dir + '/' + self.Files[file_name][1] + '/' + file_name) :
            return True
         else:
            return False

     def __remote_createFolder(self, folder_name, home_dir='~', working_dir = '~'):
         args  = self.buildSSHLogin()  
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS
         args = args + ['python', 'daemon.py','--home-dir', home_dir, '--create-sample-dir', working_dir + '/' + folder_name]
         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, '<<Success!>>')
         return boolean
     
     
     def __remote_doesFolderExist(self, folder_name, home_dir='~', working_dir = '~'):
         args  = self.buildSSHLogin()  
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS
         args += ['python', 'daemon.py','--home-dir', home_dir, '--does-sample-dir-exist', working_dir + '/' + folder_name]
         #print ' '.join(args)
         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, '<<Success!>>')
         #print  folder_name + ' folder ' + str(boolean)
         return boolean

     def __remote_doesFileExist(self, file_name, home_dir='~', working_dir = '~'):
         args  = self.buildSSHLogin()  
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS
         args += ['python', 'daemon.py','--home-dir', home_dir, '--does-file-exist', working_dir + '/' + file_name]
         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, '<<Success!>>')
         return boolean

     def __remote_doesFilePatternExist(self, file_pattern, home_dir='~', working_dir = '~'):
         args  = self.buildSSHLogin()  
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS
         args += ['python', 'daemon.py','--home-dir', home_dir, '--does-file-pattern-exist', working_dir + '/' + file_pattern]
         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, '<<Success!>>')
         return boolean

     def safeDoubleQuotes(self, string) :
         return ('\"' + string + '\"')
     def __remote_getFileNamesWithPattern(self, file_pattern, home_dir='~', working_dir = '~'):
         args  = self.buildSSHLogin()  
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS
         args += ['python', 'daemon.py','--home-dir', home_dir, '--get-files-with-pattern',  self.safeDoubleQuotes(working_dir + '/' + file_pattern)]
         #print ' '.join(args) 
         p = self.create_a_process(args)  
         result = p.communicate()
         array =  self._parse_array_results(result)
         #print array
         return array

     def __remote_isDBFormatted(self, remote_DB_dir, dbname, algorithm, working_dir = '~'):
            #Make sure DATABASESES are formatted
            dbnamePath = remote_DB_dir + '/' + dbname
            if algorithm.upper() == 'LAST': 
               suffixes = [ 'des', 'sds', 'suf', 'bck', 'prj', 'ssp', 'tis' ]
     
            if algorithm.upper() == 'BLAST': 
               suffixes = [ 'psq', 'phr', 'pin' ]
     
            for suffix in suffixes:
               filepattern = dbnamePath + '*' + '.' + suffix
                
               if not self.__remote_doesFilePatternExist(filepattern, working_dir = working_dir ):
                  return False
               
            return True
            
     #format the remote db
     def __remote_formatDB(self, dbname, algorithm, home_dir='~', working_dir= '~'):
         args  = self.buildSSHLogin()  
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS
         args += ['python', 'daemon.py','--home-dir', working_dir, '--format-database', dbname, '--algorithm', algorithm]
         #print ' '.join(args)
 
         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, '<<Success!>>')
         #print 'formatDB ' + str(boolean)
         return boolean

     
     def printError(self, result):
         if len(result[1].strip())>0:
           print "Remote Execution Error:"
           print "<-----------------------"
           for line in result:
              print "    " + line
           print "----------------------->"
        
     def printResponse(self, result):
         if len(result[0].strip())>0:
           print "Remote Execution Response:"
           print "<======================="
           for line in result:
              print "    " + line
           print "=======================>"
        
     def copy_file(self, source, target):
         user, server = self.getUserServer()
         if PATHDELIM=='\\':
             args = [self.SCP] + WIN_RSA_KEY_ARGS + [source, user+'@'+server+':'+ target]
         else:
             args  = self.buildSSHLogin('scp')  
             args += [source , user+'@'+server+':'+ target]
             
         p = self.create_a_process(args)  
         result = p.communicate()[0].strip()
         if result=='':
            return True
         else:
            return False
          
     def copy_file_back(self, source, target):
         user, server = self.getUserServer()
         if PATHDELIM=='\\':
            args = [self.SCP] + WIN_RSA_KEY_ARGS + [user+'@'+server+':'+ source, target]
         else: 
            args  = self.buildSSHLogin('scp')  
            args += [user+'@'+server+':~/'+ source, target]
     
         #print ' '.join(args)
         p = self.create_a_process(args)  
         result = p.communicate()[0].strip()
         if result=='':
            return True
         else:
            return False
     
     def __remote_Submit(self, J, working_dir='~'):
         args  = self.buildSSHLogin()  
         if PATHDELIM =='\\':
             args = args + WIN_RSA_KEY_ARGS

         args += ['python', 'daemon.py','--home-dir', working_dir, '--submit-job',\
                J.a, '--sample-name' , J.S,  '--algorithm', J.m, '--dbname', J.d,\
                '--submit-string', "\""+self.sub_string +"\"", '--submission-type', self.submission_type]
         command = ' '.join(args)

         p = self.create_a_process(args)  
         result = p.communicate()
         (boolean, message)  = self._interpret_results(result, '<<Success!>>')

         if boolean == False and self.submission_type=='0':
            self.submission_type='1' 

         return boolean, str(result), command

     
     def get_number_of_running_jobs(self, sample_name, algorithm):
         args  = self.buildSSHLogin()  
         args += ['python', 'MetaPathways/' + sample_name + '/daemon.py','--home-dir', '\'\'', '--get-number-of-running-jobs', user]
         p = self.create_a_process(args)  
         result = p.communicate()
         if result[1].strip()=='':
            try:
               return  int(result[0].strip())
            except:
               return 0
         else:
            return 0
     
     def databaseFiles(self, sourcedir, regPattern):
        files = [ re.sub(r'.*\/','',f) for f in glob(sourcedir + PATHDELIM + '*')  if regPattern.search(f) ] 
        return files
     
     #def retrieve_results(sample_name, dbname):
     
     def create_a_process(self, args):
         p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         return p
     
     
     class empty(object):
         pass
     
     def blastgrid(argv):
     
         opts = empty()
         for key, value in argv.items():
             setattr(opts, key, value)
               
     
         for dbname in opts.dbnames: 
            target = opts.sample_name +  PATHDELIM + 'blast_results' + PATHDELIM 
            source = MetaPathways + sample_name + '/' + sample_name + '.' + db_output_suffix
     
            copy_file_back(source, target)
     
         #print "Completed calculation!"
     
             

     if PATHDELIM=='/':
       SSH = 'ssh'
       SCP = 'scp'
     else:
       SSH = 'executables' + PATHDELIM + 'win' + PATHDELIM + 'bit64' + PATHDELIM + 'plink'
       SCP = 'executables' + PATHDELIM + 'win' + PATHDELIM + 'bit64' + PATHDELIM + 'pscp'
     
     def  __init__(self, gridParams):
         self.user = None
         self.service = None
         self.walltime= None
         self.remoteActiveFiles  = {}
         self.remoteActiveFolders= {}
         self.remoteFormattedDBs = {}
         self.serverUp = False 
         self.isaws=False
         self.awsparams={}
         self.remote_home_dir = '~'
         self.working_directory='~'
         self.base_output_folder= None
         self.max_parallel_jobs = 4
         self.os = 'mac'
         self.bits = 'bit64'
         self.type ='ordinary'
         self.sub_string = ""

         self.performance = Performance()
         self.A =1
         self.B =1
         self.C =1
         self.D =1

         try:
           for key, value in gridParams.iteritems():
              setattr(self, key, value) 
         except:
            print "ERROR : in  creating BlastService"
            pass

            
         if not self.user:
              self.messagelogger.write("ERROR: User for Grid service not specified\n")
              sys.exit(0)

         if not self.server:
              self.messagelogger.write("ERROR: Service for Grid service not specified\n")
              sys.exit(0)

         self.Folders = [ 'MetaPathways', 'MetaPathways/databases',  'MetaPathways/executables',  'MetaPathways/samples', \
                  'MetaPathways/samples/.qstatdir']
         self.Files = {  
                 'blastp': [ 'executables' + PATHDELIM + self.os + PATHDELIM + self.bit, 'MetaPathways/executables/'],
                 'blastn': [ 'executables' + PATHDELIM + self.os + PATHDELIM + self.bit, 'MetaPathways/executables/'],\
                 'lastal': [ 'executables' + PATHDELIM + self.os + PATHDELIM + self.bit, 'MetaPathways/executables/'],\
                 'lastdb': [ 'executables' + PATHDELIM + self.os + PATHDELIM + self.bit, 'MetaPathways/executables/'],\
                 'makeblastdb': [ 'executables' + PATHDELIM + self.os + PATHDELIM + self.bit, 'MetaPathways/executables/'],\
                 'daemon.py': [ 'libs' + PATHDELIM + 'python_scripts', '']
              }


     def isUp(self):
         user, server = self.user, self.server
         if  self.__remote_check_if_server_is_up(user , server ):
            return True
         else:
            return False

     def set_sample_name(self, sample_name):
         self.sample_name = sample_name

     def set_base_output_folder(self, base_output_folder):
         self.base_output_folder = base_output_folder

     def submitJob(self, J):
         if not self.isInputReady(J):
            self.messagelogger.write("ERROR: Input not ready (%s %s %s %s)" %(J.S, J.d, J.a, J.m, self.server))
            return False
         status = self.submit(J)
         if not status[0]:
            print "failed to submit " + "%s %s %s %s %s" %(J.S, J.d, J.a, J.m, self.server)
            print "ERROR: Remote error message " + status[1]
            print "INFO: command : " + status[2]
            return False

         #print "submitted " + "%s %s %s %s %s" %(J.S, J.d, J.a, J.m, self.server)
         return True

     def submit(self, J):
         return self.__remote_Submit(J, working_dir = self.working_directory)

     def areNewResultsAvailable(self, sample_names, algorithm):
          remote_Samples_dir  ='MetaPathways/samples/' 
          for sample_name in sample_names:
               filepattern = remote_Samples_dir + '/' + sample_name + '/' +  '*' + '.' + algorithm.upper() 
               if self.__remote_doesFilePatternExist(filepattern, working_dir = self.working_directory ):
                 return True
          return False

     def isInputReady(self, J):
         sampleDir = 'MetaPathways/samples/' + J.S 
         Folders = [ sampleDir, sampleDir + '/.qstatdir' ]
         # check if the remote folders are active
         for target in Folders:
           if not target in self.remoteActiveFolders:
             if not self.__remote_doesFolderExist(target, working_dir = self.working_directory):
                 self.__remote_createFolder(target, working_dir = self.working_directory)
                 if  self.__remote_doesFolderExist(target, working_dir = self.working_directory):
                    self.remoteActiveFolders[target] = True
                    self.messagelogger.write("SUCCESS: Successfully created remote folder \"%s\" in Server \"%s\"!\n" %(target, self.server))
                 else:
                    self.messagelogger.write("ERROR: Cannot create remote folder \"%s\" in Server \"%s\"!\n" %(target, self.server))
                    return False
             else:
                 self.messagelogger.write("OK: Found already created remote folder \"%s\" in Server \"%s\"!\n" %(target, self.server))
                 self.remoteActiveFolders[target] = True

         source = self.base_output_folder +  PATHDELIM + J.S + PATHDELIM + 'blast_results' + PATHDELIM + 'grid' + PATHDELIM + 'split_batch' + PATHDELIM + J.a
         target  =  'MetaPathways/samples/' + J.S + '/' + J.a
         if not J.a in self.remoteActiveFiles and  not self.__remote_doesFileExist(target, working_dir = self.working_directory):
             self.__remote_copy_file(source, self.working_directory + '/' + target)
             if not self.__remote_doesFileExist(target, working_dir = self.working_directory):
                self.messagelogger.write("ERROR: Cannot upload query file \"%s\" to server \"%s\"!\n" %(J.a, self.server))
                return False
             else:
                self.remoteActiveFiles[J.a]=True

         source = self.blast_db_folder +  PATHDELIM + J.d
         target  =  'MetaPathways/databases/'  + J.d
        
         if not J.d in self.remoteActiveFiles and  not self.__remote_doesFileExist(target, working_dir = self.working_directory):
             if not path.exists(source):
                self.messagelogger.write("ERROR: DB file \"%s\" for DB \"%s\" is missing, required to upload to server \"%s\"!\n" %(source, J.d, self.server))
                return False

             self.messagelogger.write("STATUS: Uploading  DB file \"%s\" to server \"%s\"!\n" %(J.d, self.server))
             self.__remote_copy_file(source, self.working_directory + '/' + target)
             if not self.__remote_doesFileExist(target, working_dir = self.working_directory):
                self.messagelogger.write("ERROR: Cannot upload DB file \"%s\" to server \"%s\"!\n" %(J.d, self.server))
                return False
             else:
                self.messagelogger.write("OK: Founde already Uploaded  DB file \"%s\" to server \"%s\"!\n" %(J.d, self.server))
                self.remoteActiveFiles[J.d]=True


         remote_DB_dir  ='MetaPathways/databases/' 
         if not J.d in self.remoteFormattedDBs: 
           
           if  not self.__remote_isDBFormatted(remote_DB_dir, J.d, J.m, working_dir = self.working_directory ):
             self.messagelogger.write("WARNING: DB \"%s\" is not formatted on server \"%s\" for algorithm \"%s\"!\n" %(J.d, self.server, J.m))
             if self.__remote_formatDB(remote_DB_dir +'/' + J.d, J.m, working_dir = self.working_directory)  and  self.__remote_isDBFormatted(remote_DB_dir, J.d, J.m, working_dir = self.working_directory ):
                self.messagelogger.write("SUCCESS: Successfully formatted  DB \"%s\" on server \"%s\" for algorithm \"%s\"!\n" %(J.d, self.server, J.m))
                self.remoteFormattedDBs[J.d]=True
             else:
                self.messagelogger.write("ERROR: Failed to format  DB \"%s\"  on server \"%s\" for algorithm \"%s\"!\n" %(J.d, self.server, J.m))
                return False
           else:
             self.messagelogger.write("OK: DB \"%s\" is already formatted on server \"%s\" for algorithm \"%s\"!\n" %(J.d, self.server, J.m))
             self.remoteFormattedDBs[J.d]=True
    
         return True
           

     def isReadyToSubmit(self, load):
         self.A =  max(self.A-1, 1)

         if self.A > 1:
            return False

         if not self.isSetUp():
            self.B = min(2*self.B, 100)
            self.A = self.B
            return False
         else:
            if self.isQFull(load):
               self.B = min(2*self.B, 100)
               self.A = self.B
               return False
            else:
               self.B = self.A
               return True
         return True
      
     def isReadyToHarvest(self, sample_names, algorithm):
         self.C =  max(self.C-1, 1)
         if self.C > 1:
            return False
         if not self.isSetUp():
            self.D = min(2*self.D, 100)
            self.C = self.D
            return False
         else:
            if not self.areNewResultsAvailable(sample_names, algorithm):
               self.D = min(2*self.D, 100)
               self.C = self.D
               return False
            else:
               self.D = self.C
               return True
         return True

     def splitResultName(self, _a_result):
         
         fields = [ x.strip() for x in _a_result.split('.') ]
         if len(fields) != 3:
            return (None, None,  None)
         return (fields[0], fields[1], fields[2])

     def  job_Was_Submitted_By_Current_Server(self, sample_name, split_name,  dbname, algorithm, list_jobs_submitted) :
         if not sample_name in list_jobs_submitted:
            return False 

         if not dbname in list_jobs_submitted[sample_name]:
            return False 

         if not split_name in list_jobs_submitted[sample_name][dbname]:
            return False 

         if not algorithm in list_jobs_submitted[sample_name][dbname][split_name]:
            return False 

         if not self.server in list_jobs_submitted[sample_name][dbname][split_name][algorithm]:
            return False 

         return True


     def harvest(self, sample_names, algorithm, list_jobs_submitted):
          _all_results = []
          remote_Samples_dir  ='MetaPathways/samples/' 
          for sample_name in sample_names:
              filepattern = remote_Samples_dir + '/' + sample_name + '/' +  '__DELIMITER__' +  algorithm 
              _results = self.__remote_getFileNamesWithPattern(filepattern, working_dir = self.working_directory )
              for _a_result in _results:
                 split_name, dbname, algorithm = self.splitResultName(_a_result)
                 if not  self.job_Was_Submitted_By_Current_Server(sample_name, split_name,  dbname, algorithm, list_jobs_submitted): 
                   continue
                 
                 source = self.working_directory + '/' + remote_Samples_dir + '/'  + sample_name + '/' + _a_result
                 localTargetFile = self.base_output_folder+ PATHDELIM + sample_name + PATHDELIM +  'blast_results' + PATHDELIM + 'grid' + PATHDELIM + 'split_results' + PATHDELIM  + _a_result
                 if self.__remote_DownloadFile(source, localTargetFile):
                    nL = self.__local_No_of_Lines(localTargetFile)
                    remoteFile = remote_Samples_dir + '/' + sample_name + '/' + _a_result
                    nR = self.__remote_No_of_Lines(remoteFile, working_dir=self.working_directory)
                    #print "Lines check ", nL, nR
                    if nR == nL:
                        #remoteFile = remote_Samples_dir + '/' + sample_name + '/' + '_a_result'
                        if self.__remote_deleteRemoteFile(remoteFile, working_dir = self.working_directory):
                           J = Job(sample_name, dbname, split_name, algorithm)
                           _all_results.append(J) 
                 else:
                     self.messagelogger.write("ERROR: Failed to download file \"%s\"  from server \"%s\"!\n" %(source,  self.server))

          return _all_results 
         
     def isQFull(self, size):
         if size >= int(self.max_parallel_jobs):
           return True
         else:
           return False



     def isSetUp(self):
         if  not self.isUp():
            self.messagelogger.write("ERROR: Cannot connect to  Server \"%s\"! Perhaps it is down!\n" %(self.server))
            return False

        # copy the daemon.py file
         daemonFile = self.MetaPathwaysDir  + PATHDELIM + 'libs' + PATHDELIM + 'python_scripts' +  PATHDELIM + 'daemon.py'
         daemonFileRemote =  'daemon.py'
         #print  self.server + ' ' +  daemonFile + '   ' +  self.remote_home_dir + '/' + daemonFileRemote
         #print self.remoteActiveFiles
         if not daemonFileRemote in self.remoteActiveFiles:
              #if not self.__remote_doesFileExist(daemonFileRemote):
              if self.__remote_copy_file(daemonFile, self.remote_home_dir + '/' + daemonFileRemote) :
                  self.messagelogger.write("SUCCESS: Successfully uploaded the daemon file to  in Server \"%s\"!\n" %(self.server))
              else:
                  self.messagelogger.write("ERROR: Could not upload daemon file to  in Server \"%s\"!\n" %(self.server))
                  return False

              self.remoteActiveFiles[daemonFileRemote] = True

        # create the folders under the working directory
         for  f in self.Folders:   
           if not f in self.remoteActiveFolders:
              self.messagelogger.write("MESSAGE: Checking remote folder \"%s\" in Server \"%s\"!\n" %(f, self.server))
              if not self.__remote_doesFolderExist(f, working_dir = self.working_directory):
                   if self.__remote_createFolder(f, working_dir = self.working_directory):
                      self.remoteActiveFolders[f] = True
                      self.messagelogger.write("SUCCESS: Successfully created remote folder \"%s\" in Server \"%s\"!\n" %(f, self.server))
                   else:
                      self.messagelogger.write("ERROR: Cannot create remote folder \"%s\" in Server \"%s\"!\n" %(f, self.server))
                      return False
              else:
                 self.messagelogger.write("STATUS: Found existing remote folder \"%s\" in Server \"%s\"!\n" %(f, self.server))
                 self.remoteActiveFolders[f] = True


         for  f in self.Files:   
           if not f in self.remoteActiveFiles: 
             if not self.__remote_doesFileExist(self.Files[f][1] + '/' + f, working_dir = self.working_directory):
                  #print self.working_directory + '/' +   self.Files[f][1] + '/' + f
                  self.messagelogger.write("STATUS: Creating/Uploading remote file \"%s\" in Server \"%s\"!\n" %(f, self.server))
                  if self.__remote_createFile(f, working_dir = self.working_directory ):
                     self.remoteActiveFiles[f] = True
                     self.messagelogger.write("SUCCESS: Successfully created/uploaded remote file \"%s\" in Server \"%s\"!\n" %(f, self.server))
                  else:
                     self.messagelogger.write("ERROR: Cannot create/upload remote file \"%s\" in Server \"%s\"!\n" %(f, self.server))
                     return False
             else:
                self.remoteActiveFiles[f] = True
    
         return True


     def deleteRemoteSampleFolders(self, samples):
         if  self.isSetUp():
            if self.__remote_Remove_Sample_Folders(samples, working_dir = self.working_directory):
                self.messagelogger.write("SUCCESS: Successfully removed sample folders in Server \"%s\"!\n" %(self.server))
            else:
                self.messagelogger.write("WARNING: Faild to remove sample folders in Server \"%s\"!\n" %(self.server))
                self.messagelogger.write("WARNING: Please DO NOT forget to remove the stale folders manually in Server \"%s\"!\n" %(self.server))
         
