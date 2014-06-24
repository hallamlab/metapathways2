#!/usr/bin/python

from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = [""]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


try:
    import subprocess
    import sys, re, inspect, time, shutil 
    from os import makedirs, sys, listdir, environ, path, remove
    from optparse import OptionParser
    from glob import glob


    from libs.python_modules.utils.metapathways_utils import printf, eprintf, fprintf, Job, Performance
    from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
    from libs.python_modules.grid.BlastService import BlastService
    from libs.python_modules.utils.utils import *
    #import curses
except:
    print """ Could not load some user defined  module functions"""
    print """ Make sure your typed \"source MetaPathwaysrc\""""
    print """ """
    traceback.print_exc(10)
    sys.exit(3)

PATHDELIM = pathDelim()
WIN_RSA_KEY_ARGS = ["-i", 'executables' + PATHDELIM + 'win' + PATHDELIM + 'bit64' + PATHDELIM + 'win_rsa.ppk', '-batch']

class BlastBroker:
      sample_name = None
      fastaFile =None
      database =None
      blastType =None
      batchSize =None
      working_directory =None
      gridParams=None
      x = 0
      y = 0

#      base_output_folder = None
#      BlastServices=[]
#      BlastServers={}

#      samples_and_inputs = {}
#      samples_and_databases = {}
#      samples_and_algorithms= {}
#      batchSize =500
#
#
#      list_splits = {}
#      list_jobs = {}
#      list_jobs_submitted = {}
#      list_jobs_completed = {}
#      serverLoads = {}
#      list_jobs_submitted_file = 'list_jobs_submitted.txt'

      def __init__(self, messagelogger):
          self.messagelogger= messagelogger
          try:
               self.base_output_folder = None
               self.last_server_submitted_to = None
               self.BlastServices=[]
               self.BlastServers={}
               self.samples_and_inputs = {}
               self.samples_and_databases = {}
               self.samples_and_algorithms= {}
               self.batchSize =500
     
               self.list_splits = {}
               self.list_jobs = {}
               self.list_jobs_submitted = {}
               self.list_jobs_completed = {}
               self.serverLoads = {}
               self.performance = Performance()
               self.list_jobs_submitted_file = 'list_jobs_submitted.txt'

               self.pending_job = None
               self.displaylen = 0

               self.migrationMatrix = {}

               self.messagelogger.write("OK: Blast Broker initialized!\n")
          except:
               self.messagelogger.write("ERROR: Failed to initialize Blast Broker!\n")
               sys.exit(0)
     
          return None
      def getAverageDelay(self, server = None):
           return self.performance.getAverageDelay(server)

      def getStdDeviationDelay(self, server = None):
           return self.performance.getStdDeviationDelay(server)

      def setLastSubmittedServerTo(self, server):
          self.last_server_submitted_to = server
          
      def getLastSubmittedServerTo(self):
          return self.last_server_submitted_to
          
           
      def setAlgorithm(self, algorithm, algorithmSuffix):
           self.algorithm=algorithm
           self.algorithmSuffix=algorithmSuffix
           return True

      def setBatchSize(self, batchSize):
         self.batchSize = batchSize
         return True

      def getBatchSize(self):
         return self.batchSize

      def setBaseOutputFolder(self, base_output_folder):
         self.base_output_folder = base_output_folder
         self.messagelogger.write("OK: Added base output folder \"%s\" \n"  %(self.base_output_folder))
         return True

      def addService(self, service):
          self.BlastServices.append(service) 
          self.BlastServers[service.server ]  = service
          self.serverLoads[service.server] = 0
          self.messagelogger.write("OK: Added blast service \"%s\" \n"  %(service.server))


          return True


      def getServices(self):
          return self.BlastServers.keys()

      def getService(self, servicename):
          if not servicename in self.BlastServers: 
             self.messagelogger.write("ERROR: No service \"%s\" registered yet!\n"  %(servicename))
          return self.BlastServers[servicename]

      def addSamples(self,  samples_and_inputs): 
          if not self.base_output_folder:
              self.messagelogger.write("ERROR: Base output folder should be set before adding sample and inputs!\n")
              return  False

          self.samples_and_inputs = samples_and_inputs
          return True

      def getAlgorithm(self, sample):
           if not sample in self.samples_and_algorithms:
              self.messagelogger.write("WARNING: Algorithm not added  to sample \"%s\"!\n" %(sample))
              return None
           return self.samples_and_algorithms[sample]

      #compute the loads of the active servers 
      def compute_server_loads(self):
          for sample in self.list_jobs_submitted:
             for db in self.list_jobs_submitted[sample]:
               for split in self.list_jobs_submitted[sample][db]:
                   for server in self.list_jobs_submitted[sample][db][split][self.algorithm]:
                       if not server in self.BlastServers: 
                          continue
                       self.serverLoads[server] += 1
                       if not sample in self.list_jobs_completed:
                          continue
                       if not db in self.list_jobs_completed[sample]:
                          continue
                       if not split in self.list_jobs_completed[sample][db]:
                          continue
                       if not server in self.list_jobs_completed[sample][db][split][self.algorithm]:
                          continue
                       self.serverLoads[server] -= 1

      def compute_performance(self):
          for sample in self.list_jobs_completed:
             for db in self.list_jobs_completed[sample]:
               for split in self.list_jobs_completed[sample][db]:
                   for self.algorithm in self.list_jobs_completed[sample][db][split]:
                     for server in self.list_jobs_completed[sample][db][split][self.algorithm]:
                       if not sample in self.list_jobs_submitted:
                          continue
                       if not db in self.list_jobs_submitted[sample]:
                          continue
                       if not split in self.list_jobs_submitted[sample][db]:
                          continue
                       if not self.algorithm in self.list_jobs_submitted[sample][db][split]:
                          continue
                       if not server in self.list_jobs_submitted[sample][db][split][self.algorithm]:
                          continue
                       #print  sample + ' ' + db + ' ' + split + ' ' + self.algorithm + ' ' + server + ' '  + str(self.list_jobs_completed[sample][db][split][self.algorithm][server])
                       #print   str(self.list_jobs_submitted[sample][db][split][self.algorithm][server]) + ' ' +  str(self.list_jobs_completed[sample][db][split][self.algorithm][server]) + ' ' +  str(self.list_jobs_submitted[sample][db][split][self.algorithm][server] - self.list_jobs_completed[sample][db][split][self.algorithm][server])
                       data =  self.list_jobs_completed[sample][db][split][self.algorithm][server] - self.list_jobs_submitted[sample][db][split][self.algorithm][server]
                       self.performance.addPerformanceData(server, data)

      def startAWS(self, server): 
         clustername = server.cluster_name
         master_name = cli.cli( ['listmaster', clustername] )
         if master_name:
              self.messagelogger.write("OK: Found an existing AWS cluster with  master %s!\n" %(master_name))
         else:
             self.messagelogger.write("WARNING: No AWS cluster is running currently!\n")
             try:
                 start_cluster = cli.cli( ['start', clustername, '-c',  server.amazon_aws_config  ] ) 
             except:
                 start_cluster = cli.cli( ['terminate', clustername ] ) 
                 start_cluster = cli.cli( ['start', clustername, '-c',  server.amazon_aws_config  ] ) 


             master_name = cli.cli( ['listmaster', clustername] )
             if master_name:
                 self.messagelogger.write("SUCCESS: Successfully created a cluster with master %s!\n" %(master_name))
             else:
                 self.messagelogger.write("FAILED: Failed to create a cluster!\n")
                 return None
#
         return master_name


      def addAlgorithm(self,  sample, algorithm): 
          if not self.samples_and_inputs:
              self.messagelogger.write("ERROR: Samples and intputs should be set before adding algorithm!\n")
              return  False

          self.samples_and_algorithms[sample]=algorithm
          self.algorithm = algorithm   # Fix me if you want to support many types of algorithm in one run
          self.messagelogger.write("OK: Added algorithm \"%s\" to sample \"%s\"!\n" %(algorithm, sample))
          return True
         
      def addDatabase(self,  sample, database): 
          if not self.samples_and_inputs:
              self.messagelogger.write("ERROR: Samples and intputs should be set before adding database!\n")
              return  False

          if not sample in self.samples_and_databases:
              self.samples_and_databases[sample]= {}
             
          self.samples_and_databases[sample][database] =True
          self.messagelogger.write("OK: Added database \"%s\" to sample \"%s\"!\n" %(database, sample))
          return True

      def getDBs(self, samples = [] ):
          dbs = []
          for sample in samples:
             if not sample in self.samples_and_databases:
                 self.messagelogger.write("ERROR: \"%s\" is not a sample in Broker to retrieve database !\n" %(sample))
             else:
                 dbs = dbs + self.samples_and_databases[sample].keys()
          return dbs


      def checkDBs(self):
           DBs = self.getDBs()
           print 'FIXME'
           for D in DBs:
               if not doesFileExist(D):
                   self.messagelogger.write("ERROR: Database \"%s\" not found!\n" %(D))
                   return False
           return True
                

      def getSamples(self):
           return self.samples_and_inputs.keys() 
         

      def get_working_folder_paths(self):
         FOLDERS = [ 'grid',  'grid' + PATHDELIM + 'split_batch',  'grid' + PATHDELIM + 'split_results']
         for servicename in self.getServices():
             FOLDERS.append( 'grid'  + PATHDELIM + servicename )

         FOLDERPATHS = [] 
         for sample in self.getSamples():
             parentDir =  self.base_output_folder + PATHDELIM + sample + PATHDELIM + 'blast_results'
             for F in FOLDERS:
                FOLDERPATHS.append(parentDir + PATHDELIM + F)

         return FOLDERPATHS

      def are_working_folders_available(self, eraseExisting = False):
         FOLDERPATHS =  self.get_working_folder_paths()
         status = True
         for f in FOLDERPATHS:
             if not path.exists(f):
               status = False
               break

         return status 
      
      def create_working_folders(self, eraseExisting = False):
         FOLDERPATHS =  self.get_working_folder_paths()
         for f in FOLDERPATHS:
             if eraseExisting:
                removeFolderIfFound(f)
             createFolderIfNotFound(f)

         return True 
      
      def createJobs(self, redo = False):
          for S in self.getSamples():
             if redo:
                self.removeJobListIfExists(S)
             self.load_job_list_file(S)
             jobListStatus = self.isValidJobList(S)

             if jobListStatus == 1:
                self.messagelogger.write("OK: Previously created job list is correct for sample \"%s\"!\n" %(S))
                continue
             else:
                # job list not found or corrupt
                if jobListStatus == 2:
                    self.messagelogger.write("WARNING: Previously created job list is incorrect for sample \"%s\"!\n" %(S))
                if self.createJobList(S, fromScratch=False):
                   self.messagelogger.write("SUCCESS: Created the job list for sample \"%s\"!\n" %(S))

          return True

      def createJobList(self, S, fromScratch=False):
          parentDir =  self.base_output_folder + PATHDELIM + S + PATHDELIM + 'blast_results' + PATHDELIM + 'grid'
          job_list_file =  parentDir + PATHDELIM + 'list_jobs.txt'

          if fromScratch:
             removeIfFileExists( job_list_file)
          try:
             if not doesFileExist(job_list_file):
                 listfile  = open(job_list_file, 'w')
                 listfile.close()
          except:
                 self.messagelogger.write("ERROR: Cannot open job list file for sample \"%s\"!\n" %(S))
                 return False

          self.load_job_list_file(S)
          DBs = self.getDBs(samples=[S])
          self.load_list_splits_for_sample(S)

          if doesFileExist(job_list_file):
              listfile  = open(job_list_file, 'a')
              for d in DBs:
                 for a in self.list_splits[S]:
                    J = Job(S, d, a, self.getAlgorithm(S) )  
                    if not self.isJobInList(J):
                       self.add_job_to_list_jobs(J, listfile)
              listfile.close()

          return True

      def add_job_to_list_jobs(self, J, listfile):
          fprintf(listfile, "%s\t%s\t%s\t%s\n" %(J.S, J.d, J.a, self.getAlgorithm(J.S)) )
          return True


      def load_job_list_file(self, S):
           
           parentDir =  self.base_output_folder + PATHDELIM + S + PATHDELIM + 'blast_results' + PATHDELIM + 'grid'
           list_jobs_file=parentDir + PATHDELIM + 'list_jobs.txt'
            

           # return if the list_jobs.txt file does not exist
           if not path.exists(list_jobs_file):
               return

           if enforce_number_of_fields_per_row(list_jobs_file, 4):
               self.messagelogger.write('WARNING: Corrupt list file for sample %s' %(S))
               self.messagelogger.write('STATUS: Successfully recovered corrupt list file for sample %s' %(S))

           # list_jobs splitname, 
           if path.exists(list_jobs_file):
                 listfile = open(list_jobs_file, 'r')
                 lines = listfile.readlines()
                 listfile.close()
                 for line in lines:
                     fields = [ x.strip() for x in line.strip().split('\t') ]
                     if len(fields) == 4:
                          if not fields[0] in  self.list_jobs:
                             self.list_jobs[fields[0]] = {}

                          if not fields[1] in  self.list_jobs[fields[0]]:
                             self.list_jobs[fields[0]][fields[1]] = {}

                          if not fields[2] in  self.list_jobs[fields[0]][fields[1]]:
                             self.list_jobs[fields[0]][fields[1]][fields[2]] = {}

                          self.list_jobs[fields[0]][fields[1]][fields[2]][fields[3]] = True

      def load_job_status_lists(self):
           status = True
           for S in self.getSamples():
              if self.load_job_status_list(S):
                  self.messagelogger.write("OK: Loaded jobs submitted and complete list for sample \"%s\"!\n" %(S))
              else:
                  self.messagelogger.write("ERROR: Failed to load jobs submitted and complete list for sample \"%s\"!\n" %(S))
                  status = False
              
           return status 
         


      def load_job_status_list(self, S):
           parentDir =  self.base_output_folder + PATHDELIM + S + PATHDELIM + 'blast_results' + PATHDELIM + 'grid'
           list_jobs_submitted_file=parentDir + PATHDELIM + 'list_jobs_submitted.txt'
           list_jobs_completed_file=parentDir + PATHDELIM + 'list_jobs_completed.txt'

           try:
              if path.exists(list_jobs_submitted_file):
                  self.messagelogger.write("OK: Found job submitted list  file for sample \"%s\"!\n" %(S))
                  load_job_status_file(list_jobs_submitted_file, self.list_jobs_submitted) 
                  self.messagelogger.write("OK: Successfully loaded status from job submitted list file for sample \"%s\"!\n" %(S))
           except:
                  self.messagelogger.write("ERROR: Cannot load job submitted file for sample \"%s\"!\n" %(S))
                  return False

           try:
              if path.exists(list_jobs_completed_file):
                  self.messagelogger.write("OK: Found job completed list  file for sample \"%s\"!\n" %(S))
                  load_job_status_file(list_jobs_completed_file, self.list_jobs_completed) 
                  self.messagelogger.write("OK: Successfully loaded status from job completed list file for sample \"%s\"!\n" %(S))
           except:
                 self.messagelogger.write("ERROR: Cannot load job completed file for sample \"%s\"!\n" %(S))
                 return False

           return True

      def load_list_splits(self):
            for S in self.getSamples():
               self.load_list_splits_for_sample(S)

      def load_list_splits_for_sample(self, S):
           parentDir =  self.base_output_folder + PATHDELIM + S + PATHDELIM + 'blast_results' + PATHDELIM + 'grid'
           self.list_splits[S]={}
           read_one_column(parentDir + PATHDELIM + 'split_list.txt', self.list_splits[S], col=0)

      def load_job_lists(self):

         for S in self.getSamples():
            self.load_job_list_file(S)

         for S in self.getSamples():
            for d in self.getDBs( samples = [S] ):
               if not d in self.samples_and_databases[S]:
                  try: 
                     del self.list_jobs[S][d]
                  except KeyError:
                     pass 

      def isJobInList(self, J):
           if  J.S in self.list_jobs:
              if J.d in self.list_jobs[J.S]:
                 if  J.a in self.list_jobs[J.S][J.d]:
                    if  J.m in self.list_jobs[J.S][J.d][J.a]:
                        return True

           return False

      def isValidJobList(self, S):
            parentDir =  self.base_output_folder + PATHDELIM + S + PATHDELIM + 'blast_results' + PATHDELIM + 'grid'
            splits= {}
           
            read_one_column(parentDir + PATHDELIM + 'split_list.txt', splits, col=0)
            status = 0 # flag for no databases
            for d in self.getDBs( samples = [ S ] ):
               for a in  self.list_splits[S]:
                  status = 1 # found algorithm
                  J = Job(S,d,a, self.getAlgorithm(S))

                  if not self.isJobInList(J):
                     status = 2 # not in job list
                     return status
            return status


      def doesValidSplitExist(self, s):
            if not s in self.samples_and_inputs:
              self.messagelogger.write("WARNING: Does not have sample \"%s\" in Broker\n!" % (s))
              return False

            parentDir =  self.base_output_folder + PATHDELIM + s + PATHDELIM + 'blast_results' + PATHDELIM + 'grid'
            t = 0
            if doesFolderExist(parentDir + PATHDELIM + 'split_batch'):
               t+=1

            if doesFolderExist(parentDir + PATHDELIM + 'split_results'):
               t+=1

            if doesFileExist( parentDir + PATHDELIM + 'split_list.txt'):
               t+=1

            if t < 3:
               return False

            F = {}
            read_one_column(parentDir + PATHDELIM + 'split_list.txt', F, col=0)
      
            for f in F:
                if not doesFileExist(parentDir + PATHDELIM + 'split_batch' + PATHDELIM + f):
                   return False
             
            c1 = 0
            for f in F:
               c1 += countNoOfSequencesInFile(parentDir + PATHDELIM + 'split_batch' + PATHDELIM + f)
            c2 = countNoOfSequencesInFile(self.samples_and_inputs[s])

            if c1!=c2:
               return False

        
            return True  

      def splitInput(self, sample, force = False):
         parentDir1 =  self.base_output_folder + PATHDELIM + sample + PATHDELIM + 'blast_results' +\
                       PATHDELIM + 'grid' 

         parentDir2 =  parentDir1 + PATHDELIM + 'split_batch'
         if force:  
              clearFolderIfExists(parentDir)
         
         self.messagelogger.write("STATUS: Creating splits for file \"%s\" for sample \"%s\"!\n" % (self.samples_and_inputs[sample], sample))
         if create_splits(parentDir2, parentDir1 + PATHDELIM + 'split_list.txt', self.samples_and_inputs[sample], self.getBatchSize(),  self.batchSize, splitPrefix='split', splitSuffix=''):
            self.messagelogger.write("SUCCESS: Created splits for file \"%s\" for sample \"%s\"!\n" % (self.samples_and_inputs[sample], sample))
         else:
            self.messagelogger.write("ERROR: Cannot create splits for file \"%s\" for sample \"%s\"!\n" % (self.samples_and_inputs[sample], sample))
           
         return True
         sys.exit(0)


         checkOrCreateFolder(working_directory + PATHDELIM + 'grid')
         
         grid_folder = working_directory + PATHDELIM + 'grid' 
         batch_split_folder = working_directory + PATHDELIM + 'grid' + PATHDELIM + 'batch_splits'

         checkOrCreateFolder(batch_split_folder)

         self.list_blocks_filename = grid_folder +  PATHDELIM + 'list_blocks.txt' 

         # if it deos not exist then potentially it has not been split and vice versa
         if not path.exists(self.list_blocks_filename):
            self.create_blast_splits(batch_split_folder, self.list_blocks_filename, maxSize = batchSize )

         for gridparam in self.gridParams:
            try:
               gridservice= BlastService(gridparam, self.base_output_folder)
               gridservice.set_sample_name(self.sample_name)
               self.BlastServices.append(gridservice)
               #checkOrCreateFolder(working_directory + PATHDELIM + 'grid' + PATHDELIM + gridparam.serviceAddress)
            except:
               pass

         self.list_submitted_file = grid_folder +  PATHDELIM + 'list_submitted.txt' 
         self.completed_list_filename = grid_folder +  PATHDELIM + 'list_completed.txt' 

      def is_complete(self):
         samples_dictionary={}
         read_list(self.list_blocks_filename, samples_dictionary, col=0)
         completed_dictionary={}
         read_list(self.completed_list_filename, completed_dictionary, col=0)
         if len(samples_dictionary) == len(completed_dictionary):
            return True
         else:
            return False

      def submit_jobs(self):   
          if not path.exists(self.list_submitted_file):
              try:
                 file = open(self.list_submitted_file,'w')
                 file.close()
              except IOError:
                 print "Cannot remove file " + self.list_submitted_file + " !"
                 print "Or cannot create file " + self.list_submitted_file + " !"
                 sys.exit(0)

          if not path.exists(self.completed_list_filename):
              try:
                 file = open(self.completed_list_filename,'w')
                 file.close()
              except IOError:
                 print "Cannot remove file " +  self.completed_list_filename + " !"
                 print "Or cannot create file " + self.completed_list_filename + " !"
                 sys.exit(0)

 
          while not self.is_complete(): 
             time.sleep(1)
             for grid in self.BlastServices:
                 print "Connecting to " + grid.serviceAddress
                 grid.submit_job(self.fastaFile, self.database)
            #     grid.isUp()
                 time.sleep(1)
             return
         

      # this splits the input fasta formatted file into smaller fasta files each 
      # with size bounds maxSize (number of sequences) or maxBytes (in bytes)
      def create_blast_splits(self, target_folder, blocks_list_filename, maxSize=500, maxBytes = 40000000):
          blockno = 0 
          currblocksize = 0 
          currblockbyteSize = 0 

          fastareader = FastaReader(self.fastaFile)
          # Read sequences from sorted sequence file and write them to block files
          try:
             blocklistfile = open(blocks_list_filename, 'w')
          except:
             print "ERROR:  Cannot open " + blocks_list_filename
             sys.exit(0)

          sample_name = 'split'  
          fragments = []
          for name in fastareader:
                fragments.append(fastareader.seqname) 
                fragments.append(fastareader.sequence)
     
                if currblocksize >= maxSize -1 or currblockbyteSize >= maxBytes:
                    #TODO adjust the 000 to match the format
                    blockfile = open(target_folder +  PATHDELIM + sample_name + '.000' + str(blockno) + '.fasta', 'w')
                    fprintf(blockfile, "%s",'\n'.join(fragments))
                    fragments=[]
                    blockfile.close()
                     # Add this block name to the blocklistfile
                    #TODO adjust the 000 to match the format
                    fprintf(blocklistfile, "%s\n", sample_name + ".000" + str(blockno))
                    blockno += 1
                    currblocksize = 0 
                    currblockbyteSize = 0 
                else: 
                    currblocksize += 1
                    currblockbyteSize += len(fastareader.sequence)
     
     
          
          if fragments:
             #TODO adjust the 000 to match the format
             blockfile = open(target_folder +  PATHDELIM + sample_name + '.000' + str(blockno) + '.fasta', 'w')
             fprintf(blockfile, "%s",'\n'.join(fragments))
             blockfile.close()
             fragments = []
             #TODO adjust the 000 to match the format
             fprintf(blocklistfile, "%s\n", sample_name + ".000" + str(blockno))
             blockno += 1
     
          #Add this block name to the blocklistfile
          blocklistfile.close()
          currblocksize = 0 
          currblockbyteSize = 0 

      def getSplitResults(self, P):
          _split_list = []
          if not P[0] in self.list_jobs_completed:
             return _split_list 

          if not P[1] in self.list_jobs_completed[P[0]]:
             return _split_list 

          for a in self.list_jobs[P[0]][P[1]]:
             if a in self.list_jobs_completed[P[0]][P[1]]:
               _split_list.append( a + "." + P[1] + '.' + self.algorithm)
             else:  
               return []
          
          return _split_list



      def completed_Sample_DB_pairs(self): 
         _completed_pairs = [] 
         for S in self.list_jobs: 
           for d in self.list_jobs[S]:
              completed= True
              for a in self.list_jobs[S][d]:
                if not S in self.list_jobs_completed:
                   completed=False
                   break
                if not d in self.list_jobs_completed[S]:
                   completed=False
                   break
                if not a in self.list_jobs_completed[S][d]:
                   completed=False
                   break
                if not self.algorithm in self.list_jobs_completed[S][d][a]:
                   completed=False
                   break

              if completed:
                _completed_pairs.append( (S, d) )
         return _completed_pairs


      
      def completed_Samples(self): 
         _completed_samples = [] 
         for S in self.list_jobs: 
           completed= True
           for d in self.list_jobs[S]:
              for a in self.list_jobs[S][d]:
                if not S in self.list_jobs_completed:
                   completed=False
                   break
                if not d in self.list_jobs_completed[S]:
                   completed=False
                   break
                if not a in self.list_jobs_completed[S][d]:
                   completed=False
                   break
                if not self.algorithm in self.list_jobs_completed[S][d][a]:
                   completed=False
                   break

              if completed==False:
                 break

           if completed:
              _completed_samples.append(S)

         return _completed_samples


      def incomplete_Samples(self): 
          incompleteSamples = {}
          for S in self.getSamples():
             parentDir =  self.base_output_folder + PATHDELIM + S + PATHDELIM + 'blast_results'
             for d in self.getDBs( samples = [S]):
                blastoutfile = parentDir + PATHDELIM + S + '.' + d + '.' + self.algorithm +"out"

                if path.exists(blastoutfile):
                    count = countNoOfSequencesInFile( blastoutfile)
                    if count ==0: 
                       self.messagelogger.write("WARNING: Empty \"%s\" output results for sample \"%s\" against database \"%s\" found!\n" %(self.algorithm, S, d))
                else:
                    incompleteSamples[S]=True
          return incompleteSamples.keys()

      def isCompletelySubmitted(self, S):
         if not S in self.list_jobs: 
            return False
         for d in self.list_jobs[S]:
            if not d in  self.samples_and_databases[S]:
               continue
            for a in self.list_jobs[S][d]:
              if not S in self.list_jobs_submitted:
                  return False
              if not d in self.list_jobs_submitted[S]:
                  return False
              if not a in self.list_jobs_submitted[S][d]:
                  return False
         return True



      def get_A_Job(self, S):

         "get a job from sample S"
         if not S in self.list_jobs: 
            return None

         for d in self.list_jobs[S]:
            if not d in  self.samples_and_databases[S]:
               continue
            for a in self.list_jobs[S][d]:
               if not S in self.list_jobs_submitted:
                  job = Job(S,d,a, self.getAlgorithm(S))
                  return job
               if not d in self.list_jobs_submitted[S]:
                  job = Job(S,d,a, self.getAlgorithm(S))
                  return job
               if not a in self.list_jobs_submitted[S][d]:
                  job = Job(S,d,a, self.getAlgorithm(S))
                  return job

         return None


      def submittedSuccessfully(self, J):
          serverRanks = []

          # rank the blast services C is a service as potential candidates to submit
          for C in self.BlastServices:
              delay = self.avgWait(C.server)
              serverRanks.append( (C, delay) )

          serverRanks.sort(key = lambda tup: tup[1])

          #for T in serverRanks:
          #    print T[0].server + ' ' + str(T[1])

          
          # try to submit the job in hand to an available service

          for T in serverRanks:
             if J.server!= None and  (T[0].server==J.server or T[1] > self.avgWait(J.server) ) : 
                continue
             C = T[0]
             if C.isReadyToSubmit(self.serverLoads[C.server]):
                if C.submitJob(J):
                   self.addToSubmittedList(C.server, J)
                   self.incrementLoad(C)
                   self.setLastSubmittedServerTo(C.server)
                   return True

          #print 'Failed to Submit ' + J.S + ' ' + J.d + ' ' + J.a + ' ' + J.m + ' ' + C.server
          return None 

    
      def decrementLoad(self,  C):
         self.serverLoads[C.server] -= 1


      def incrementLoad(self,  C):
         self.serverLoads[C.server] += 1

      def __addToStatusList(self, server, J, list_file_name, list_to_add_to):
         parentDir =  self.base_output_folder + PATHDELIM + J.S + PATHDELIM + 'blast_results' + PATHDELIM + 'grid'
         list_jobs_stats_file=parentDir + PATHDELIM + list_file_name
         try:
            if not doesFileExist(list_jobs_stats_file):
                self.messagelogger.write("WARNING: Cannot file  \"%s\" for sample \"%s\"!\n" %(list_file_name, J.S))
                self.messagelogger.write("SUCCESS: Create file  \"%s\" for sample \"%s\"!\n" %(list_file_name, J.S))
                listfile  = open(list_jobs_stats_file, 'w')
                listfile.close()
         except:
            self.messagelogger.write("ERROR: Cannot open job list %s file for sample \"%s\"!\n" %(list_file_name, J.S))
            print "ERROR: Cannot open job list %s file for sample \"%s\"!\n" %(list_file_name, J.S)
            sys.exit(1)

         try:
            listfile  = open(list_jobs_stats_file, 'a')
            eventTime = int(time.time())
            fprintf(listfile, "%s\t%s\t%s\t%s\t%s\t%s\n" %(J.S, J.d, J.a, J.m, server, str(eventTime)) )
            listfile.close()
         except:
            self.messagelogger.write("ERROR: Cannot open job list %s file for sample \"%s\"!\n" %(list_file_name, J.S))
            print "ERROR: Cannot open job list %s file for sample \"%s\"!\n" %(list_file_name, J.S)
            sys.exit(1)


         
         if not J.S in list_to_add_to:
            list_to_add_to[J.S] = {}

         if not J.d in list_to_add_to[J.S]:
            list_to_add_to[J.S][J.d] = {}
        
         if not J.a in list_to_add_to[J.S][J.d]:
            list_to_add_to[J.S][J.d][J.a] = {}

         if not J.m in list_to_add_to[J.S][J.d][J.a]:
            list_to_add_to[J.S][J.d][J.a][J.m] = {}

         list_to_add_to[J.S][J.d][J.a][J.m][server] = eventTime
         return  True

      def addToSubmittedList(self, server, J):
            self.__addToStatusList(server, J, 'list_jobs_submitted.txt', self.list_jobs_submitted)

      def addToCompletedList(self, server, J):
            self.__addToStatusList(server, J, 'list_jobs_completed.txt', self.list_jobs_completed)

      def harvest(self):
          for service_name in self.getServices(): 
              C = self.getService(service_name) 
              if C.isReadyToHarvest(self.getSamples(), self.algorithm):
                 _Results = C.harvest(self.getSamples(), self.algorithm, self.list_jobs_submitted )
                 for  J in _Results:
            #       print 'Completed ' + J.S + ' ' + J.d + ' ' + J.a + ' ' + J.m + ' ' + C.server
                   self.addToCompletedList(C.server, J)
                   self.decrementLoad(C)

      def consolidateSplitResults(self, P, split_results):
          sourceParentDir =  self.base_output_folder + PATHDELIM + P[0] + PATHDELIM + 'blast_results' + PATHDELIM + 'grid' + PATHDELIM + 'split_results'
          targetParentDir =  self.base_output_folder + PATHDELIM + P[0] + PATHDELIM + 'blast_results' 
          targetFileName =  targetParentDir + PATHDELIM + P[0] + '.' + P[1] + '.' + self.algorithm +"out"

          try:
             targetfile = open( targetFileName, 'w')
          except:
             self.messagelogger.write("ERROR: Cannot create consolidated search results file %s!\n" %(targetFileName ))
             sys.exit(0)
          for filename in  split_results:
             sourceFileName = sourceParentDir + PATHDELIM + filename
             try:
                sourcefile = open(sourceFileName, 'r')
                resultLines = sourcefile.readlines()
                sourcefile.close()
             except:
                self.messagelogger.write("ERROR: Cannot create consolidated search results file %s!\n" %(sourceFileName ))
                sys.exit(0)

             try:
                for line in resultLines:
                    fprintf(targetfile, "%s", line)
             except:
                self.messagelogger.write("ERROR: Cannot write result from file %s to the consolidated file!\n" %(sourceFileName ))
                sys.exit(0)

          self.messagelogger.write("SUCCESS: Successfully consolidated search results into file %s!\n" %(targetFileName ))
          targetfile.close()

          """ Now delete the consolidates split_files files """ 
          for filename in  split_results:
             sourceFileName = sourceParentDir + PATHDELIM + filename
             os.remove(sourceFileName)


      def consolidateHarvest(self): 
           S_DB_pairs = self.completed_Sample_DB_pairs()
           for P in S_DB_pairs:  
              split_results = self.getSplitResults( P)
              print split_results
              if split_results:
                 self.consolidateSplitResults(P, split_results) 
                 self.messagelogger.write("SUCCESS: Consolidated %s and %s search results!\n"  %(P[0], P[1]))
              else: 
                 self.messagelogger.write("WARNING: Not split results to Consolidate %s and %s search results!\n"  %(P[0], P[1]))


      def isJobCompleted(self, J):
          if not J.S in self.list_jobs_completed:
             return False

          if not J.d in self.list_jobs_completed[J.S]:
             return False

          if not J.a in self.list_jobs_completed[J.S][J.d]:
             return False

          if not J.m in self.list_jobs_completed[J.S][J.d][J.a]:
             return False

          return True


      def isSomeJobPending(self):         
          "Checks if any job is pending"
          now = time.time()
          J = Job(None, None, None, None)
          time_limit = self.getAverageDelay() + 0.5*self.getStdDeviationDelay()
          time_limit = 40
          for sample in self.list_jobs_submitted:
             for db in self.list_jobs_submitted[sample]:
               if not db in  self.samples_and_databases[sample]:
                   continue
               for split in self.list_jobs_submitted[sample][db]:
                  if self.algorithm  in self.list_jobs_submitted[sample][db][split]:
                     min_submission_time = 0
                     for server in self.list_jobs_submitted[sample][db][split][self.algorithm]:
                       submission_time = self.list_jobs_submitted[sample][db][split][self.algorithm][server]
                       #print time_limit
                       if submission_time > min_submission_time:
                           min_submission_time = submission_time
                           last_submitted_server = server

                     
                     J.setValues(sample, db, split, self.algorithm, min_submission_time, last_submitted_server)
                    
                     if not self.isJobCompleted(J):
                       if now - min_submission_time  > time_limit :
                         self.pending_job = J
                         #print 'diff ' + str(now - min_submission_time)
                         #print str(now) + ' ' + str(min_submission_time) + '  ' + str(time_limit)
                         #print J.S + ' ' + J.d + ' ' + J.a + ' ' + J.m
                         return True

          #print 'no pending ' + str(time_limit)
          return False

      def get_pending_and_slow_job(self):
          if self.pending_job:
               J =  self.pending_job 
               self.pending_job = None
               return J

          return None
              

      def _count_jobs_in_list(self, S, d, list):
          if not  S in  list:
             return 0
          if not  d  in list[S]:
             return 0
          num = 0
          splits  = list[S][d].keys()
          for split in  splits:
              if  self.algorithm in list[S][d][split]:
                   num += 1
          return num

      def numTotalJobs(self, S, d):
          return self._count_jobs_in_list(S, d, self.list_jobs)

      def numSubmittedJobs(self, S, d):
          return self._count_jobs_in_list(S, d, self.list_jobs_submitted)

      def numCompletedJobs(self, S, d):
          return self._count_jobs_in_list(S, d, self.list_jobs_completed)


      def numCompletedJobsServer(self, server):
          return self._count_jobs_in_server(server, self.list_jobs_completed) 
      def numSubmittedJobsServer(self,  server):
          return self._count_jobs_in_server(server, self.list_jobs_submitted)

      def avgWait(self, service_name):
          service = self.getService(service_name)
          load = self.serverLoads[service.server]
          avg = 1
          #print ' '.join( [service.server, str(service.A), str(service.C) ])
          delay = load*avg + service.A + service.C
          return delay

      def _count_jobs_in_server(self, server,  list):
        num = 0
        for S in list:
          for  d  in list[S]:
             splits  = list[S][d].keys()
             for split in  splits:
                if self.algorithm in list[S][d][split]:
                   if server in list[S][d][split][self.algorithm]:
                     num += 1
        return num


   
      def display_stats(self):
          Samples = self.getSamples()
          allComp = 0
          allTot = 0
          allSub = 0
          allRun = 0

          samplewisedisplayStr = ''
          for S in Samples:
             sampComp = 0
             sampTot = 0
             sampSub = 0
             sampRun = 0
             dbwisedisplayStr = ''
             DBs = self.getDBs(samples = [ S ] )
             for d in DBs:
                 numTot = self.numTotalJobs(S, d)
                 numSub = self.numSubmittedJobs(S, d)
                 numComp = self.numCompletedJobs(S, d)
                 numRun = numSub - numComp

                 sampTot += numTot
                 sampSub += numSub
                 sampComp += numComp
                 sampRun += numRun

                 dbwisedisplayStr +=   ' %29s | %8s | %8s | %8s' %(d, str(numSub), str(numRun), str(numComp)) + '\n'
             samplewisedisplayStr +=  '  %28s | %8s | %8s | %8s' %(S, str(sampSub), str(sampRun), str(sampComp)) + '\n' + dbwisedisplayStr
             allTot += sampTot
             allSub += sampSub
             allComp += sampComp
             allRun += allSub - allComp
          self.displayStr = '%30s %s' %('Total jobs:', str(allTot)) + '\n'
          self.displayStr += '%30s | %8s | %8s | %8s' %('Stats', '#submted', '#running', '#comptd') + '\n' + samplewisedisplayStr
          self.displayStr += '%30s | %8s | %8s | %8s' %('All Samples', str(allSub), str(allRun), str(allComp)) + '\n' + samplewisedisplayStr

          Services = self.getServices()
          self.displayStr += '\n'
          for server in Services:
              numComp = self.numCompletedJobsServer(server)
              numSub = self.numSubmittedJobsServer(server)
              numRun= numSub - numComp
          #    print ' '.join( [server, str(numSub), str(numComp)] )
              avgwait = self.avgWait(server)
              self.displayStr += '%30s | %8s | %8s | %8s  Delay: %s' %(server[0:8], str(numSub), str(numRun), str(numComp), str(avgwait)) + '\n'
         
          self.displayStr += '\n'
       #   self.displayStr += "  %20s" %('Migration')  + ' ' +  str(self.migration) +  '\n'

          migrationMatrixStr = self.getMigrationMatrix()
          self.displayStr += migrationMatrixStr
        

          #self.displayStr = 'Total  all' 
          #sys.stdout.write("\x1b[2J\x1b[H")
          #curses.move(self.y, self.x); 
          #sys.stdout.write("\b"*self.displaylen)
          sys.stdout.write(self.displayStr)
          sys.stdout.flush()
          self.displaylen = len(self.displayStr)
          #curses.getyx(stdscr, self.y, self.x); 
           
      def getMigrationMatrix(self):
          serverList = self.migrationMatrix.keys()
          outstr = '%30s' %(' ')
          for server in serverList:
             outstr+= '%30s' %(server[0:8])
          outstr += '\n'

          for server1 in serverList:
             outstr += '%30s' %(server1[0:8])
             for server2 in serverList:
                outstr += '%30s' %(str(self.migrationMatrix[server1][server2]))
             outstr += '\n'
          return outstr

      def setupStatsVariables(self):
          serversList =  self.getServices()
          for fromServer in serversList:
            for toServer in serversList:
               if not fromServer in self.migrationMatrix:
                  self.migrationMatrix[fromServer] = {}
               self.migrationMatrix[fromServer][toServer] = 0

          Services = self.getServices()
          for server in Services:
             numSub = self.numSubmittedJobsServer(server)
             self.migrationMatrix[server][server] = numSub


      def updateMigrationCount(self, fromServer, toServer):
          if fromServer and toServer:
            self.migrationMatrix[fromServer][toServer] += 1

      def incrementTransitionMatrix(self, server) :
          if server:
            self.migrationMatrix[server][server] += 1

      def Do_Work(self):

          _A = self.incomplete_Samples()
          A = {}
          for a in _A:
             A[a] = True

          while A:
            for S in A:         
                while not self.isCompletelySubmitted(S):
                   J = self.get_A_Job(S)
                   if J == None:
                       continue

                   self.display_stats()
                   self.harvest() # harvest here before you risk getting stuck in the while loop
                   #try to submit a job
                   while not self.submittedSuccessfully(J):
                      self.display_stats()
                      self.harvest()  # now harvest
                      time.sleep(1)
                   self.incrementTransitionMatrix(self.getLastSubmittedServerTo())
                   self.updateMigrationCount(J.server, self.getLastSubmittedServerTo() )
 
            while self.isSomeJobPending():         
                J = self.get_pending_and_slow_job()  # get a slow job to migrate
                if J == None:
                    continue

                self.display_stats()
                self.harvest() # harvest here before you risk getting stuck in the while loop
                while not self.submittedSuccessfully(J): # try to submit a job
                   self.display_stats()
                   self.harvest()
                   time.sleep(1)
                self.incrementTransitionMatrix(self.getLastSubmittedServerTo())
                self.updateMigrationCount(J.server, self.getLastSubmittedServerTo() )

            self.harvest()
            self.consolidateHarvest() 
            C = self.completed_Samples()
          #  self.display_stats()
            for T in C:
              if T in A:
                del A[T]
            #print 'sleeping'
            time.sleep(1)

          return True


      def Delete_Remote_Directories(self):
          samples = self.getSamples()
          Services = self.getServices()
          for service in Services:
             server =  self.getService(service) 
             server.deleteRemoteSampleFolders(samples) 
