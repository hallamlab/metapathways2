#!/usr/bin/python

# The LowClust daemon script (that runs either remotely or locally)

# File created on 27 Jan 2012.
from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, The LowClust Project"
__credits__ = [""]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

import subprocess
import sys
from os import remove, makedirs, sys, listdir, environ, path
import re 
import inspect
from commands import getstatusoutput
from optparse import OptionParser
import shutil 
import traceback
from glob import glob
import time
import math



#config = load_config()

script_info={}
script_info['script_usage'] = []
usage= """./lowclustdaemon.py  --home-dir cwd   """
parser = OptionParser(usage)
parser.add_option("--home-dir", dest="home_dir", default ='\'\'',
                  help='home dir [REQUIRED]')
parser.add_option("--does-sample-dir-exist", dest="does_sample_dir_exist", default='',
                  help='does sample dir exist')

parser.add_option("--does-file-exist", dest="does_file_exist", default='',
                  help='does file dir exist')

parser.add_option("--create-sample-dir", dest="create_sample_dir",  default='', 
                  help='create sample dir')

parser.add_option("--cluster-blocks", dest="cluster_blocks",  default='', 
                  help='cluster the first block by comparing with the rest from BLAST results')

parser.add_option("--similarity", dest="similarity",  default=40, type='int', 
                  help='similarity cutoff [default 40%]')

parser.add_option("--create-blocks", dest="create_blocks",  default='', 
                  help='creates block of batch size')

parser.add_option("--block-size", dest="block_size",  default=10000, 
                  help='block size in #sequences')

parser.add_option("--block-mb", dest="block_mb",  default=2, 
                  help='block size in MBytes')


parser.add_option("--submit-job", dest="submit_job",  default='',   
                  help='submit job')

parser.add_option("--submit-num-jobs", dest="submit_num_jobs",  default=1, type='int',   
                  help='submit number of jobs')

parser.add_option("--max-parallel-jobs", dest="max_parallel_jobs",  default=200, type='int',   
                  help='max number of jobs that are allowed to run in parallel')

parser.add_option("--executables", dest="executables_dir",  default='',   
                  help='submit job')

parser.add_option("--location", dest="location",  default='local',   choices=['remote', 'local'], 
                  help='slave location')

parser.add_option("--is-complete", dest="is_complete",  default='',   
                  help='is current round of blasting complete')

parser.add_option("--resubmit", dest="resubmit",  default='', choices = ['', 'all', 'stalled'],   
                  help='resubmit the stalled or all the jobs [DEFAULT : stalled')

parser.add_option("--is-previous-calculation-incomplete", dest="is_previous_calculation_incomplete",  default='',   
                  help='was the previous calculation incomplete')

parser.add_option("--get-number-of-existing-blocks", dest="get_number_of_existing_blocks",  default='',   
                  help='gets the number of existing blocks by counting the items in list_blocks.txt')

parser.add_option("--get-number-of-running-jobs", dest="get_number_of_running_jobs", default='', 
                  help='get the number of running jobs')

parser.add_option("--get-number-of-completed", dest="get_number_of_completed", default='', 
                  help='get the number of completed  samples')

parser.add_option( "--algorithm", dest="algorithm", choices = ['BLAST', 'LAST'], default = "BLAST",
                   help='the algorithm used for computing homology [DEFAULT: BLAST]')

parser.add_option( "--format-database", dest="database_file", default = '', 
                   help='format the BLAST/LAST db for search for computing homology [DEFAULT: sorted.block.0]')


def fprintf(file, fmt, *args):
    file.write(fmt % args)
  
def printf(fmt, *args):
     sys.stdout.write(fmt % args)

def get_sequence_name(line):
     fields = re.split(' ', line)
     name = re.sub('>','',fields[0])
     return name

class FastaReader:
    """Parses a GenBank record from a string or file."""
    stop = False
    START_PATTERN = re.compile(r'^>')
    name = None
    future_name =None
    sequence=""
    seqname=""
    def __init__(self, fasta_filename):
        try:
            self.file = open(fasta_filename, 'r')
        except IOError:
            print "Cannot open fasta file " + fasta_filename

    def __iter__(self):
        return self

 
    def next(self):
        if self.stop:
          raise StopIteration

        try:
           if not self.name: 
               self.name = self.file.readline().strip()
           line = self.file.readline().strip()
        except:
           line = None


        if not line:
           self.stop = True
           raise StopIteration


        fragments = []
        while line and not self.START_PATTERN.search(line):
            fragments.append(line.strip()) 
            line = self.file.readline()

       # print line
        if self.future_name:
            self.name = self.future_name

        if line:
          self.future_name = line.strip()

        self.sequence =''.join(fragments)
        self.seqname = self.name.replace('>','')

        self.commentline = self.name.replace('>','')
        self.name = self.commentline.split()[0]
        
        return self.name


def get_number_of_completed(outputdir):
     completed_dictionary={}
     read_list(outputdir + '/' +  'list_completed.txt', completed_dictionary)
     print str(len(completed_dictionary))
     return 

#     try:
#        completedlistfile = open(outputdir + '/list_completed.txt', 'a')
#     except IOError:
#        print "Cannot read file " + outputdir + '/list_completed.txt' + " !"
#        sys.exit(0)
 
#     submitted_dictionary={}
#     read_list(outputdir + '/' +  'list_submitted.txt', submitted_dictionary, col=0)
#     for file in submitted_dictionary:
#         if path.exists( outputdir + '/' + file + '.blastout' ):
#              if not file in  completed_dictionary:
#                  fprintf(completedlistfile,"%s\n",file) 
#                  completed_dictionary[file]=True
#     completedlistfile.close()

#     print str(len(completed_dictionary))

def  get_number_of_running_jobs(outputdir):
     submitted_dictionary={}
     read_list(outputdir + '/' +  'list_submitted.txt', submitted_dictionary, col=0)
     completed_dictionary={}
     read_list(outputdir + '/' +  'list_completed.txt', completed_dictionary, col=0)
     print str(len(submitted_dictionary)-len(completed_dictionary))


def  remove_sample_dir(sample_dir):
     files = glob(sample_dir + '*')
     for  f in files:
        remove(f)
     if path.exists(sample_dir):
        shutil.rmtree(sample_dir)

def remove_files(dir, filenames):
     for file in filenames:
        try:
           remove(dir + '/' + file)
        except IOError:
           print "Cannot remove file  " + dir + '/' + file + " !"

def clear_files(dir, filenames):
     for file in filenames:
        try:
           remove(dir + '/' + file)
           file = open(dir + '/' + file, 'w')
           file.close()
        except IOError:
           print "Cannot remove file  " + dir + '/' + file + " !"

# Read the cluster file and add all its sequence names to the clustered_before list
def  read_clustered_sequence_names(cluster_file, clustered_before):
     try:
         clusterfile = open(cluster_file, 'r')
     except IOError:
         print "Cannot open the cluster.sequences file "
         sys.exit(0)

     lines = clusterfile.readlines()
     clusterfile.close()
     for line in lines:
         parts = [ x.strip() for x in line.strip().split('\t') ]
         if len(parts) > 0 :
            clustered_before[parts[0]] = True
     

# (Re)create the sequence blocks along with the necessary log files 
def create_blocks(outputdir, block_file_name, maxMBytes,   maxSize):

     maxBytes = 1024*1024*maxMBytes
     # list_submitted.txt records all submitted blast jobs for the current round (just the name of the query block file)
     list_submitted_file = outputdir + '/' + 'list_submitted.txt'
     try:
        if path.exists(list_submitted_file):
           remove(list_submitted_file)
        file = open(list_submitted_file,'w')
        file.close()
     except IOError:
        print "Cannot remove file " + list_submitted_file + " !"
        print "Or cannot create file " + list_submitted_file + " !"
        sys.exit(0)

     # list_completed.txt records all completed blast jobs for the current round (just the name of the query block file)
     list_completed_file = outputdir + '/' + 'list_completed.txt'
     try:
        if path.exists(list_completed_file):
           remove(list_completed_file)
        file = open(list_completed_file,'w')
        file.close()
     except IOError:
        print "Cannot remove file " + list_completed_file + " !"
        print "Or cannot create file " + list_completed_file + " !"
        sys.exit(0)

     # blocklistfile records all the blocks in the current round. Use it to delete all blocks from previous rounds.
     try:
        if path.exists(outputdir + '/' + block_file_name):
           blocklistfile = open(outputdir + '/' + block_file_name, 'r')
           blockfilenames = [ x.strip() for x in blocklistfile.readlines() ]
          
           remove_files(outputdir, blockfilenames)
           blocklistfile.close()
     except IOError:
        print "Cannot read file " + outputdir + '/' + block_file_name + " !"
        sys.exit(0)

     # clustered_before is the list of already clustered sequences
     clustered_before = {}
     if path.exists(outputdir + '/clustered.sequences.txt'):
         read_clustered_sequence_names(outputdir + '/clustered.sequences.txt', clustered_before)

     try:
        blocklistfile = open(outputdir + '/' + block_file_name, 'w')
     except IOError:
        print "Cannot read file " + outputdir + '/' + block_file_name + " !"
        sys.exit(0)


     fragments= []
     seq_beg_pattern = re.compile(">")
     blockno = 0
     currblocksize = 0
     currblockbyteSize = 0

     fastareader = FastaReader(outputdir + '/sorted.fasta')
     # Read sequences from sorted sequence file and write them to block files

     for name in fastareader:
        if not name in clustered_before:
           fragments.append('>' + fastareader.seqname) 
           fragments.append(fastareader.sequence)

           if currblocksize >= maxSize -1 or currblockbyteSize >= maxBytes:
               blockfile = open(outputdir +  '/sorted.block.' + str(blockno), 'w')
               fprintf(blockfile, "%s",'\n'.join(fragments))
               fragments=[]
               blockfile.close()
                # Add this block name to the blocklistfile
               fprintf(blocklistfile, "%s\n", 'sorted.block.' + str(blockno))
               blockno += 1
               currblocksize = 0
               currblockbyteSize = 0
           else: 
               currblocksize += 1
               currblockbyteSize += len(fastareader.sequence)


     if fragments:
        blockfile = open(outputdir +  '/sorted.block.' + str(blockno), 'w')
        fprintf(blockfile, "%s",'\n'.join(fragments))
        blockfile.close()
        fragments = []
        fprintf(blocklistfile, "%s\n", 'sorted.block.' + str(blockno))
        blockno += 1

     #Add this block name to the blocklistfile
     currblocksize = 0
     currblockbyteSize = 0

     blocklistfile.close()
     print blockno

def  read_one_column(listfilename, dictionary, col=0) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields) > col:
          dictionary[fields[col]] = True
    listfile.close()
  except:
    traceback.print_exception()

# Read the contents of a file into a dictionary (col begin with 0)
def read_list(listfilename, dictionary, col=0) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields) > col:
        dictionary[fields[0]] = fields[col]
    listfile.close()
  except:
    traceback.print_exception()
         
def add_to_listfile(listfilename, key ):
    listfile = open(listfilename, 'a')
    fprintf(listfile, "%s\n", key);
    listfile.close()

# returns the number of existing blocks by couting the number of 
# distinct items in list_blocks.txt
def get_number_of_existing_blocks(outputdir):
     if not path.exists( outputdir + '/' +  'list_blocks.txt'):
        print str(0)
        return

     samples_dictionary={}
     read_list(outputdir + '/' +  'list_blocks.txt', samples_dictionary, col=0)

     print str(len(samples_dictionary))
     return 


# Returns True if there the size of list of blocks submitted is 
# larger than in the list in complete
def is_previous_calculation_incomplete(outputdir):
     if not path.exists( outputdir + '/' +  'list_blocks.txt') or not path.exists( outputdir + '/' +  'list_completed.txt'):
        print 'no'
        return False

     samples_dictionary={}
     read_list(outputdir + '/' +  'list_blocks.txt', samples_dictionary, col=0)
     completed_dictionary={}
     read_list(outputdir + '/' +  'list_completed.txt', completed_dictionary, col=0)

     if len(samples_dictionary) != len(completed_dictionary):
         print 'yes'
         return True
     else:
         print 'no'
         return False


# Returns True if all the necessary blast jobs for this round have been completed
def is_complete(outputdir):
     samples_dictionary={}
     read_list(outputdir + '/' +  'list_blocks.txt', samples_dictionary, col=0)
     completed_dictionary={}
     read_list(outputdir + '/' +  'list_completed.txt', completed_dictionary, col=0)
     if len(samples_dictionary) == len(completed_dictionary):
         print 'yes'
         return True
     else:
         print 'no'
         return False

# Resubmits the stalled or all the jobs 
def resubmit_jobs(outputdir, type):
     print '     ' + type
     completed_dictionary={}
     read_list(outputdir + '/' +  'list_completed.txt', completed_dictionary, col=1)
     if path.exists( outputdir + '/' + 'list_submitted.txt' ):
        try:
           remove(outputdir + '/' + 'list_submitted.txt')
        except:
            print "Cannot delete file " + file 

     if type =='stalled':
        submittedfile = open(outputdir + '/list_submitted.txt', 'w')
        for file, jobid in completed_dictionary.iteritems():
           fprintf(submittedfile, "%s\n", file + '\t' +  jobid )
        submittedfile.close()

     print 'yes'
     return True

# Format the target database for blast
def format_database(outputdir, database_file, algorithm, location, executables_dir ):
    if algorithm =='BLAST':
            args = [ executables_dir +'/formatdb', '-p', 'T', '-i', (outputdir + '/' + database_file)  ]

    if algorithm =='LAST':
            args = [ executables_dir +'/lastdb', '-p', '-c', (outputdir + '/' + database_file),  (outputdir + '/' + database_file) ]

    #print ' '.join(args)
    p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = p.communicate()
    if result[1].strip()=='':
        return True
    return False


# Submit a blast job locally
def local_submit(outputdir, query, target, algorithm, executables_dir):
    if algorithm=='BLAST':
          blastcmdargs = [  executables_dir + '/blastp', '-num_threads', '10', '-max_target_seqs', 
                            '1000',  '-outfmt', '6', '-db', outputdir + '/' + target, 
                            '-query', outputdir + '/' + query, '-evalue', '0.001', '-out', ( outputdir + '/' + query + ".blastout") 
                         ]
    if algorithm=='LAST':
          blastcmdargs = [  executables_dir + '/lastal', '-o', (outputdir + '/' + query + ".blastout"), '-f', '2', 
                            outputdir + '/' + target,  outputdir + '/' + query
                         ]

    #print ' '.join(blastcmdargs)
    p = subprocess.Popen(blastcmdargs, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = p.communicate()

   # print result
    if result[0].strip()=='':
      return True
    return False

# Submit a blast job remotely
def grid_submit(outputdir, query, target, algorithm, executables_dir): 
    commandfile = open('lowclust_job_qsub.txt', 'w')
    if algorithm=='BLAST':
         blastcmdargs = [  executables_dir + '/blastp', '-num_threads', '1', '-max_target_seqs', 
                          '5',  '-outfmt', '6', '-db', outputdir + '/' + target, 
                         '-query', outputdir  + '/' + query, '-evalue', '0.001', '-out', (outputdir +'/' + query + ".blastout") 
                        ]

    if algorithm=='LAST':
          blastcmdargs = [  executables_dir + '/lastal', '-o', (outputdir + '/' + query + ".blastout"), '-f', '2', 
                            outputdir + '/' + target,  outputdir + '/' + query
                         ]

    blastCommand =  ' '.join(blastcmdargs) 
    commandfile = open(outputdir + '/lowclust_job_qsub.txt', 'w')
    fprintf(commandfile, "%s\n",blastCommand)
    fprintf(commandfile, "%s\n","echo \"#EOF\" >> " + (outputdir +'/' + query + ".blastout"))
    fprintf(commandfile, "%s\n","echo \"#EOF\" >> " + outputdir +'/.qstatdir/' + query + '.job')
    commandfile.close()

#    args = [ 'qsub', '-l', 'walltime=10:00:00', '-m',  'ea',  '-j',  'eo',  
#             '-e' ,outputdir +'/.qstatdir/' + query + '.job', '-l',  'procs=1',\
#             '-l', 'mem=10gb', outputdir + '/lowclust_job_qsub.txt'
#          ]

    args = [ 'qsub',   
           #  '-e' ,outputdir +'/.qstatdir/' + query + '.job', 
              '-e', outputdir + '/.qstatdir/error.log', 
              '-o', outputdir + '/.qstatdir/output.log', 
              outputdir + '/lowclust_job_qsub.txt'
          ]

    p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = p.communicate()
    #print "result =", result
    returncode = p.returncode
    if returncode == 0:
       print "adding to submitted list"
       add_to_listfile(outputdir + '/' + 'list_submitted.txt', query + '\t' + query + '.job' + '\t' + str(int(time.time())))
    else:
       print "not adding to submitted list"
       return False

    return True

def read_completed_task(statsDir, jobid_dictionary):
    files = glob(statsDir +'*')
    for file in files:
       shortfile = re.sub(r'^.*/','',file)
       jobid_dictionary[shortfile] = 1
       #print 'a' + file 
    #   remove(file)
    #print jobid_dictionary
    return


def remove_stats_file(statsDir, blockname):
     statFilePath =statsDir +'/' + blockname + '.job'
     if path.exists(statFilePath):
        remove(statFilePath)

def is_Homology_Search_Complete(boutFileName):
    try: 
        if path.exists(boutFileName):
            filehandle = open(boutFileName, 'r')
            lines = filehandle.readlines()
            if len(lines) > 0 and  re.search(r'#EOF', lines[len(lines) -1]):
               return True
            else:
               return False
    except:
         return False
     

# Perform one round of clustering, taking new representatives from the first block and recruiting
# members from all of the blocks (based on similarity from blast output)
def cluster_blocks(outputdir, similarity, algorithm):
     samples_dictionary={}
     read_list(outputdir + '/' +  'list_blocks.txt', samples_dictionary, col=0)

     try:
        block0file = open(outputdir + '/sorted.block.0', 'r')
     except IOError:
        print "Cannot read file " + outputdir + '/sorted.block.0' + " !"
        sys.exit(0)

     # clusters maps cluster representatives to a list of their cluster members
     clusters = {}
     # clustered_before is the list of already clustered sequences
     clustered_already = {}
     if path.exists(outputdir + '/clustered.sequences.txt'):
         read_clustered_sequence_names(outputdir + '/clustered.sequences.txt', clustered_already)

     # get_blast_results maps each sequence in block0file to a list of tuples containing matching sequences and their similarity score
     print "reading blast results...."
     get_blast_results = {}
     for file in samples_dictionary.keys():
         retrieve_blast_results(outputdir, file+".blastout", get_blast_results, algorithm)
     print "done"
     
     seq_beg_pattern = re.compile(">")

     # Iterate over sequences in block0file
     print "recruiting cluster members...."
     for line in block0file.readlines():
        line = line.strip()
     
        if seq_beg_pattern.search(line):
            name=get_sequence_name(line)
        else:
            continue

        # If sequence has already been clustered, skip it
        if name in clustered_already:
           continue

        # ...Otherwise, select it as a new representative
        clusters[name] =[]
        
        # Add all the matching sequences as its cluster members
        if name in get_blast_results:
           for result in get_blast_results[name]:
              if (not result[0] in clustered_already) and  (similarity < result[1]):
                  clusters[name].append(result)
                  clustered_already[result[0]] =True
     block0file.close()
     print "done"
     
     # At this point we have already formed the new clusters      
     try:
        clusterfile = open(outputdir + '/clustered.sequences.txt', 'a')
     except IOError:
        print "Cannot read file " + outputdir + '/clusters.txt' + " !"
        sys.exit(0)

     # Write the clusters to the clustered sequences file
     for cluster in clusters:
        fprintf(clusterfile, "%s\tC\t%d\n",cluster, 100)
        for memberAndSim in clusters[cluster]:
          fprintf(clusterfile, "%s\tM\t%d\n",memberAndSim[0], memberAndSim[1])
     clusterfile.close()
     print "done"
     
     # At this point we can safely delete the blast output files
     delete_files_in_folder_with_suffix(outputdir, "blastout")
     # Also delete the list_submitted.txt and list_completed.txt files 
     clear_files(outputdir, ['list_submitted.txt', 'list_completed.txt'])


def delete_files_in_folder_with_suffix(outputdir, suffix):
     files = glob(outputdir + '/' + '*' + suffix)
     try:
        for  f in files:
           remove(f)
     except:
        pass

#    for file in completed_dictionary:
#         if path.exists( outputdir + '/' + file + '.blastout' ):
#             try:
#               remove(outputdir + '/' + file + '.blastout')
#             except:
#               print "Cannot delete file " + file 
#               sys.exit(0)

  
                  
        
def  retrieve_blast_results(outputdir, blastoutput, get_blast_results, algorithm):
     comment_pattern = re.compile(r"^#")

     try:
        blastoutfile = open(outputdir + '/' + blastoutput, 'r')
     except IOError:
        print "Cannot read file " + outputdir + '/' + blastoutput + " !"
        sys.exit(0)

     for line in blastoutfile.readlines():
         if comment_pattern.search(line):
            continue

         fields = [ x.strip() for x in line.strip().split('\t')]
         if (algorithm=='BLAST' and len(fields) < 11)  or (algorithm=='LAST' and  len(fields) != 3) :
            continue

         source = fields[0] 
         target = fields[1] 
         similarity = float(fields[2])

         if not target in get_blast_results:
            get_blast_results[target] = []

         if source==target:
              continue
         get_blast_results[target].append( (source, similarity) )
     blastoutfile.close()

def harvest_completed_jobs_local(outputdir):
    submitted_dictionary={}
    read_list(outputdir + '/'  + 'list_submitted.txt', submitted_dictionary, col=0)
    #print submitted_dictionary

    completed_dictionary={}
    read_list(outputdir + '/' +  'list_completed.txt', completed_dictionary, col=0)

    for file in submitted_dictionary:
        if path.exists( outputdir + '/' + file + '.blastout' ):
            if not file in  completed_dictionary:
               add_to_listfile(outputdir + '/' + 'list_completed.txt', file)


def harvest_completed_jobs_remote(outputdir):
    submitted_dictionary={}
    read_list(outputdir + '/'  + 'list_submitted.txt', submitted_dictionary, col=1)
    #print submitted_dictionary
    statsDir = outputdir +  '/.qstatdir/'
    jobid_dictionary={}
    read_completed_task(statsDir, jobid_dictionary)

    for block in submitted_dictionary:
       if submitted_dictionary[block]  in jobid_dictionary:
          if is_Homology_Search_Complete(outputdir + '/' + block + '.blastout'): 
              remove_stats_file(statsDir, block)
              add_to_listfile(outputdir + '/' + 'list_completed.txt', block + '\t' + submitted_dictionary[block] + '\t' + str(int(time.time())))


def submit_job(outputdir, location, algorithm, max_parallel_jobs,  executables_dir, submit_num_jobs):
   if location=='remote':
      harvest_completed_jobs_remote(outputdir)

   if location=='local':
      harvest_completed_jobs_local(outputdir)

   try:
     samples_dictionary={}
     read_list(outputdir + '/' +  'list_blocks.txt', samples_dictionary, col=0)
     submitted_dictionary={}
     read_list(outputdir + '/' +  'list_submitted.txt', submitted_dictionary, col=2)
     completed_dictionary={}
     read_list(outputdir + '/' +  'list_completed.txt', completed_dictionary, col=2)

     numrunning_jobs = len(submitted_dictionary.keys()) - len(completed_dictionary.keys())
     if numrunning_jobs >= max_parallel_jobs:
        return

    # if blastdb_format(outputdir, 'sorted.block.0', location) == False:
    #     print "Failed to formatdb file " + key
    #     return

     jobs_submitted = 0
     for key in samples_dictionary:
       if not key  in submitted_dictionary:
           if location=='local':
               if local_submit(outputdir, key, 'sorted.block.0', algorithm, executables_dir) == False:
                 print "Failed to blast locally file " + key
               add_to_listfile(outputdir + '/' + 'list_submitted.txt', key)

           if location=='remote':
              if grid_submit(outputdir, key, 'sorted.block.0', algorithm,  executables_dir) == False:
                 print "Failed to blast grid file " + key

           jobs_submitted+=1
           if jobs_submitted >= submit_num_jobs:
             break

     if location=='remote': 
        key = get_one_pending_job(submitted_dictionary, completed_dictionary )
        if key and grid_submit(outputdir, key, 'sorted.block.0', algorithm,  executables_dir) == False:
             print "Failed to blast grid file " + key

   except:
       traceback.print_exception()

   print "yes"


def get_one_pending_job(submitted_dictionary, completed_dictionary):
     time_intervals = []
     for key in submitted_dictionary:
       if  key  in completed_dictionary:
         time_intervals.append(int(completed_dictionary[key]) - int(submitted_dictionary[key]) )
     
     if not time_intervals: 
         return None

     avgtsq = 0
     avgt = 0
     for t in time_intervals:  
         avgtsq += t*t
         avgt += t
     avgtsq = avgtsq/len(time_intervals)
     avgt = avgt/len(time_intervals)
     stdev = math.sqrt(avgtsq - avgt*avgt)
        
     currtime = int( time.time())
     for key in submitted_dictionary:
       if not key  in completed_dictionary:
           if currtime - int(submitted_dictionary[key]) - avgt > 10*stdev:
              return key

     return None
# checks if the supplied arguments are adequate
def isValid(opts, args):
    if (opts.home_dir==None):
       return True
    else:
       return False

def doesFileExist(file):
    if path.exists(file):
       return True
    else:
       return False

def doesFilePattExist(filepatt):
    files = glob(filepatt)
    if len(files)>0:
       return True
    else:
       return False

def doesSampleDirExists(sampledir):
    if path.exists(sampledir):
       return True
    else:
       return False


def extendPath( folders):
    newdir = ''
    for folder in folders:
       if len(folder)!=0:
        newdir =  newdir + folder  + '/'
    return newdir

# Check that folders and files are set up, then call the appropriate function depending on the arguments passed in
def main(argv):
    (opts, args) = parser.parse_args()
    if isValid(opts, args):
       print usage
       sys.exit(0)

    # initialize the input directory or file
    if len(opts.does_sample_dir_exist):
        sampledir = extendPath([opts.home_dir, opts.does_sample_dir_exist])
        if not doesSampleDirExists(sampledir):
           print 'no'
        else:
           print 'yes'

    # create the sample directory
    if len(opts.create_sample_dir)>0:
        try:
            sampledir = extendPath([opts.home_dir, opts.create_sample_dir])
            makedirs(sampledir)
            print 'yes'
        except:
            print 'no'

    # create the sample directory
    if len(opts.does_file_exist)>0:
        try:
            #file = extendPath([ opts.home_dir, opts.create_sample_dir,opts.dies_file_exist] )
            if doesFileExist(opts.does_file_exist): 
               print 'yes'
            else:
               print 'no'
        except:
            print 'no'

    if len(opts.cluster_blocks) > 0 :
       try:
          clustetednum = cluster_blocks(opts.cluster_blocks, opts.similarity, opts.algorithm)
       except:
          return 0

   
    if len(opts.create_blocks) > 0:
       try:
          batch_count = create_blocks(opts.home_dir, 'list_blocks.txt', float(opts.block_mb), int(opts.block_size))
       except:
          return 0

    # format the first bock command
    if len(opts.executables_dir) > 0 and len(opts.database_file) > 0:
       try:
          format_database(opts.home_dir, opts.database_file, opts.algorithm, opts.location, opts.executables_dir )
       except:
          return 0

    if len(opts.executables_dir) > 0 and len(opts.submit_job) > 0 and opts.location in ['remote', 'local']:
       try:
          submit_job(opts.submit_job, opts.location, opts.algorithm, opts.max_parallel_jobs, opts.executables_dir, opts.submit_num_jobs)
       except:
          return 0

    #get the number of running jobs
    if len(opts.get_number_of_running_jobs)>0:
       try:
          get_number_of_running_jobs(opts.get_number_of_running_jobs)
       except:
          return 0

    #get number of completed samples 
    if len(opts.get_number_of_completed)>0:
       try:
          get_number_of_completed(opts.get_number_of_completed)
       except:
          return 0

    # check is completed 
    if len(opts.is_complete) :
      try:
         is_complete(opts.is_complete)
      except:
         return 0

    # check is there an incomplete calculation from a previous round 
    if len(opts.is_previous_calculation_incomplete) :
      try:
         is_previous_calculation_incomplete(opts.home_dir)
      except:
         return 0

    # check is there an incomplete calculation from a previous round 
    if len(opts.get_number_of_existing_blocks) :
      try:
         get_number_of_existing_blocks(opts.home_dir)
      except:
         return 0

    # resubmit the jobs
    if len(opts.resubmit) :
      try:
         resubmit_jobs(opts.home_dir, opts.resubmit)
      except:
         return 0

# the main function of lowclustdaemon
if __name__ == "__main__":
    main(sys.argv[1:])    

