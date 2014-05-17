#!/usr/bin/python
# File created on 27 Jan 2012.
from __future__ import division

__author__ = "Kishori M Konwar, Niels W Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = [""]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar, Niels W Hanson"
__status__ = "Release"

import subprocess
import sys
from os import remove, makedirs, sys, listdir, environ, path, system
import re 
import inspect
from optparse import OptionParser
import shutil 
import traceback
import string, random
from glob import glob



#config = load_config()
metapaths_config = """template_config.txt""";

script_info={}
script_info['script_usage'] = []
usage= """./daemon  --home-dir cwd   """
parser = OptionParser(usage)
parser.add_option("--home-dir", dest="home_dir", default ='\'\'',
                  help='home dir [REQUIRED]')
parser.add_option("--does-sample-dir-exist", dest="does_sample_dir_exist", default='',
                  help='does sample dir exist')
parser.add_option("--does-file-exist", dest="does_file_exist", default='',
                  help='does file dir exist')

parser.add_option("--does-file-pattern-exist", dest="does_file_pattern_exist", default='',
                  help='does file with the pattern exist')

parser.add_option("--number-of-lines-in-file", dest="number_of_lines_in_file", default='', 
                  help='counts the number of lines in the file')

parser.add_option("--get-files-with-pattern", dest="get_files_with_pattern", default='',
                  help='get file names with the pattern supplied')

parser.add_option("--get-result-filenames-with-suffix_in_dir", dest="get_result_filenames_with_suffix_in_dir", default='',
                  help='get file names with the suffix')

parser.add_option("--get-result-filenames", dest="get_result_filenames", default='',
                  help='get the result file names')


parser.add_option("--eof-tag", dest="eof_tag", default='',
                  help='the end of file tag')

parser.add_option("--suffix", dest="suffix", default='',
                  help='the suffix at the end of the filename')

parser.add_option("--sample-name", dest="sample_name",  default='', 
                  help='sample name')

parser.add_option("--create-sample-dir", dest="create_sample_dir",  default='', 
                  help='create sample dir')

parser.add_option("--remove-sample-dir", dest="remove_sample_dir",  default='', 
                  help='remove sample directory')

parser.add_option("--remove-sample-dirs", dest="remove_sample_dirs",  default=[], action = "append", 
                  help='remove sample directories')

parser.add_option("--remove-file", dest="remove_file",  default='', 
                  help='remove file')

parser.add_option("--number-of-sequences-in-file", dest="number_of_sequences_in_file",  default='', 
                  help='count the number of sequences in file')

parser.add_option("--split-into-batches", dest="split_into_batches",  default='', 
                  help='count the number of sequences in file')

parser.add_option("--format-database", dest="format_database",  default='', 
                  help='formats the database')

parser.add_option("--batch-size", dest="batch_size",  default=500, 
                  help='batch size')

parser.add_option("--os-type", dest="os_type", action= 'store_true',  default=False,   
                  help='return OS type')

parser.add_option("--cpu-type", dest="cpu_type",  default='',   
                  help='return CPU type')

parser.add_option("--algorithm", dest="algorithm", choices = ['BLAST', 'LAST'], default = "BLAST",
                  help='the algorithm used for computing homology [DEFAULT: BLAST]')

parser.add_option("--submit-job", dest="submit_job",  default='',   
                  help='submit job')

parser.add_option("--submission-type", dest="submission_type",  choices=['0', '1'],  default='0',   
                  help='submission type qsub = 0, userscritpt  = 1')

parser.add_option("--submit-string", dest="sub_string",  default='',   
                  help='system dependent header, such as #PBS -l')

parser.add_option("--memory", dest="memory",  default='10gb',   
                  help='memory size request')

parser.add_option("--walltime", dest="walltime",  default='10:00:00',   
                  help='wall time request')

parser.add_option("--database-file", dest="database_file",  default='',   
                  help='database file name')

parser.add_option("--database-files", dest="database_files",  default=[], action='append',   
                  help='database file names')

parser.add_option("--dbname", dest="dbname",  default='',   
                  help='dbname')

parser.add_option("--dbnames", dest="dbnames",  default=[], action='append',   
                  help='dbnames, an array')

parser.add_option("--is-complete", dest="is_complete",  default='',   
                  help='is consolidate')

parser.add_option("--consolidate", dest="consolidate",  default='',   
                  help='consolidate with dbname and sample name')

parser.add_option("--get-number-of-running-jobs", dest="get_number_of_running_jobs", default='', 
                  help='get the number of running jobs')

parser.add_option("--get-number-of-samples", dest="get_number_of_samples", default='', 
                  help='get the number of samples')

parser.add_option("--get-number-of-completed", dest="get_number_of_completed", default='', 
                  help='get the number of completed  samples')


def fprintf(file, fmt, *args):
    file.write(fmt % args)
  
def printf(fmt, *args):
     sys.stdout.write(fmt % args)

def get_sequence_name(line):
     fields = re.split(' ', line)
     name = re.sub('>','',fields[0])
     #print name
     return name

def get_number_of_completed(sample_name, algorithm):
     namePrefix = 'MetaPathways/' + sample_name + '/' + sample_name +'_' 
     samples_dictionary={}
     read_list(namePrefix + 'completed.txt', samples_dictionary, filtercol=2, filterword=algorithm)
     print str(len(samples_dictionary))

def get_number_of_samples(sample_name, algorithm):
     namePrefix = 'MetaPathways/' + sample_name + '/' + sample_name +'_' 
     samples_dictionary={}
     read_list(namePrefix + 'samples.txt', samples_dictionary, col=1, filtercol=4, filterword=algorithm)
     print str(len(samples_dictionary))
     return len(samples_dictionary)

def _get_number_of_lines_in_file(filename):
    try:
       listfile = open(filename, 'r')
       lines = listfile.readlines()
       listfile.close()
       return len(lines)
    except:
       return 0


def  get_number_of_running_jobs(login):
     args = [ 'qstat', '-u', login]
     p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
     result = p.communicate() 
     lines = result[0].strip().split('\n')
     num_running_jobs = 0 
     loginPattern = re.compile( login )

     for line in lines:
        if loginPattern.match(line): 
           num_runnning_jobs += 1 
         
     if num_running_jobs==0:   
        num_running_jobs = str(len(lines)-4) 

     print str(num_running_jobs)


def  remove_file(filename, working_dir='~'):
     try:
        if path.exists(filename):
          remove(filename)
        print '<<Success!>>' 
     except:
        print '<<Failure!>>' 

def  remove_sample_dirs(sample_dirs, working_dir='~'):
    try:
      for sample_dir in sample_dirs:
        if path.exists(sample_dir):
          files = glob(sample_dir + '/*')
          for  f in files:
             remove(f)
          if path.exists(sample_dir):
            shutil.rmtree(sample_dir)
      print '<<Success!>>' 
    except:
        print '<<Failure!>>' 
            


def  remove_sample_dir(sample_dir):
     files = glob(sample_dir + '/*')
     for  f in files:
        remove(f)
     if path.exists(sample_dir):
        shutil.rmtree(sample_dir)

def create_the_sequence_files(sequence_file_name, samplefilename, sample_name, database_file, dbname,  size, algorithm):
     try:
        sequencefile = open(sequence_file_name, 'r')
     except IOError:
        print "Cannot read file " + sequence_file_name + " !"

     sequence_lines = sequencefile.readlines()

     sequencefile.close()
     fragments= []
     name=""

     seq_dictionary={}
     seq_beg_pattern = re.compile(">")
     for line in sequence_lines:
        line = line.strip()
        if seq_beg_pattern.search(line):
          if len(name) > 0:
             sequence=''.join(fragments)
             seq_dictionary[name]=sequence
             fragments = []
          name=get_sequence_name(line)
        else:
          fragments.append(line)

     if len(name) > 0:
        sequence=''.join(fragments)
        seq_dictionary[name]=sequence


     samplefile = open(samplefilename , 'w')

     count =0
     filecount=0
     smallfilename = sample_name +'_' + str(filecount) +'.faa'
     fprintf(samplefile,'%s_%s_%s_%s\t%s\t%s\t%s\t%s\n',algorithm, dbname, database_file, smallfilename, dbname,\
             database_file, smallfilename, algorithm)

     smallfile = open(smallfilename , 'w')
     for name in seq_dictionary:
       if count %size==0  and count>0:
         smallfile.close()
         filecount+=1
         smallfilename = sample_name +'_' + str(filecount) +'.faa'
         smallfile = open(smallfilename , 'w')
         fprintf(samplefile,'%s_%s_%s_%s\t%s\t%s\t%s\t%s\n',algorithm, dbname, database_file, smallfilename,\
                 dbname, database_file, smallfilename, algorithm)

       fprintf(smallfile,'>%s\n',name)
       fprintf(smallfile,'%s\n',seq_dictionary[name])
       count+=1

     smallfile.close()
     samplefile.close()

def already_split_for_dbname(filename, dbname, database_filename, algorithm): 

     dbnames_file_algorithm=[]
     
     read_list_tuples(filename, dbnames_file_algorithm, [1, 2, 4])
     for db, dbfile, alg in  dbnames_file_algorithm:
        if db==dbname and dbfile==database_filename and algorithm == alg:
           return True

     return False


def split_into_batches(sequence_file_name, database_files,  dbnames, size, algorithm):
     sample_name = re.sub(r'[.]qced[.]faa','',sequence_file_name)
     submittedfilename = sample_name +'_submitted.txt'
     completedfilename = sample_name +'_completed.txt'
     samplefilename = sample_name +'_samples.txt'

     if not path.exists(samplefilename):
        if len(dbnames) > 0 and len(database_files)>0:
           create_the_sequence_files(sequence_file_name, samplefilename, sample_name, database_files[0], dbnames[0],size, algorithm)
           submittedfile = open(submittedfilename , 'w')
           submittedfile.close()
           completedfile = open(completedfilename , 'w')
           completedfile.close()


     if not path.exists(submittedfilename): 
        submittedfile = open(submittedfilename , 'w')
        submittedfile.close()

     if not path.exists(completedfilename):
        completedfile = open(completedfilename , 'w')
        completedfile.close()


     for database_file, dbname in zip(database_files, dbnames):
         if already_split_for_dbname(samplefilename, dbname, database_file, algorithm): 
            print "already split for " + dbname
            continue

         samples_filename_dictionary={}
         
         read_one_column(samplefilename, samples_filename_dictionary, col=3)

         samplefile = open(samplefilename , 'a')
         for smallfilename in samples_filename_dictionary:
            fprintf(samplefile,'%s_%s_%s_%s\t%s\t%s\t%s\t%s\n',algorithm, dbname, database_file, smallfilename, dbname, database_file, smallfilename, algorithm)
            #printf('%s_%s_%s\t%s\t%s\t%s\n',dbname, database_file, smallfilename, dbname, database_file, smallfilename)
         samplefile.close()


     return _get_number_of_lines_in_file(samplefilename)

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
    traceback.print_exc(1)

# col begin with 0
def read_filtered_list(listfilename, dictionary,  col=1, filter_cols=[], filter_by_words=[]):
  try:
    if not path.exists(listfilename):
      return

    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       status = True
       for filcol, words  in zip(filter_cols, filter_by_words):
          if fields[filcol]!=words:
             status = False
       if status == False:
           continue

       if len(fields) > col:
        dictionary[fields[0]] = fields[col]
    listfile.close()
  except:
    print 'traceback'
    traceback.print_exc(1)


def read_list(listfilename, dictionary, col=1, filtercol=-1, filterword ='') :
  try:
    if not path.exists(listfilename):
      return

    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if filtercol >-1 :
          if len(fields) <= filtercol or  fields[filtercol]!=filterword:
               continue
       if len(fields) > col:
        dictionary[fields[0]] = fields[col]
    listfile.close()
  except:
    print 'traceback'
    traceback.print_exc(1)

# col begin with 0
def  read_list_tuples(listfilename, list_return, cols) :
  unique =dict( (col, True) for col in cols)
  maxcol = max(cols)
  
  try:
    if not path.exists(listfilename):
      return
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields) <= maxcol:
             continue
       tuples = [] 
       for col in cols:
          tuples.append(fields[col] ) 
       list_return.append(tuple(tuples))
    listfile.close()
  except:
    traceback.print_exc(1)


def  read_list_reverse(listfilename, dictionary) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields)==1:
          dictionary[fields[0]] = 1
       if len(fields)==2:
          dictionary[fields[1]] = fields[0]
    listfile.close()
  except:
    traceback.print_exc(1)
         
def add_to_listfile(listfilename, key, id ):
    listfile = open(listfilename, 'a')
    fprintf(listfile, "%s\t%s\n", key, id);
    listfile.close()

def read_completed_task(statsDir, jobid_dictionary):
    files = glob(statsDir +'*')
    for file in files:
       shortfile = re.sub(r'^.*/','',file)
       jobid_dictionary[shortfile] = 1
       #print 'a' + file 
       remove(file)
    #print jobid_dictionary
    return

def  number_of_lines_in_file(filename, working_dir ='~'):
     try:  
        file = open(filename, 'r')
        lines = file.readlines()
        file.close()
        size = len(lines)
     except:   
        return 0
     print str(size)                                                             
     return size

def format_database(database, algorithm, working_dir = '~'):
     if algorithm=='LAST':
       dbformatter = working_dir + '/' + 'MetaPathways/executables/lastdb' 
       args = [ dbformatter, '-p', '-c', working_dir + '/' +  database, working_dir + '/' + database ]

     if algorithm=='BLAST':
       dbformatter = working_dir + '/' + 'MetaPathways/executables/makeblastdb' 
       args = [ dbformatter, '-input_type', 'fasta', '-in', working_dir  + '/' + database, '-dbtype', 'prot' ]

     p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
     result = p.communicate()

     if result[1].strip()=='':
        print '<<Success!>>' 
        return True
     else:
        print '<<Failure!>>' 
        return False


def  _create_dictionary_of_arrays(listfilename, dictionary, col1=1, col2=2) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields) > col1 and len(fields) > col2:
           if not fields[col1] in  dictionary:
              dictionary[fields[col1]] = []
           dictionary[fields[col1]].append(fields[col2])
    listfile.close()
  except:
    traceback.print_exc(1)

def  append_file_content_from_to(sourcefilename, targetfilename):
    sourcefile = open(sourcefilename, 'r')
    targetfile = open(targetfilename, 'a')
    sourcelines = sourcefile.readlines()
    for line in sourcelines:
        fprintf(targetfile, "%s\n", line.strip())
    sourcefile.close()
    targetfile.close()


def consolidate(sample_name, dbname, algorithm):
     namePrefix = 'MetaPathways/' + sample_name + '/' + sample_name +'_' 
     by_dbnames={}
     _create_dictionary_of_arrays(namePrefix + 'samples.txt', by_dbnames, col1=1, col2=3)

     if algorithm == 'BLAST':
        db_output_suffix =  dbname + '.blastout'
     if algorithm == 'LAST':
        db_output_suffix =  dbname + '.lastout'

     consolidatedfile = 'MetaPathways/' + sample_name + '/' + sample_name +'.' + db_output_suffix

     if path.exists(consolidatedfile):
       remove(consolidatedfile)
     for smallfile in by_dbnames[dbname]:
         #print smallfile + ' ' + consolidatedfile
          try:
              append_file_content_from_to(smallfile + '.' + db_output_suffix, consolidatedfile)
          except:
              pass
              #traceback.print_exc(1)

     print '<<Success!>>' 
     return True

def is_complete(sample_name, dbname, algorithm):
     namePrefix = 'MetaPathways/' + sample_name + '/' + sample_name +'_' 
     samples_dictionary={}
     read_filtered_list(namePrefix + 'samples.txt', samples_dictionary, col=1, filter_cols =  [1, 4], filter_by_words=[dbname, algorithm])
     completed_dictionary={}
     read_filtered_list(namePrefix + 'completed.txt', completed_dictionary, col=1, filter_cols= [2, 3], filter_by_words = [algorithm, dbname])
     if len(samples_dictionary) == len(completed_dictionary):
         print '<<Success!>>' 
         return True
     else:
         print '<<Failure!>>' 
         return False


def submit_job(sample_name, split_file, dbname,  algorithm, working_dir = '~', sub_string=None, submission_type = '0'):

   #working_dir = '/home/sgeadmin'
   try:
     sampleDir = 'MetaPathways/samples/' + sample_name 
     commandfileName =  sampleDir+ '/'  + 'metapaths_job_qsub.txt'
     #commandfileName = path.abspath(__commandfileName)

     try:
        if path.exists(commandfileName):
            remove(commandfileName)
        commandfile = open(commandfileName, 'w+')
     except IOError, e:
         print "I/O error({0}): {1}".format(e.errno, e.strerror)
        

     shortoutputfilename =   split_file+ "." + dbname +  "." + algorithm
     outputfile =  sampleDir + '/' +  shortoutputfilename
     if algorithm=='BLAST': 
          command = working_dir + '/' + 'MetaPathways/executables/blastp -num_threads 1  -max_target_seqs 5  -outfmt 6' +\
                        ' -db ' +  ('MetaPathways/databases/' + dbname)  +\
                        ' -query ' + sampleDir + '/' + split_file  + ' -evalue 0.000001 '+\
                        ' -out '  +  outputfile

     if algorithm=='LAST':
          command = working_dir + '/' +  'MetaPathways/executables/lastal' +\
                        ' -o '  +  outputfile +\
                        ' -f 0 '  +  ('MetaPathways/databases/' + dbname)  +\
                        ' ' + sampleDir + '/' + split_file


     #fprintf(commandfile, "%s\n","#$ -wd " + working_dir)
    
     
     if sub_string!=None:
         fprintf(commandfile, "#%s\n",re.sub(',', ' ',sub_string))

     fprintf(commandfile, "%s\n",command)
     fprintf(commandfile, "%s\n","echo \"   \" >> " + outputfile)
     fprintf(commandfile, "%s\n","echo \"#EOF\" >> " + outputfile)
     fprintf(commandfile, "%s\n","echo \"#EOF\" >> " + sampleDir +'/.qstatdir/' + shortoutputfilename + '.job')
     commandfile.close()

     
     random_id= id_generator()
     args = [ "qsub",
   #            "-wd", working_dir,   #  '-e' ,outputdir +'/.qstatdir/' + query + '.job', 
               "-e", sampleDir + "/.qstatdir/" + shortoutputfilename + random_id +  "error.log",
               "-o", sampleDir + "/.qstatdir/" + shortoutputfilename + random_id + "output.log",
                commandfileName
             ]

     #print ' '.join(args) 
     if submission_type=='0':
        p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result = p.communicate()
        returncode = p.returncode
        sendResultBack(result, returncode)

     if submission_type=='1':
        new_job_file = open("MetaPathways/new_job.txt", 'w+')
        new_job_file.close()
        sendResultBack( ['',''], 0)

   except:
      pass

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

def sendResultBack(result, returncode):
   try:
     stdout_msg = result[0]
     stdout_err = result[1]
     if returncode == 0:
        print '<<Success!>>',stdout_msg,stdout_err
     else:  
        print '<<Failure!>>',stdout_msg,stdout_err
   except:
        print '<<Failure!>>',stdout_msg, stdout_err, traceback.print_exc(1)

def os_type():
    try:
       from platform import machine, platform
    except:
       print "unknown"
       return
    
    if platform():
       print platform()

def number_of_sequences_in_file(sequence_file_name):
     try:
        sequencefile = open(sequence_file_name, 'r')
     except IOError:
        print "Cannot read file " + sequence_file_name + " !"
     sequence_lines = sequencefile.readlines()
     sequencefile.close()
     name=""
     count =0
     seq_beg_pattern = re.compile(">")
     for line in sequence_lines:
        line = line.strip()
        if seq_beg_pattern.search(line):
          count+=1
     return count


# checks if the supplied arguments are adequate
def isValid(opts, args):
    if (opts.home_dir==None):
       return True
    else:
       return False


# in the folder and then delete the folder too
def removeSampleDir( origFolderName):
    folderName = origFolderName + '/*'
    files = glob(folderName)
    for  f in files:
       remove(f)
    if path.exists(origFolderName):
      shutil.rmtree(origFolderName)


def doesFileExist(_file):
    fileName = path.abspath(re.sub(r'~/', '', _file))
    if path.exists(fileName):
       return True
    else:
       return False

def doesFilePatternExist(filepatt):
    files = glob(filepatt)
    if len(files)>0:
       return True
    else:
       return False

def  isValidResultFile(FileName, suffix, eof_tag):
     eof_PATTERN = re.compile(eof_tag)
     try:
        if path.exists(FileName):
           filehandle = open(boutFileName, 'r')
           lines = filehandle.readlines()
           filehandle.close()
            
           if len(lines) > 0 and  eof_PATTERN.search(lines[len(lines) -1]):
              return True
           else:
              return False
     except:
         return False


def getResultFileNames(folder, suffix, eof_tag ): 
    try:
        resultnamesfile = open(folder + '/' + 'list_results.txt')
        result_lines = resultnamesfile.readlines()
        resultnamesfile.close()
        names = []
        for line in lines:
           names += [ x.strip() for x in  line.split(',')]
    except:
        print "Cannot open  file %s" %(folder + '/' + 'list_results.txt')

    try:
        result_names = []
        for name in names:
          if isValidResultFile(folder + '/' + name, suffix, eof_tag):
             result_names.append(name)
        print ' '.join(names)
    except:
        print "Cannot read some result file in %s" %(folder)


def getFileNamesWithSuffix(folder, suffix, eof_tag) :
    files = glob(folder + '/' + '\*' + suffix)
    basenames = [ path.basename(x) for x in files]
    print ','.join(basenames)


def getFilesWithPattern(filepatt):
    fields = filepatt.split('__DELIMITER__') 
    folder = path.abspath(fields[0]).replace('~','')
    suffix = fields[1]
    files = glob(folder + '//*' + suffix)
    basenames = [ path.basename(x) for x in files]
    print ','.join(basenames)


def doesSampleDirExists(sampledir):
    #dir= path.abspath(re.sub(r'~/', '', sampledir))
    dir = sampledir
    if path.exists(dir):
       print "<<Success!>>" + dir
    else:
       print "<<Failure!>>" + dir


def extendPath( folders):
    newdir = ''
    for folder in folders:
       if len(folder)!=0:
        newdir =  newdir + folder  + '/'
    return newdir

def main(argv):
    (opts, args) = parser.parse_args()
    if isValid(opts, args):
       print usage
       sys.exit(0)


    # initialize the input directory or file
    if len(opts.does_sample_dir_exist):
        sampledir = opts.does_sample_dir_exist
        #sampledir = extendPath([opts.home_dir, opts.does_sample_dir_exist])
        doesSampleDirExists(sampledir)

    # create the sample directory
    if len(opts.create_sample_dir)>0:
        try:
            #sampledir = extendPath([opts.home_dir, opts.create_sample_dir])
            sampledir = opts.create_sample_dir
            makedirs(sampledir)
            print '<<Success!>>' 
        except:
            print '<<Failure!>>' 

    # remove the sample directory
    if len(opts.remove_sample_dir)>0:
        try:
            sampledir = opts.remove_sample_dir
            #sampledir = extendPath([opts.home_dir, opts.remove_sample_dir])
            remove_sample_dir(sampledir)
            print '<<Success!>>' 
        except:
            print '<<Failure!>>' 

    # remove the sample directories
    if len(opts.remove_sample_dirs)>0:
        remove_sample_dirs(opts.remove_sample_dirs,  working_dir = opts.home_dir )


    # remove the sample directory
    if len(opts.remove_file)>0:
        remove_file(opts.remove_file, working_dir=opts.home_dir)

    # create the sample directory
    if len(opts.does_file_exist)>0:
        try:
            #file = extendPath([ opts.home_dir, opts.create_sample_dir,opts.dies_file_exist] )
            if doesFileExist(opts.does_file_exist): 
               print '<<Success!>>' 
            else:
               print '<<Failure!>>' 
        except:
             print '<<Failure!>>' 
     
    # does file with the pattern exist
    if len(opts.does_file_pattern_exist)>0:
        try:
            if doesFilePatternExist(opts.does_file_pattern_exist): 
               print '<<Success!>>' 
            else:
               print '<<Failure!>>' 
        except:
             print '<<Failure!>>' 

    if len(opts.get_files_with_pattern)>0:
        try:
            getFilesWithPattern(opts.get_files_with_pattern)
        except:
            print ' '

    if len(opts.get_result_filenames) > 0: 
        try:
            getResultFileNames(opts.get_result_filenames, opts.filename_suffix,  opts.eof_tag) 
        except:
            print ' '


    if len(opts.get_result_filenames_with_suffix_in_dir) > 0: 
        try:
            getFileNamesWithSuffix(opts.get_result_filenames_with_suffix_in_dir, opts.filename_suffix,  opts.eof_tag) 
        except:
            print ' '

    if len(opts.number_of_sequences_in_file) > 0:
        count = number_of_sequences_in_file(opts.number_of_sequences_in_file)
        print str(count)
    #print opts


    if len(opts.split_into_batches) > 0:
       try:
          batch_count = split_into_batches(opts.split_into_batches, opts.database_files,\
                         opts.dbnames,  int(opts.batch_size), opts.algorithm)
          print batch_count
       except:
          return 0

    if len(opts.submit_job) > 0 and len(opts.dbname) and len(opts.sample_name) and  len(opts.algorithm):
       try:
          submit_job(opts.sample_name, opts.submit_job, opts.dbname, opts.algorithm, working_dir = '~', sub_string = opts.sub_string, submission_type = opts.submission_type)
       except:
          return 0

    if opts.os_type==True:
       try:
          os_type()
       except:
          return 0

    # format database
    if len(opts.format_database) > 0:
       try:
          format_database(opts.format_database, opts.algorithm, working_dir =opts.home_dir)
       except:
          return 0

    # get the number of lines in a file
    if len(opts.number_of_lines_in_file) > 0:
       try:
          number_of_lines_in_file(opts.number_of_lines_in_file, working_dir =opts.home_dir)
       except:
          return 0

    #get the number of running jobs
    if len(opts.get_number_of_running_jobs)>0:
       try:
          get_number_of_running_jobs(opts.get_number_of_running_jobs)
       except:
          return 0

    #get number of samples 
    if len(opts.get_number_of_samples)>0:
       try:
          get_number_of_samples(opts.get_number_of_samples, opts.algorithm)
       except:
          return 0

    #get number of completed samples 
    if len(opts.get_number_of_completed)>0:
       try:
          get_number_of_completed(opts.get_number_of_completed, opts.algorithm)
       except:
          return 0

    # check is completed 
    if len(opts.is_complete) :
      try:
         is_complete(opts.is_complete, opts.dbname, opts.algorithm)
      except:
         return 0

    # check is completed 
    if len(opts.consolidate)>0 and len(opts.dbname)>0: 
      try:
         consolidate(opts.consolidate, opts.dbname, opts.algorithm)
      except:
         return 0

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])    

