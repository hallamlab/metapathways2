#!/usr/bin/python

"""This script run the pathologic """

try:
   import optparse, sys, re, csv, traceback
   from os import path, _exit
   import logging.handlers
   from glob import glob
   import multiprocessing

   from libs.python_modules.utils.sysutil import pathDelim
   from libs.python_modules.utils.metapathways_utils  import fprintf, printf, eprintf,  exit_process
   from libs.python_modules.utils.sysutil import getstatusoutput

   from libs.python_modules.utils.pathwaytoolsutils import *

except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     print traceback.print_exc(10)
     sys.exit(3)


PATHDELIM= pathDelim()

def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

def files_exist( files , errorlogger = None):
    status = True    
    for file in files:
       if not path.exists(file):
          if errorlogger:
             errorlogger.write( 'ERROR\tCould not find ptools input  file : ' +  file )
          status = False
    return not status



usage = sys.argv[0] + """ -c <contigs> -o <output> -r <reads>  -O <orfgff> --rpkmExec <rpkmexec> """
parser = None
def createParser():
    global parser

    epilog = """This script computes the RPKM values for each ORF, from the BWA 
                recruits. 
             """

    epilog = re.sub(r'\s+', ' ', epilog)

    parser = optparse.OptionParser(usage=usage, epilog = epilog)

    # Input options


    parser.add_option('-c', '--contigs', dest='contigs', default=None,
                           help='the contigs file')

    parser.add_option('-o', '--output', dest='output', default=None,
                           help='orfwise RPKM file')

    parser.add_option('--stats', dest='stats', default=None,
                           help='output stats for ORFs  into file')

    parser.add_option('-r', '--rpkmdir', dest='rpkmdir', default=None,
                           help='the directory that should have the read files')

    parser.add_option('-O', '--orfgff', dest='orfgff', default=None,
                           help='folder of the PGDB')

    parser.add_option('-s', '--sample_name', dest='sample_name', default=None,
                           help='name of the sample')

    parser.add_option('--rpkmExec', dest='rpkmExec', default=None,
                           help='RPKM Executable')

    parser.add_option('--bwaExec', dest='bwaExec', default=None,
                           help='BWA Executable')

    parser.add_option('--bwaFolder', dest='bwaFolder', default=None,
                           help='BWA Folder')


def getSamFiles(readdir, sample_name):
   '''This function finds the set of fastq files that has the reads'''

   samFiles = []
   _samFile = glob(readdir + PATHDELIM + sample_name + '.sam')

   if _samFile:
      samFiles += _samFile

   _samFiles = glob(readdir + PATHDELIM + sample_name + '_[0-9].sam')

   if _samFiles:
     samFiles += _samFiles

   return samFiles



def getReadFiles(readdir, sample_name):
   '''This function finds the set of fastq files that has the reads'''

   fastqFiles = []

   _fastqfiles = glob(readdir + PATHDELIM + sample_name + '_[Rr][0-9].[fF][aA][Ss][Tt][qQ]')
   if _fastqfiles:
      fastqFiles += _fastqfiles

   _fastqfiles = glob(readdir + PATHDELIM + sample_name + '_[0-9].[fF][aA][Ss][Tt][qQ]')
   if _fastqfiles:
      fastqFiles += _fastqfiles


   _fastqfiles = glob(readdir + PATHDELIM + sample_name + '_[0-9].[fF][qQ]')
   if _fastqfiles:
      fastqFiles += _fastqfiles

   _fastqfiles = glob(readdir + PATHDELIM + sample_name + '_[Rr][0-9].[fF][qQ]')
   if _fastqfiles:
      fastqFiles += _fastqfiles

   _fastqfiles = glob(readdir + PATHDELIM + sample_name + '.[fF][aA][Ss][Tt][qQ]')
   if _fastqfiles:
      fastqFiles += _fastqfiles

   _fastqfiles = glob(readdir + PATHDELIM + sample_name + '.[fF][qQ]')
   if _fastqfiles:
      fastqFiles += _fastqfiles

   return fastqFiles



def indexForBWA(bwaExec, contigs,  indexfile):
    cmd = "%s index -p %s %s"  %(bwaExec, indexfile, contigs, )

    result = getstatusoutput(cmd)

    if result[0]==0:
       return True

    return False


def runUsingBWA(bwaExec, sample_name, indexFile,  readFiles, bwaFolder) :

    if len(readFiles) > 2:
       return False

    num_threads =  int(multiprocessing.cpu_count()*0.8)
    if num_threads < 1:
       num_threads = 1

    bwaOutput = bwaFolder + PATHDELIM + sample_name + '.sam'

    if len(readFiles) == 2:
       cmd = "%s mem -t %d -o %s %s %s %s"  %(bwaExec, num_threads,  bwaOutput, indexFile,  readFiles[0], readFiles[1])

    if len(readFiles) == 1:
       cmd = "%s mem -t %d -p -o %s  %s %s "  %(bwaExec, num_threads,   bwaOutput, indexFile,  readFiles[0])

    result = getstatusoutput(cmd)

    if result[0]==0:
       return True

    return False



def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser

    options, args = parser.parse_args(argv)
    if not (options.contigs!=None and  path.exists(options.contigs)):
       parser.error('ERROR\tThe contigs file is missing')
       return 255

    if not (options.rpkmExec !=None and path.exists(options.rpkmExec) ) :
       parser.error('ERROR\tThe RPKM executable is missing')
       return 255 

    if not (options.bwaExec !=None and path.exists(options.bwaExec) ) :
       parser.error('ERROR\tThe BWA executable is missing')
       return 255 

    if not (options.rpkmdir !=None and path.exists(options.rpkmdir) ):
       parser.error('ERROR\tThe RPKM directory is missing')
       return 255 

    if not (options.bwaFolder !=None and path.exists(options.bwaFolder) ):
       parser.error('ERROR\tThe BWA directory is missing')
       return 255 

    if  options.sample_name==None :
       parser.error('ERROR\tThe sample name is missing')
       return 255 


    # read the input sam and fastq  files
    samFiles = getSamFiles(options.rpkmdir, options.sample_name)
    readFiles = getReadFiles(options.rpkmdir, options.sample_name)

    if not samFiles and readFiles:
        if not readFiles:
           eprintf("ERROR\tCannot find the read files not found for sample %s!\n", options.sample_name)
           eprintf("ERROR\tMetaPathways need to have the sample names in the format %s.fastq or (%s_1.fastq and %s_2.fastq) !\n", options.sample_name, options.sample_name, options.sample_name)
           if errorlogger:
             errorlogger.eprintf("ERROR\tCannot find the read files not found for sample %s!\n", options.sample_name)
             errorlogger.eprintf("ERROR\tMetaPathways need to have the sample names in the format %s.fastq or (%s_1.fastq and %s_2.fastq) !\n", options.sample_name, options.sample_name, options.sample_name)
             return 255
    
        # index for BWA
        bwaIndexFile = options.bwaFolder + PATHDELIM + options.sample_name
        indexSuccess = indexForBWA(options.bwaExec, options.contigs, bwaIndexFile) 

        if not indexSuccess:
           eprintf("ERROR\tCannot index the preprocessed file %s!\n", options.contigs)
           if errorlogger:
              errorlogger.eprintf("ERROR\tCannot index the preprocessed file %s!\n", options.contigs)
           return 255
           #exit_process("ERROR\tMissing read files!\n")
    
    
        #print 'running'

        bwaRunSuccess = runUsingBWA(options.bwaExec, options.sample_name,  bwaIndexFile, readFiles, options.bwaFolder) 
        print 'run success', bwaRunSuccess
        if not bwaRunSuccess:
           eprintf("ERROR\tCannot successfully run BWA for file %s!\n", options.contigs)
           if errorlogger:
              errorlogger.eprintf("ERROR\tCannot successfully run BWA for file %s!\n", options.contigs)
           return 255
           #exit_process("ERROR\tFailed to run BWA!\n")


    # is there a RPKM executable installed
    print 'rpkm running', bwaRunSuccess
    if not path.exists(options.rpkmExec):
       eprintf("ERROR\tRPKM executable %s not found!\n", options.rpkmExec)
       if errorlogger:
          errorlogger.printf("ERROR\tRPKM executable %s not found!\n",  options.rpkmExec)
       return 255
       #exit_process("ERROR\tRPKM executable %s not found!\n" %(options.rpkmExec))


    # command to build the RPKM
    command = "%s --c %s"  %(options.rpkmExec, options.contigs)
    command += " --multireads " 
    if options.output:
       command += " --ORF-RPKM %s" %(options.output)
       command += " --stats %s" %(options.stats)

    if options.orfgff:
       command += " --ORFS %s" %(options.orfgff)

    samFiles = getSamFiles(options.bwaFolder, options.sample_name)

    if not samFiles:
       return 0

    for samfile in samFiles:
        command += " -r " + samfile

    try:
       status  = runRPKMCommand(runcommand = command) 
    except:
       status = 1
       pass

    if status!=0:
       eprintf("ERROR\tRPKM calculation was unsuccessful\n")
       return 255
       #exit_process("ERROR\tFailed to run RPKM" )

    return status

def runRPKMCommand(runcommand = None):
    if runcommand == None:
      return False

    result = getstatusoutput(runcommand)
    return result[0]


# this is the portion of the code that fixes the name

def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

     return attributes


# this is the function that fixes the name
def  fix_pgdb_input_files(pgdb_folder, pgdbs = []):
     pgdb_list = glob(pgdb_folder + '/*/input/organism.dat')     

     for pgdb_organism_file in pgdb_list:
        process_organism_file(pgdb_organism_file)


def fixLine(line, id):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[0]+'\t' + id
     

def getID(line):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[1]
     
def process_organism_file(filel):
     patternsToFix = [ re.compile(r'NAME\tunclassified sequences'), re.compile(r'ABBREV-NAME\tu. sequences') ]
     patternID =  re.compile(r'^ID\t.*')
     try:
         orgfile = open(filel,'r')
     except IOError:
         print "ERROR : Cannot open organism file" + str(filel)
         return 

     lines = orgfile.readlines()
     newlines = []
     needsFixing = False

     id = None
     for line in lines:
         line = line.strip()
         if len(line)==0:
            continue
         flag = False

         result = patternID.search(line)
         if result:   
             id = getID(line)
          
         for patternToFix in patternsToFix:
             result = patternToFix.search(line)
             if result and id:
                 newline = fixLine(line, id)
                 newlines.append(newline)
                 flag= True
                 needsFixing = True

         if flag==False:
            newlines.append(line)

     orgfile.close()
     if needsFixing:
       write_new_file(newlines, filel)


def write_new_file(lines, output_file):
    
    print "Fixing file " + output_file 
    try:
       outputfile = open(output_file,'w')
       pass
    except IOError:
         print "ERROR :Cannot open output file "  + output_file
   
    for line in lines:
       fprintf(outputfile, "%s\n", line)

    outputfile.close()


def MetaPathways_rpkm(argv, extra_command = None, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\tRPKM_CALCULATION\n")
    createParser()
    returncode = main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (returncode,'')

if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

