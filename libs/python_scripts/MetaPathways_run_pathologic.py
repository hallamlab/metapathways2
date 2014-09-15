#!/usr/bin/python

"""This script run the pathologic """

try:
   import optparse, sys, re, csv, traceback
   from os import path, _exit
   import logging.handlers
   from glob import glob

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

PATHDELIM=pathDelim()



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



usage = sys.argv[0] + """ -i input_folder -p pgdb_dir --ptoolsExec pathwaytools_executable """
parser = None
def createParser():
    global parser

    epilog = """The pathway prediction algorithm, Pathologic in the Pathway Tools software, is run with the folder ptools as the input. The result of this step is an ePGDB (environmental pathway genome database).
The resulting ePGDB is in the ~/ptools-local/pgdbs/user folder. They can be viewed using the Pathway Tools software."""

    epilog = re.sub(r'\s+', ' ', epilog)

    parser = optparse.OptionParser(usage=usage, epilog = epilog)

    # Input options


    parser.add_option('-i', '--input', dest='inputfolder', default=None,
                           help='runs pathologic on the input folder')

    parser.add_option('-s', '--sample', dest='sample_name', default=None,
                           help='sample name')

    parser.add_option('-r', '--reactions', dest='reactions_list', default=None,
                           help='creates the metacyc reaction lists extracted from PGDB')

    parser.add_option('-p', '--pgdb', dest='pgdbdir', default=None,
                           help='folder of the PGDB')

    parser.add_option('--ptoolsExec', dest='ptoolsExec', default=None,
                           help='PathoLogic Executable')

    parser.add_option('--no-taxonomic-pruning', dest='no_taxonomic_pruning', default=True,
                           help='Option to stop taxonomic pruning')

    parser.add_option('--no-web-cel-overview', dest='no_web_cel_overview', default=True,
                           help='Option to turn off cellular overview')




def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser

    options, args = parser.parse_args(argv)
    if options.inputfolder ==None:
       parser.error('ERROR\tInput folder for Pathologic not found')
    else:
      # required files to be able to build ePGDB
      files = [ 
                options.inputfolder + PATHDELIM + '0.pf',
              #  options.inputfolder + PATHDELIM + '0.fasta',
                options.inputfolder + PATHDELIM + 'genetic-elements.dat',  
                options.inputfolder + PATHDELIM + 'organism-params.dat'
              ]

      if files_exist( files , errorlogger = errorlogger):
        exit_process("ERROR\tCannot find all inputs for Pathologic in folder %s : "  %(options.inputfolder) )

    # is there a pathwaytools executable installed
    if not path.exists(options.ptoolsExec):
       eprintf("ERROR\tPathwayTools executable %s not found!\n", options.ptoolsExec)
       if errorlogger:
          errorlogger.printf("ERROR\tPathwayTools executable %s not found!\n",  options.ptoolsExec)
       exit_process("ERROR\tPathwayTools executable %s not found!\n" %(options.ptoolsExec))


    # command to build the ePGDB
    command = "%s -patho %s"  %(options.ptoolsExec, options.inputfolder)
    if options.no_taxonomic_pruning:
       command += " -no-taxonomic-pruning "

    if options.no_web_cel_overview:
       command += " -no-web-cel-overview"

    command += " -api"

    status =0

    
    fix_pgdb_input_files(options.pgdbdir, pgdbs = [])

    if not path.exists(options.pgdbdir): 
      status  = runPathologicCommand(runcommand = command) 
      fix_pgdb_input_files(options.pgdbdir, pgdbs = [])


    if status!=0:
       eprintf("ERROR\tFailed to run Pathologic on input %s : \n" %(options.inputfolder))
       eprintf("INFO\tKill any other PathwayTools instance running on the machine and try again\n")
       if errorlogger:
          errorlogger.write("ERROR\tFailed to run Pathologic on input %s : " %(options.inputfolder))
          errorlogger.write("INFO\tKill any other PathwayTools instance running on the machine and try again")
          errorlogger.write("     : " + command)
       exit_process("ERROR\tFailed to run Pathologic on input %s : "  %(options.inputfolder) )
    
    try:
        pythonCyc = PythonCyc()
        pythonCyc.setOrganism(options.sample_name.lower())
        pythonCyc.setPToolsExec(options.ptoolsExec)
        pythonCyc.startPathwayTools()

        resultLines = pythonCyc.getReactionListLines()
        pythonCyc.stopPathwayTools()

        reaction_list_file = open(options.reactions_list, 'w')
        for line in resultLines:
           fprintf(reaction_list_file,"%s\n",line.strip())

        reaction_list_file.close()
    except:
        eprintf("ERROR\tFailed to run extract pathways for %s : \n" %(options.sample_name))
        eprintf("INFO\tKill any other PathwayTools instance running on the machine and try again")
        if errorlogger:
            errorlogger.write("ERROR\tFailed to run extract pathways for %s : " %(options.sample_name))
            errorlogger.write("INFO\tKill any other PathwayTools instance running on the machine and try again\n")
        pass 


def runPathologicCommand(runcommand = None):
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


def MetaPathways_run_pathologic(argv, extra_command = None, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\tBUILD_PGDB\n")
    createParser()
    main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (0,'')

if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

