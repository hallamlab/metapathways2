#!/usr/bin/python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


try:
   from optparse import make_option
   from os import makedirs,  path, listdir, remove, rename, _exit
   import os, sys, errno, shutil, re
   from glob import glob
   from datetime import date
   #from metapaths_utils  import pars[s._command_line_parameters

   from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
   #from libs.python_modules.utils.utils import *, hasInput, createFolderIfNotFound
   from libs.python_modules.utils.utils import *
   from libs.python_modules.parsers.parse  import parse_metapaths_parameters
   from libs.python_modules.pipeline.metapathways_pipeline import print_commands,  execute_tasks
   from libs.python_modules.pipeline.MetaPathways_gather_run_stats import MetaPathways_gather_run_stats
   from libs.python_modules.utils.metapathways_utils import fprintf, printf, eprintf, remove_existing_pgdb, exit_process, WorkflowLogger, generate_log_fp

   from libs.python_modules.pipeline.sampledata import *
   from libs.python_modules.pipeline.jobscreator import *
   from libs.python_modules.pipeline.commands import *
   import libs.python_scripts


except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)



PATHDELIM = pathDelim()


def copyFile(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

def dry_run_status( commands ):
    for command in commands:
        printf("%s", command[0])
        if command[4] == True:
           printf("%s", " Required")
        else:
           printf("%s", " Not Required")
    printf("\n")


def get_refdb_name( dbstring ):
    dbstring = dbstring.rstrip()
    dbstring = dbstring.lstrip()
    dbstring = dbstring.lower() 
    return dbstring



def format_db(formatdb_executable, seqType, raw_sequence_file, formatted_db,  algorithm):
     _temp_formatted_db  =  formatted_db+ "__temp__"

     """ format with 4GB file size """
     if algorithm=='BLAST':
         cmd='%s -dbtype %s --max_file_sz 4294967296  -in %s -out %s' %(formatdb_executable, seqType, raw_sequence_file, _temp_formatted_db)

     if algorithm=='LAST':
         # dirname = os.path.dirname(raw_sequence_file)    
         cmd='%s -s 4G -p -c %s  %s' %(formatdb_executable, _temp_formatted_db, raw_sequence_file)

     result= getstatusoutput(cmd)
     temp_fileList = glob(_temp_formatted_db + '*') 
     try:
        for tempFile in temp_fileList:
           file = re.sub('__temp__','', tempFile)
           rename( tempFile, file);

     except:
        return False

     if result[0]==0:
        return True
     else:
        return False


def create_quality_check_cmd(min_length, type,  log, contig_lengths_file, input, output, mapping, config_settings):
    cmd = "%s  --min_length %d --log_file %s  -i %s -o  %s -M %s -t %s -L %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['PREPROCESS_FASTA']), min_length, log, input,  output,  mapping, type,   contig_lengths_file) 
    #cmd = "%s  --min_length %d --log_file %s  -i %s -o  %s -M %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['PREPROCESS_FASTA']), min_length, log, input,  output, mapping) 
    return cmd


def create_orf_prediction_cmd(options, log, input, output, config_settings):
    cmd = "%s %s -i %s -o %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['ORF_PREDICTION']), options, input, output) 
    return cmd

# convert an input gbk file to fna faa and gff file
def  convert_gbk_to_fna_faa_gff(input_gbk, output_fna, output_faa, output_gff, config_settings):
    cmd = "%s  -g %s --output-fna %s --output-faa %s --output-gff %s" %((config_settings['METAPATHWAYS_PATH'] \
                 + config_settings['GBK_TO_FNA_FAA_GFF']), input_gbk, output_fna, output_faa, output_gff) 
    return cmd

# convert an input gff file to fna faa and gff file
def  convert_gff_to_fna_faa_gff(inputs, outputs,  config_settings):
    cmd = "%s " %(config_settings['METAPATHWAYS_PATH']+ config_settings['GFF_TO_FNA_FAA_GFF'])
    for  source, target in zip(inputs, outputs):
       cmd += ' --source ' + source + ' --target ' + target
    return cmd


#converts the orfs into 
#def create_orf_converter_cmd(options, input, output, config_settings):
#    cmd = "%s  %s  %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['ORF_CONVERTER']), options,  input) 
#    return cmd

def create_aa_orf_sequences_cmd(input_gff, input_fasta, output_fna, output_faa, output_gff,config_settings):
    cmd = "%s -g  %s  -n %s --output_nuc %s --output_amino %s --output_gff %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['GFF_TO_FASTA']), input_gff, input_fasta, output_fna, output_faa, output_gff) 
    return cmd

# create genbank files
def create_genebank_file_cmd(output, input_orf, input_fas, sample_name, config_settings) :
    cmd = "%s  --out %s --orphelia %s --fasta %s  --sample-name %s"\
           %((config_settings['METAPATHWAYS_PATH'] + config_settings['CREATE_GENBANK']), output, input_orf, input_fas, sample_name) 
    return cmd

# extract from the genbank files a fasta file
def create_genebank_2_fasta_cmd(options, input, output, config_settings):
    cmd = "%s %s --o %s %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['GENBANK_TO_FASTA']), options, output, input) 
    return cmd

def create_create_filtered_amino_acid_sequences_cmd(input, type, output, orf_lengths_file,  logfile, min_length, config_settings):
    cmd = "%s  --min_length %d --log_file %s  -i %s -o  %s -L %s -t %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['PREPROCESS_FASTA']), min_length, logfile, input,  output, orf_lengths_file, type) 
    return cmd 

#compute refscores
def create_refscores_compute_cmd(input, output, config_settings, algorithm):
    if algorithm == 'BLAST':
       cmd = "%s  -F %s -B %s -o %s -i %s -a  %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['COMPUTE_REFSCORE']), \
             (config_settings['METAPATHWAYS_PATH'] + config_settings['FORMATDB_EXECUTABLE']), \
             (config_settings['METAPATHWAYS_PATH'] + config_settings['BLASTP_EXECUTABLE']), \
              output, input, algorithm) 

    elif algorithm == 'LAST':
       cmd = "%s  -F %s -B %s -o %s -i %s -a %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['COMPUTE_REFSCORE']), \
             (config_settings['METAPATHWAYS_PATH'] + config_settings['LASTDB_EXECUTABLE']), \
             (config_settings['METAPATHWAYS_PATH'] + config_settings['LAST_EXECUTABLE']), \
              output, input, algorithm) 
    else:
        eprintf("Unknown choice for annotation algorithm : " + "\"" + algorithm + "\"" + "\n Please check your parameter file\n")
        exit_process("Unknown choice for annotation algorithm : " + "\"" + algorithm + "\"" + "\n Please check your parameter file\n")

    return cmd

def formatted_db_exists(dbname, suffixes):
    for suffix in suffixes:
       allfileList = glob(dbname + '*.' + suffix) 
       fileList = []
       tempFilePattern = re.compile(r''+ dbname + '\d*.' + suffix +'$');

       for aFile in allfileList:
           searchResult =  tempFilePattern.search(aFile)
           if searchResult:
             fileList.append(aFile)

       if len(fileList)==0 :
          eprintf("ERROR :  if formatted correctely then expected the files with pattern %s\n", dbname + suffix)
          return False

    return True

def check_if_raw_sequences_exist(filename):
    return path.exists(filename)


# Makes.re that the ref database is formatted for blasting 
def check_an_format_refdb(dbname, seqType,  config_settings, params, globallogger = None): 



    algorithm=  get_parameter( params,'annotation','algorithm').upper()
    
    suffixes=[]
    
    # we do not use LAST for searchingin the taxonomic databas. e.g., greengenes, silva, etc
    # if the db formatting request is done with nucl and LAST, we switch to BLAST-based formatting
    if algorithm == 'LAST' and seqType == 'nucl':
       algorithm = 'BLAST'
    
    if algorithm == 'LAST' and seqType == 'prot':
        suffixes = [ 'des', 'sds', 'suf', 'bck', 'prj', 'ssp', 'tis' ]
    
    if algorithm == 'BLAST':
      if seqType=='prot':
        suffixes = ['phr', 'psq', 'pin']
    
      if seqType=='nucl':
        suffixes = ['nhr', 'nsq', 'nin']
        
    # formatted DB directories
    taxonomic_formatted = config_settings['REFDBS'] + PATHDELIM + 'taxonomic' + PATHDELIM + 'formatted'
    functional_formatted = config_settings['REFDBS'] + PATHDELIM + 'functional' + PATHDELIM + 'formatted'
    # check if formatted folder exis. if not create it
    for d in [taxonomic_formatted, functional_formatted]:
        if not createFolderIfNotFound(d):
            eprintf("WARNING : Creating formatted subdirectory in blastDB folder.\n")
    
    # formatted database output paths
    if seqType == 'nucl':
       seqPath= config_settings['REFDBS'] + PATHDELIM + 'taxonomic' + PATHDELIM +  dbname
       formattedDBPath = taxonomic_formatted + PATHDELIM +  dbname
    elif seqType == 'prot':
       seqPath = config_settings['REFDBS'] + PATHDELIM + 'functional' + PATHDELIM +  dbname
       formattedDBPath = functional_formatted + PATHDELIM +  dbname
    else:
       eprintf("ERROR : Undefined sequnce type for %s!\n", dbname) 
       if globallogger!=None:
          globallogger.write("ERROR \t Undefined sequnce type for %s!\n" %( dbname) )
       exit_process()
    
    # database formatting executables paths
    if algorithm == 'LAST' and seqType =='prot':
      executable  = config_settings['METAPATHWAYS_PATH'] + PATHDELIM + config_settings['LASTDB_EXECUTABLE']
    else: # algorithm == 'BLAST':
      executable  = config_settings['METAPATHWAYS_PATH'] + PATHDELIM + config_settings['FORMATDB_EXECUTABLE']

    if not (formatted_db_exists(formattedDBPath,  suffixes) ):
        eprintf("WARNING : You do not seem to have Database %s formatted!\n", dbname)
        if globallogger!=None:
          globallogger.write("WARNING\t You do not seem to have Database %s formatted!\n" %(dbname) )
        if check_if_raw_sequences_exist(seqPath):
            eprintf("          Found raw sequences for  Database %s in folder %s!\n", dbname, seqPath)
            eprintf("          Trying to format on the fly .... for %s!\n", algorithm )
            if globallogger!=None:
               globallogger.write("WARNING\t Found raw sequences for  Database %s in folder %s!\n" %(dbname, seqPath) )
               globallogger.write("Trying to format on the fly .... for %s!\n" %(algorithm ) )

            result =format_db(executable, seqType, seqPath, formattedDBPath, algorithm)
            if result ==True:
                eprintf("          Formatting successful!\n")
                return 
            else:
                eprintf("          Formatting failed! Please consider formatting manually or do not try to annotate with this database!\n")
                if globallogger!=None:
                  globallogger.write("ERROR\tFormatting failed! Please consider formatting manually or do not try to annotate with this database!\n")
                exit_process()

        eprintf("ERROR : You do not even have the raw sequence for Database %s to format!\n", dbname)
        eprintf("        in the folder %s\n", seqPath)
        eprintf("        Please put the appropriate files in folder \"blastDB\"\n")
        if globallogger!=None:
            globallogger.write("ERROR \t You do not even have the raw sequence for Database %s to format!\n" %( dbname) )
            globallogger.write("in the folder %s\n" %(seqPath))
            globallogger.write("Please put the appropriate files in folder \"blastDB\"\n")
        exit_process()
  

#    fullRefDbMapName = config_settings['METAPATHWAYS_PATH'] + PATHDELIM + config_settings['REFDBS'] +PATHDELIM + dbname + '-names.txt'
#    if not doFilesExist( [fullRefDbMapName ] ):
#        s.exit(0)


def  make_sure_map_file_exists(config_settings, dbname, globallogger = None):
    dbmapFile = config_settings['REFDBS'] + PATHDELIM + 'functional' + PATHDELIM + 'formatted' + PATHDELIM + dbname + "-names.txt"
    seqFilePath = config_settings['REFDBS'] + PATHDELIM + 'functional' + PATHDELIM + dbname
    if not doFilesExist( [dbmapFile ] ):
         eprintf("WARNING: Trying to create database map file for %s\n", dbname)
         if globallogger!= None:
            globallogger.write("WARNING: Trying to create database map file for %s\n" %( dbname) )

         if not doFilesExist( [seqFilePath] ):
            eprintf("ERROR : You do not even have the raw sequence for Database  %s to format!\n", dbname)
            eprintf("      : Make sure you have the file %s\n", seqFilePath)

            if globallogger!= None:
               globallogger.write("ERROR \t You do not even have the raw sequence for Database  %s to format!\n" %( dbname))
               globallogger.write("Make sure you have the file %s\n" %( seqFilePath))

            exit_process()

         mapfile = open(dbmapFile,'w')
         seqFile = open(seqFilePath,'r')
         for line in seqFile:
             if re.match(r'>', line):
                 fprintf(mapfile, "%s\n",line.strip())
         seqFile.close()
         mapfile.close()

    return dbmapFile

# creates the command to blast the sample orf sequences against the reference databas.for the 
# purpose of annotation
def create_blastp_against_refdb_cmd(input, output, output_dir, sample_name,
        dbfilename, config_settings, params,  run_command, algorithm, db_type, globallogger = None): 
    max_evalue = float(get_parameter(params, 'annotation', 'max_evalue', default=0.000001))
    system =    get_parameter(params,  'metapaths_steps', 'BLAST_REFDB', default='yes')
    max_hits =    get_parameter(params,  'annotation', 'max_hits', default=5)


    dbname = get_refdb_name(dbfilename);
    if run_command:
         check_an_format_refdb(dbfilename, 'prot',  config_settings, params, globallogger = globallogger)

    if system=='grid':
       batch_size = get_parameter(params,  'grid_engine', 'batch_size', default=500)
       max_concurrent_batches = get_parameter(params,  'grid_engine', 'max_concurrent_batches', default=500)
       user = get_parameter(params,  'grid_engine', 'user', default='')
       server = get_parameter(params,  'grid_engine', 'server', default='')

       mem = get_parameter(params,  'grid_engine', 'RAM', default='10gb')
       walltime= get_parameter(params,  'grid_engine', 'walltime', default='10:10:10')

       cmd = {'sample_name':output_dir,  'dbnames':[dbname], 'database_files':[dbfilename],\
               'run_type':'overlay', 'batch_size':batch_size,'max_parallel_jobs':max_concurrent_batches,\
               'mem':mem, 'walltime':walltime,\
               'user':user, 'server':server, 'algorithm': algorithm } 
    else:
       if algorithm == 'BLAST':
          cmd="%s -num_threads 16  -max_target_seqs %s  -outfmt 6  -db %s -query  %s -evalue  %f  -out %s"\
            %((config_settings['METAPATHWAYS_PATH'] + config_settings['BLASTP_EXECUTABLE']), max_hits,\
             config_settings['REFDBS'] + PATHDELIM + db_type + PATHDELIM + 'formatted' + PATHDELIM + dbfilename, input, max_evalue, output) 
       if algorithm == 'LAST':
          cmd="%s -o %s -f 0 %s %s"\
            %((config_settings['METAPATHWAYS_PATH'] + config_settings['LAST_EXECUTABLE']), \
              output, config_settings['REFDBS'] + PATHDELIM + db_type + PATHDELIM + 'formatted' + PATHDELIM + dbfilename, input) 
    return cmd


# Command for pars[s.ng the blast flie snd create the parse blast files
#     input -- blastoutput
#     output -- parseed files
#     refscorefile   -- refscore file
#     min_bsr   -- minimum bsr ratio for accepting into annotation
#     max_evalue  -- max evalue
#     min_score   -- min score
#     min_length  -- min_length in amino acids, typcically 100 amino acids.ould be minimum
#     dbmane      -- name of the refrence databases
def create_parse_blast_cmd(input, refscorefile, dbname, dbmapfile,  config_settings, params, algorithm): 
    min_bsr = float(get_parameter(params, 'annotation', 'min_bsr', default=0.4))
    min_score = float(get_parameter(params, 'annotation', 'min_score', default=0.0))
    min_length = float(get_parameter(params, 'annotation', 'min_length', default=100))
    max_evalue = float(get_parameter(params, 'annotation', 'max_evalue', default=1000))

    cmd="%s -d %s  -b %s -m %s  -r  %s  --min_bsr %f  --min_score %f --min_length %f --max_evalue %f"\
         %((config_settings['METAPATHWAYS_PATH'] + config_settings['PARSE_BLAST']), dbname, input, dbmapfile,  refscorefile, min_bsr, min_score, min_length, max_evalue);

    if algorithm == 'LAST':
       cmd = cmd + ' --algorithm LAST'

    if algorithm == 'BLAST':
       cmd = cmd + ' --algorithm BLAST'

    return cmd

# this function creates the command that is required to make the annotated genbank file
def create_annotate_genebank_cmd(sample_name, mapping_txt, input_unannotated_gff, output_annotated_gff,\
           blast_results_dir,  options, refdbs, output_comparative_annotation,\
           config_settings, algorithm):
    cmd="%s --input_gff  %s -o %s  %s --output-comparative-annotation %s \
            --algorithm %s" %((config_settings['METAPATHWAYS_PATH'] +
                config_settings['ANNOTATE']),input_unannotated_gff,\
                output_annotated_gff,  options, output_comparative_annotation,\
                algorithm)

    for refdb in refdbs:
        if algorithm == "LAST":
            cmd = cmd + " -b " + blast_results_dir + PATHDELIM + sample_name + "." + refdb+ "." + algorithm + "out.parsed.txt -d " + refdb + " -w 1 "
        else:
            cmd = cmd + " -b " + blast_results_dir + PATHDELIM + sample_name + "." + refdb+ "." + algorithm + "out.parsed.txt -d " + refdb + " -w 1 "

    cmd = cmd + " -m " +  mapping_txt
    return cmd

def create_genbank_ptinput_sequin_cmd(input_annotated_gff, nucleotide_fasta, amino_fasta, outputs, config_settings, ncbi_params_file, ncbi_sequin_sbt):
    cmd="%s -g %s -n %s -p %s " %((config_settings['METAPATHWAYS_PATH'] + config_settings['GENBANK_FILE']), input_annotated_gff, nucleotide_fasta, amino_fasta) ; 

    if 'gbk' in outputs:
       cmd += ' --out-gbk ' + outputs['gbk']

    if ncbi_params_file and  'sequin' in outputs:
       cmd += ' --out-sequin ' + outputs['sequin']
       cmd += ' --sequin-params-file ' + ncbi_params_file
       cmd += ' --ncbi-sbt-file ' +  ncbi_sequin_sbt

    if 'ptinput' in outputs:
       cmd += ' --out-ptinput ' + outputs['ptinput']


    return cmd

def  create_report_files_cmd(dbs, input_dir, input_annotated_gff,  sample_name, output_dir, config_settings, algorithm, verbose):
    db_argument_string = ''
    for db in dbs:
        dbname = get_refdb_name(db)
        db_argument_string += ' -d ' + dbname
        db_argument_string += ' -b ' + input_dir + sample_name +'.' + dbname
        if algorithm == "LAST":
            db_argument_string += '.' + algorithm + 'out.parsed.txt'
        else:
            db_argument_string += '.' + algorithm + 'out.parsed.txt'
    basefun = config_settings['REFDBS'] + PATHDELIM + 'functional_categories' + PATHDELIM
    basencbi = config_settings['REFDBS'] + PATHDELIM + 'ncbi_tree' + PATHDELIM
    # construct command        
    #cmd="%s %s --input-annotated-gff %s  --input-kegg-maps  %s  --input-cog-maps %s --output-dir %s --ncbi-taxonomy-map %s --seed2ncbi-file %s"\
    cmd="%s %s --input-annotated-gff %s  --input-kegg-maps  %s  --input-cog-maps %s --input-seed-maps %s --output-dir %s --ncbi-taxonomy-map %s "\
           %((config_settings['METAPATHWAYS_PATH'] + config_settings['CREATE_REPORT_FILES']),\
              db_argument_string, input_annotated_gff,\
              basefun + 'KO_classification.txt', basefun + 'COG_categories.txt',  basefun + 'SEED_subsystems.txt',  output_dir,\
               basencbi + 'ncbi_taxonomy_tree.txt')
    if verbose:
       cmd += " -v"
    return cmd


def create_genbank_2_sequin_cmd(input_gbk_file, output_sequin_file, config_settings):
    cmd="%s -i %s -o %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['GENBANK_2_SEQUIN']),input_gbk_file, output_sequin_file); 
    return cmd


def create_reduce_genebank_cmd(input_gbk_file, output_reduced_gbk_file, config_settings):
    cmd="%s -i %s -o %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['REDUCE_ANNOTATED_FILE']),input_gbk_file, output_reduced_gbk_file ); 
    return cmd

#creates the fasta and the pf files to be processed by pathologic and  a subsequent step will 
# create the PGDB
def create_fasta_and_pf_files_cmd(sample_name, input_reduced_annotated_gbk_file_pat, output_fasta_pf_dir, config_settings):
    cmd="%s -s %s -p %s -a %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['PATHOLOGIC_INPUT']), sample_name, output_fasta_pf_dir, input_reduced_annotated_gbk_file_pat)
    return cmd

def create_create_genetic_elements_cmd(output_fasta_pf_dir, config_settings):
    cmd="%s -p %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['CREATE_GENETIC_ELEMENTS']),output_fasta_pf_dir)
    return cmd

# this function creates an organism file for pathway tools
def create_create_organism_params_cmd(yaml_file, output_fasta_pf_dir, pgdb_ID, config_settings):
    cmd="%s --yaml-file %s --ptools-dir %s %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['CREATE_ORGANISM_PARAM']), yaml_file, output_fasta_pf_dir, pgdb_ID)
    return cmd

# this is where the command for PGDB creation is called, by using 
def create_pgdb_using_pathway_tools_cmd(output_fasta_pf_dir, taxonomic_pruning_flag, config_settings):
    cmd="%s -patho %s"  %(config_settings['PATHOLOGIC_EXECUTABLE'], output_fasta_pf_dir +  PATHDELIM)
    if taxonomic_pruning_flag=='no':
        cmd= cmd + " -no-taxonomic-pruning "
        #cmd= cmd + " -no-taxonomic-pruning -no-cel-overview -no-web-cel-ov"
    cmd= cmd + " -no-web-cel-overview"
    return cmd
 

def create_scan_rRNA_seqs_cmd(input_fasta, output_blast, refdb, config_settings,params, command_Status, db_type, globallogger = None):

    if command_Status:
       check_an_format_refdb(refdb, 'nucl',  config_settings, params, globallogger = globallogger)
    else:
       s.exit(0)


    cmd="%s -outfmt 6 -num_threads 16  -query %s -out %s -db %s -max_target_seqs 5"%((config_settings['METAPATHWAYS_PATH'] + config_settings['BLASTN_EXECUTABLE']),input_fasta,output_blast, config_settings['REFDBS'] + PATHDELIM + db_type + PATHDELIM + 'formatted' + PATHDELIM + refdb)
    return cmd


# create the command to do rRNA scan statististics
def create_rRNA_scan_statistics(blastoutput, refdbname, bscore_cutoff, eval_cutoff, identity_cutoff, config_settings, rRNA_stat_results):
    cmd= "%s -o %s -b %s -e %s -s %s"  %((config_settings['METAPATHWAYS_PATH'] + config_settings['SCAN_rRNA']), rRNA_stat_results, bscore_cutoff, eval_cutoff, identity_cutoff)

    cmd =  cmd +  " -i "  + blastoutput + " -d " + config_settings['REFDBS'] + "/taxonomic/" +  refdbname 

    return cmd

# create the command to do tRNA scan statististics
def create_tRNA_scan_statistics(in_file,stat_file, fasta_file,  config_settings):
    cmd= "%s -o %s -F %s  -i %s -T %s  -D %s"  %((config_settings['METAPATHWAYS_PATH'] + config_settings['SCAN_tRNA']) ,\
           stat_file, fasta_file, in_file,\
          (config_settings['METAPATHWAYS_PATH'] + config_settings['RESOURCES_DIR'])+ PATHDELIM + 'TPCsignal',\
          (config_settings['METAPATHWAYS_PATH'] + config_settings['RESOURCES_DIR'])+ PATHDELIM + 'Dsignal')
    return cmd

# create the command to make the MLTreeMap Images
def create_MLTreeMap_Imagemaker(mltreemap_image_output, mltreemap_final_outputs, config_settings):
    executable_path = config_settings['MLTREEMAP_IMAGEMAKER']
    if not path.isfile( executable_path):
        executable_path = config_settings['METAPATHWAYS_PATH'] + executable_path
    cmd= "%s -i %s -o %s -m a"  %(executable_path, mltreemap_final_outputs, mltreemap_image_output) 
    return cmd

# create the command to do mltreemap calculations
def create_MLTreeMap_Calculations(mltreemap_input_file, mltreemap_output_dir, config_settings):
    executable_path = config_settings['MLTREEMAP_CALCULATION']
    if not path.isfile( executable_path):
        executable_path = config_settings['METAPATHWAYS_PATH'] + executable_path
    cmd= "%s -i %s -o %s -e %s"  %(executable_path, mltreemap_input_file, mltreemap_output_dir, config_settings['EXECUTABLES_DIR'] ) 
    return cmd

# gather mltreemap calculations 
def create_MLTreeMap_Hits(mltreemap_output_dir, output_folder, config_settings):
    cmd= "%s -i %s -o %s"  %(config_settings['METAPATHWAYS_PATH'] + config_settings['MLTREEMAP_HITS'], mltreemap_output_dir, output_folder +PATHDELIM + 'sequence_to_cog_hits.txt') 
    return cmd


#gets the parameter value from a category as.ecified in the 
# parameter file
def get_parameter(params, category, field, default = None):
    if params == None:
      return default

    if category in params:
        if field in params[category]:
            return params[category][field]
        else:
            return default
    return default


# parameter file
def get_make_parameter(params,category, field, default = False):
    if category in params:
        if field in params[category]:
            return params[category][field]
        else:
            return default
    return default

def get_pipeline_steps(steps_log_file):
    try:
       logfile = open(steps_log_file, 'r')
    except IOError:
       eprintf("Did not find %s!\n", logfile) 
       eprintf("Try running in \'complete\' run-type\n")
    else:
       lines = logfile.readlines()

    pipeline_steps = None
    return pipeline_steps


def write_run_parameters_file(fileName, parameters):
    try:
       paramFile = open(fileName, 'w')
    except IOError:
       eprintf("Cannot write run parameters to file %s!\n", fileName)
       exit_process("Cannot write run parameters to file %s" %(fileName) )

#       16s_rRNA      {'min_identity': '40', 'max_evalue': '0.000001', 'min_bitscore': '06', 'refdbs': 'silva_104_rep_set,greengenes_db_DW'}
    paramFile.write("\nRun Date : " + str(date.today()) + " \n")

    paramFile.write("\n\nNucleotide Quality Control parameters[s.n")
    paramFile.write( "  min length" + "\t" + str(parameters['quality_control']['min_length']) + "\n")

    paramFile.write("\n\nORF prediction parameters[s.n")
    paramFile.write( "  min length" + "\t" + str(parameters['orf_prediction']['min_length']) + "\n")
    paramFile.write( "  algorithm" + "\t" + str(parameters['orf_prediction']['algorithm']) + "\n")


    paramFile.write("\n\nAmino acid quality control and annotation parameters[s.n")
    paramFile.write( "  min bit score" + "\t" + str(parameters['annotation']['min_score']) + "\n")
    paramFile.write( "  min seq length" + "\t" + str(parameters['annotation']['min_length']) + "\n")
    paramFile.write( "  annotation reference dbs" + "\t" + str(parameters['annotation']['dbs']) + "\n")
    paramFile.write( "  min BSR" + "\t" + str(parameters['annotation']['min_bsr']) + "\n")
    paramFile.write( "  max evalue" + "\t" + str(parameters['annotation']['max_evalue']) + "\n")

    paramFile.write("\n\nPathway Tools parameters[s.n")
    paramFile.write( "  taxonomic pruning " + "\t" + str(parameters['ptools_settings']['taxonomic_pruning']) + "\n")

    paramFile.write("\n\nrRNA search/match parameters[s.n")
    paramFile.write( "  min identity" + "\t" + str(parameters['rRNA']['min_identity']) + "\n")
    paramFile.write( "  max evalue" + "\t" + str(parameters['rRNA']['max_evalue']) + "\n")
    paramFile.write( "  rRNA reference dbs" + "\t" + str(parameters['rRNA']['refdbs']) + "\n")

    paramFile.close()


# checks if the necessary files, directories  and executables really exis.or not
def check_config_settings(config_settings, file, globalerrorlogger = None):
   essentialItems= ['METAPATHWAYS_PATH', 'EXECUTABLES_DIR', 'RESOURCES_DIR']
   missingItems = []

   for key, value in  config_settings.items():
      # make sure  MetaPathways directory is present
      if key in ['METAPATHWAYS_PATH' ]:
         if not path.isdir( config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 1.Currently it is set to \"%s\"\n",  config_settings[key] )  

            if globalerrorlogger!=None:
               globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n"  %(key, file))  
               globalerrorlogger.write("       Currently it is set to \"%s\"\n" %(config_settings[key] )  )
            missingItems.append(key) 
         continue


      # make sure  REFDB directories are present
      if key in [ 'REFDBS' ]:
         if not path.isdir( config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 2.Currently it is set to \"%s\"\n", config_settings[key] )  
            if globalerrorlogger!=None:
                globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key,file))
                globalerrorlogger.write("Currently it is set to \"%s\"\n" %( config_settings[key]) )  
            missingItems.append(key) 
         continue

      # make sure EXECUTABLES_DIR directories are present
      if key in [ 'EXECUTABLES_DIR']:
         if not path.isdir( config_settings['METAPATHWAYS_PATH'] + PATHDELIM +  config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 3.Currently it is set to \"%s\"\n", config_settings[key] )  
            if globalerrorlogger!=None:
               globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key, file))  
               globalerrorlogger.write("Currently it is set to \"%s\"\n" %( config_settings[key] )) 
            missingItems.append(key) 
         continue

      # make sure RESOURCES_DIR directories are present
      if key in [ 'RESOURCES_DIR']:
         if not path.isdir( config_settings['METAPATHWAYS_PATH'] + PATHDELIM + config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 4.Currently it is set to \"%s\"\n",  config_settings['METAPATHWAYS_PATH'] + PATHDELIM + config_settings[key] )  
            print  config_settings['METAPATHWAYS_PATH'], config_settings[key] 
            if globalerrorlogger!=None:
               globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key, file))
               globalerrorlogger.write("Currently it is set to \"%s\"\n" %( config_settings[key]))  
            missingItems.append(key) 
         continue

      # make sure  MetaPaths directory is present
      if key in ['PYTHON_EXECUTABLE' , 'PATHOLOGIC_EXECUTABLE' ]:
         if not path.isfile( config_settings[key]) :
            eprintf("ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
            eprintf("ERROR: 5.Currently it is set to \"%s\"\n", config_settings[key] )  
            if globalerrorlogger!=None:
               globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key, file)) 
               globalerrorlogger.write("Currently it is set to \"%s\"\n" %( config_settings[key] ) )
            missingItems.append(key) 
         continue

      # ignore pgdb folder for now
      if key in ['PGDB_FOLDER' ]:
          continue
      
      # check if the desired file exists. if not, then print a message
      if not path.isfile( config_settings['METAPATHWAYS_PATH'] + PATHDELIM +  value)\
        and  not path.isfile( config_settings['METAPATHWAYS_PATH'] + PATHDELIM + config_settings['EXECUTABLES_DIR'] + PATHDELIM + value ) :
           eprintf("ERROR:Path for \"%s\" is NOT set properly in configuration file \"%s\"\n", key, file)  
           eprintf("5.Currently it is set to \"%s\"\n", config_settings['METAPATHWAYS_PATH'] + value ) 
           if globalerrorlogger!=None:
              globalerrorlogger.write("ERROR\tPath for \"%s\" is NOT set properly in configuration file \"%s\"\n" %(key, file) )
              globalerrorlogger.write("Currently it is set to \"%s\"\n" %(config_settings['METAPATHWAYS_PATH'] + value)) 
           missingItems.append(key) 
           continue
     
   stop_execution = False
   for item in missingItems:
      if item in essentialItems:
         eprintf("ERROR\t Essential field in setting %s is missing in configuration file!\n", item)
         if globalerrorlogger!=None:
            globalerrorlogger.write("ERROR\tEssential field in setting %s is missing in configuration file!\n" %(item))
         stop_execution = True

   if stop_execution ==True:
      eprintf("ERROR: Terminating execution due to missing essential  fields in configuration file!\n")
      if globalerrorlogger!=None:
         globalerrorlogger.write("ERROR\tTerminating execution due to missing essential  fields in configuration file!\n")
      exit_process()

   

# This function reads the pipeline configuration file and sets the 
# paths to differenc scripts and executables the pipeline call
def read_pipeline_configuration( file, globallogger ):
    patternKEYVALUE = re.compile(r'^([^\t\s]+)[\t\s]+\'(.*)\'')
    try:
       configfile = open(file, 'r')
    except IOError:
       eprintf("ERROR :Did not find pipeline config %s!\n", file) 
       globalerrorlogger.write("ERROR\tDid not find pipeline config %s!\n" %(file)) 
    else:
       lines = configfile.readlines()

    config_settings = {}
    for line in lines:
        if not re.match("#",line) and len(line.strip()) > 0 :
           line = line.strip()
           result = patternKEYVALUE.search(line)
           
           try:
              if len(result.groups()) == 2:
                 fields = result.groups()
              else:
                 eprintf("     The following line in your config settings files is not set up yet\n")
                 eprintf("     Please rerun the pipeline after setting up this line\n")
                 eprintf("     Error in line : %s\n", line)
                 globalerrorlogger(
                      "WARNING\t\n"+\
                      "     The following line in your config settings files isn not set up yet\n"+\
                      "     Please rerun the pipeline after setting up this line\n"+\
                      "     Error in line : %s\n" %(line))

                 exit_process()
           except:
                 eprintf("     The following line in your config settings files is not set up yet\n")
                 eprintf("     Please rerun the pipeline after setting up this line\n")
                 eprintf("     Error ine line : %s\n", line)
                 globalerrorlogger(
                      "WARNING\t\n"+\
                      "     The following line in your config settings files is not set up yet\n"+\
                      "     Please rerun the pipeline after setting up this line\n"+\
                      "     Error in line : %s\n" %(line))
                 exit_process()
              
           if PATHDELIM=='\\':
              config_settings[fields[0]] = re.sub(r'/',r'\\',fields[1])   
           else:
              config_settings[fields[0]] = re.sub(r'\\','/',fields[1])   

           
    config_settings['METAPATHWAYS_PATH'] = config_settings['METAPATHWAYS_PATH'] + PATHDELIM
    config_settings['REFDBS'] = config_settings['REFDBS'] + PATHDELIM
    
    check_config_settings(config_settings, file, globallogger);
    config_settings['configuration_file'] = file

    return config_settings

#check for empty values in parameter settings 
def  checkMissingParam_values(params, choices, logger = None):
     reqdCategoryParams = { 
                            'annotation': {'dbs': False}, 
                            'orf_prediction':{}, 
                            'rRNA':{},
                            'metapaths_steps':{}
                         }

     success  = True
     for category in choices:
       for parameter in choices[category]:
         if (not params[category][parameter]) and\
            ( (category in reqdCategoryParams) and\
               (parameter in reqdCategoryParams[category]) and   reqdCategoryParams[category][parameter]) :
            print category, parameter
            print reqdCategoryParams
            print reqdCategoryParams[category]
            eprintf('ERROR: Empty parameter %s of type %s\n'  %(parameter, category))
            eprintf('Please select at least one database for %s\n'  %(category))
            if logger!=None:
               logger.write('ERROR\tEmpty parameter %s of type %s\n'  %(parameter, category))
               logger.write('Please select at least one database for %s\n'  %(category))
            success = False

     return success

# check if all of the metapaths_steps have 
# settings from the valid list [ yes, skip stop, redo]

def  checkParam_values(allcategorychoices, parameters, runlogger = None):
     for category in allcategorychoices:
        for choice in allcategorychoices[category]:
           if choice in parameters: 

             if not parameters[choice] in allcategorychoices[category][choice]:
                 logger.write('ERROR\tIncorrect setting in your parameter file')
                 logger.write('for step %s as %s' %(choice, parameters[choices]))
                 eprintf("ERROR: Incorrect setting in your parameter file" +\
                         "       for step %s as %s", choice, parameters[choices])
                 exit_process()

def checkMetapathsteps(params, runlogger = None):
     choices = { 'metapaths_steps':{}, 'annotation':{}, 'INPUT':{} }

     choices['INPUT']['format']  = ['fasta', 'gbk_unannotated', 'gbk_annotated', 'gff_unannotated', 'gff_annotated']

     choices['annotation']['algorithm'] =  ['last', 'blast'] 

     choices['metapaths_steps']['PREPROCESS_FASTA']   = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['ORF_PREDICTION']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['GFF_TO_AMINO']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['FILTERED_FASTA']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['COMPUTE_REFSCORE']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['BLAST_REFDB'] = ['yes', 'skip', 'stop', 'redo', 'grid']
     choices['metapaths_steps']['PARSE._BLAST'] = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['SCAN_rRNA']   = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['STATS_rRNA']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['ANNOTATE']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['PATHOLOGIC_INPUT']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['GENBANK_FILE']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['CREATE_SEQUIN_FILE']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['CREATE_REPORT_FILES']  = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['SCAN_tRNA']   = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['MLTREEMAP_CALCULATION']   = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['MLTREEMAP_IMAGEMAKER']    = ['yes', 'skip', 'stop', 'redo']
     choices['metapaths_steps']['PATHOLOGIC']  = ['yes', 'skip', 'stop', 'redo']


     if params['metapaths_steps']:
        checkParam_values(choices, params['metapaths_steps'], runlogger)

     checkparams = {}
     checkparams['annotation'] = []
     checkparams['annotation'].append('dbs') 

     if not checkMissingParam_values(params, checkparams, runlogger):
        exit_process("Missing parameters")


def  copy_fna_faa_gff_orf_prediction( source_files, target_files, config_settings) :

     for source, target in zip(source_files, target_files):  

         sourcefile = open(source, 'r')
         targetfile = open(target, 'w')
         sourcelines = sourcefile.readlines()
         for line in sourcelines:
            fprintf(targetfile, "%s\n", line.strip())

         sourcefile.close()
         targetfile.close()


#################################################################################
########################### BEFORE BLAST ########################################
#################################################################################
def run_metapathways(samplesData, output_dir, all_samples_output_dir, globallogger,\
                     command_line_params, params, metapaths_config, status_update_callback,\
                     config_file, run_type, config_settings = None, block_mode = False, runid = ""):

    jobcreator = JobCreator(params, config_settings)

    for input_file in samplesData.keys():
      s =  samplesData[input_file]
      jobcreator.addJobs(s, block_mode = block_mode)



    if block_mode:
       for input_file in samplesData.keys():
         s =  samplesData[input_file]
         s.stepslogger.printf("\n\n==============  BEGIN RUN " + s.sample_name + " " + runid + " ================\n")
         sample_name_banner = "PROCESSING INPUT " + input_file
         eprintf("==============  BEGIN RUN " + s.sample_name + " STEPS BLOCK 0 ================\n")
         eprintf('#'*len(sample_name_banner) + "\n")
         eprintf( '\n' + sample_name_banner + '\n')
         execute_tasks(s, verbose = command_line_params['verbose'], block = 0)    


       for input_file in samplesData.keys():
         s =  samplesData[input_file]
         sample_name_banner = "PROCESSING INPUT " + input_file
         eprintf("==============  BEGIN RUN " + s.sample_name + " STEPS BLOCK 1 ================\n")
         eprintf('#'*len(sample_name_banner) + "\n")
         eprintf( '\n' + sample_name_banner + '\n')
         execute_tasks(s, verbose = command_line_params['verbose'], block = 1)    

       for input_file in samplesData.keys():
         s =  samplesData[input_file]
         sample_name_banner = "PROCESSING INPUT " + input_file
         eprintf("==============  BEGIN RUN " + s.sample_name + " STEPS BLOCK 2 ================\n")
         eprintf('#'*len(sample_name_banner) + "\n")
         eprintf( '\n' + sample_name_banner + '\n')
         execute_tasks(s, verbose = command_line_params['verbose'], block = 2)    

    else:
       for input_file in samplesData.keys():
         s =  samplesData[input_file]
         s.stepslogger.printf("\n\n==============  BEGIN RUN " + s.sample_name + " " + runid + "  ==================\n")
         sample_name_banner = "PROCESSING INPUT " + input_file
         eprintf('#'*len(sample_name_banner) + "\n")
         eprintf( '\n' + sample_name_banner + '\n')
         execute_tasks(s, verbose = command_line_params['verbose'], block = 0)    

    return


