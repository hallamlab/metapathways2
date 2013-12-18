#!/usr/bin/env python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


"""Contains general utility code for the metapaths project"""


from shutil import rmtree
from StringIO import StringIO
from os import getenv, makedirs
from operator import itemgetter
from os.path import abspath, exists, dirname, join, isdir
from collections import defaultdict
from optparse import make_option
import sys 
import os
import math
#from copy import deepcopy
#from numpy import min, max, median, mean
#import numpy
#from numpy.ma import MaskedArray
#from numpy.ma.extras import apply_along_axis
#from numpy import array, zeros, argsort, shape, vstack,ndarray, asarray, \
#        float, where, isnan
#from cogent import LoadSeqs, Sequence
#from cogent.cluster.procrustes import procrustes
#from cogent.core.alignment import Alignment
#from cogent.core.moltype import MolType, IUPAC_DNA_chars, IUPAC_DNA_ambiguities,\
#    IUPAC_DNA_ambiguities_complements, DnaStandardPairs, ModelDnaSequence
#from cogent.data.molecular_weight import DnaMW
#from cogent.core.sequence import DnaSequence
#from cogent.app.blast import Blastall
#from cogent.app.util import get_tmp_filename
#from cogent.parse.blast import BlastResult
#from cogent.parse.fasta import MinimalFastaParser
#from cogent.util.misc import remove_files
#from cogent.util.dict2d import Dict2D
#from cogent.app.formatdb import build_blast_db_from_fasta_path,\
#    build_blast_db_from_fasta_file
#from cogent import LoadSeqs
#from cogent.util.misc import (parse_command_line_parameters, 
#                                     create_dir, 
#                                     handle_error_codes)
#
#
#from qiime.parse import (parse_otu_table,
#                         parse_qiime_config_files,
#                         parse_coords,
#                         parse_newick,
#                         PhyloNode,
#                         parse_mapping_file)
#from qiime.format import format_otu_table
#

class Performance:
   def  __init__(self):
       self.sum = {}
       self.sqsum = {}
       self.num = {}

   def getAverageDelay(self, server= None):
      if server==None:
         avg = 0
         num = 0
         for server in self.sum:
            avg += self.sum[server]
            num += self.num[server]
         if num > 0: 
            return avg/num
         else:
            return 0
         
      if self.num[server]==0:
         return 0
      avg = self.sum[server]/self.num[server]
      return avg
      
   def getStdDeviationDelay(self, server= None):
      if server==None:
         avg = 0
         avgsq = 0
         num = 0
         for server in self.sum:
            avg +=  self.sum[server]
            avgsq += self.sqsum[server]
            num += self.num[server]
         if num == 0: 
            return 0
         
      var = avgsq/num - avg*avg/(num*num)
      std  = math.sqrt(var)
      return std
      



   def addPerformanceData(self, server, data):
      if not server in self.sum: 
         self.sum[server] = 0
         self.sqsum[server] = 0
         self.num[server] = 0

      self.sum[server] += data
      self.sqsum[server] += data*data
      self.num[server] += 1
      return True


   def getExpectedDelay(self):
      return 20


class Job:
   def  __init__(self, S, d, a, m, server=None):
      self.S = S  # sample
      self.d = d  # database
      self.a = a  # split
      self.m = m  # algorithm
      self.server = None  # server
      return None

   def setValues(self, S, d, a, m, t, server=None):
      self.S = S
      self.d = d
      self.a = a
      self.m = m
      self.submission_time = t
      self.server=server
      return True


def parse_command_line_parameters(script_info, argv):
    print script_info 
    print argv
    opts = []
    return opts
class TreeMissingError(IOError):
    """Exception for missing tree file"""
    pass

class OtuMissingError(IOError):
    """Exception for missing OTU file"""
    pass

class AlignmentMissingError(IOError):
    """Exception for missing alignment file"""
    pass

class MissingFileError(IOError):
    pass

def make_safe_f(f, allowed_params):
    """Make version of f that ignores extra named params."""
    def inner(*args, **kwargs):
        if kwargs:
            new_kwargs = {}
            for k, v in kwargs.items():
                if k in allowed_params:
                    new_kwargs[k] = v
            return f(*args, **new_kwargs)
        return f(*args, **kwargs)
    return inner

class FunctionWithParams(object):
    """A FunctionWithParams is a replacement for the function factory.
    
    Specifically, the params that will be used in the __call__ method are
    available in a dict so you can keep track of them with the object
    itself.
    """
    Application = None
    Algorithm = None
    Citation = None
    Params = {}
    Name = 'FunctionWithParams' #override in subclasses
    _tracked_properties = []    #properties tracked like params

    def __init__(self, params):
        """Return new FunctionWithParams object with specified params.
        
        Note: expect params to contain both generic and per-method (e.g. for
        cdhit) params, so leaving it as a dict rather than setting
        attributes. 
        
        Some standard entries in params are:

        [fill in on a per-application basis]
        """
        self.Params.update(params)
        self._tracked_properties.extend(['Application','Algorithm','Citation'])

    def __str__(self):
        """Returns formatted key-value pairs from params."""
        res = [self.Name + ' parameters:']
        for t in self._tracked_properties:
            res.append(t + ':' + str(getattr(self, t)))
        for k, v in sorted(self.Params.items()):
            res.append(str(k) + ':' + str(v))
        return '\n'.join(res)

    def writeLog(self, log_path):
        """Writes self.Params and other relevant info to supplied path."""
        f=open(log_path, 'w')
        f.write(str(self))
        f.close()

    def getResult(self, *args, **kwargs):
        """Gets result in __call__. Override in subclasses."""
        return None

    def formatResult(self, result):
        """Formats result as string (for whatever "result" means)."""
        return str(result)

    def writeResult(self, result_path, result):
        """Writes result to result_path. May need to format in subclasses."""
        f = open(result_path, 'w')
        f.write(self.formatResult(result))
        f.close()

    def getOtuTable(self, otu_source):
        """Returns parsed OTU table from putative OTU source."""
        
        #if we have a string starting with #, assume it's an OTU file,
        #otherwise assume it's a path
        # if 4-tuple, just return it
        if type(otu_source) == type((1,3,4,44)):
            return otu_source
        if hasattr(otu_source, 'startswith') and otu_source.startswith('#'):
            try:
                return parse_otu_table(StringIO(otu_source))
            except (TypeError, ValueError), e:
                raise OtuMissingError, \
                    "Tried to read OTUs from string starting with # but got "+e
        else:
            try:
                otu_file = open(otu_source, 'U')
            except (TypeError, IOError):
                raise OtuMissingError, \
                    "Couldn't read OTU file at path: %s" % otu_source
            result = parse_otu_table(otu_file)
            otu_file.close()
            return result

    def getTree(self, tree_source):
        """Returns parsed tree from putative tree source"""
        if isinstance(tree_source, PhyloNode):
            tree = tree_source    #accept tree object directly for tests
        elif tree_source:
            try:
                f = open(tree_source, 'U')
            except (TypeError, IOError):
                raise TreeMissingError, \
                    "Couldn't read tree file at path: %s" % tree_source
            tree = parse_newick(f, PhyloNode)
            f.close()
        else:
            raise TreeMissingError, str(self.Name) + \
                " is a phylogenetic metric, but no tree was supplied."
        return tree

    def getData(self, data_source):
        """Returns data from putative source, which could be a path"""
        if isinstance(data_source, str):
            try:
                return eval(data_source)
            except (NameError, SyntaxError):
                try:
                    data_f = open(data_source, 'U')
                    data = data_f.read()
                    data_f.close()
                    try:
                        return eval(data)
                    except (NameError, SyntaxError, TypeError):
                        pass
                    return data
                except (IOError, NameError, TypeError):
                    pass
        #if we got here, either we didn't get a string or we couldn't read
        #the data source into any other kind of object
        return data_source

    def getAlignment(self, aln_source):
        """Returns parsed alignment from putative alignment source"""
        if isinstance(aln_source, Alignment):
            aln = aln_source
        elif aln_source:
            try:
                aln = LoadSeqs(aln_source, Aligned=True)
            except (TypeError, IOError, AssertionError):
                raise AlignmentMissingError, \
                    "Couldn't read alignment file at path: %s" % aln_source
        else:
            raise AlignmentMissingError, str(self.Name) + \
                " requires an alignment, but no alignment was supplied."
        return aln

    def __call__ (self, result_path=None, log_path=None,\
        *args, **kwargs):
        """Returns the result of calling the function using the params dict.
        
        Parameters:
        [fill in on a per-application basis]
        """
        print """Function with parameters"""
        result = self.getResult(*args, **kwargs)
        if log_path:
            self.writeLog(log_path)
        if result_path:
            self.writeResult(result_path, result)
        else:
            return result

def get_qiime_project_dir():
    """ Returns the top-level QIIME directory
    """
    # Get the full path of util.py
    current_file_path = abspath(__file__)
    # Get the directory containing util.py
    current_dir_path = dirname(current_file_path)
    # Return the directory containing the directory containing util.py
    return dirname(current_dir_path)
    
def get_qiime_scripts_dir():
    """ Returns the QIIME scripts directory 
    
        This value must be stored in qiime_config if the user
        has installed qiime using setup.py. If it is not in
        qiime_config, it is inferred from the qiime_project_dir.
    
    """
    qiime_config = load_qiime_config()
    qiime_config_value = qiime_config['qiime_scripts_dir']
    if qiime_config_value != None:
        result = qiime_config_value
    else:
        result = join(get_qiime_project_dir(),'scripts')
    
    #assert exists(result),\
    # "qiime_scripts_dir does not exist: %s." % result +\
    # " Have you defined it correctly in your qiime_config?"
    
    return result
    
def load_qiime_config():
    """Return default parameters read in from file"""
    
    qiime_config_filepaths = []
    qiime_project_dir = get_qiime_project_dir()
    qiime_config_filepaths.append(\
     qiime_project_dir + '/qiime/support_files/qiime_config')
    
    qiime_config_env_filepath = getenv('QIIME_CONFIG_FP')
    if qiime_config_env_filepath:
        qiime_config_filepaths.append(qiime_config_env_filepath)
    
    home_dir = getenv('HOME')
    if home_dir:
        qiime_config_home_filepath = home_dir + '/.qiime_config'
        qiime_config_filepaths.append(qiime_config_home_filepath)
    
    qiime_config_files = []
    for qiime_config_filepath in qiime_config_filepaths:
        if exists(qiime_config_filepath):
            qiime_config_files.append(open(qiime_config_filepath))
    
    return parse_qiime_config_files(qiime_config_files)

# The qiime_blast_seqs function should evetually move to PyCogent,
# but I want to test that it works for all of the QIIME functionality that
# I need first. -Greg

def extract_seqs_by_sample_id(seqs, sample_ids, negate=False):
    """ Returns (seq id, seq) pairs if sample_id is in sample_ids """
    sample_ids = {}.fromkeys(sample_ids)

    if not negate:
        def f(s):
            return s in sample_ids
    else:
        def f(s):
            return s not in sample_ids

    for seq_id, seq in seqs:
        sample_id = seq_id.split('_')[0]
        if f(sample_id):
            yield seq_id, seq
            
def split_fasta_on_sample_ids(seqs):
    """ yields (sample_id, seq_id, seq) for each entry in seqs 
    
        seqs: (seq_id,seq) pairs, as generated by MinimalFastaParser
    
    """
    for seq_id, seq in seqs:
        yield (seq_id.split()[0].rsplit('_',1)[0], seq_id, seq)
    return
        
def split_fasta_on_sample_ids_to_dict(seqs):
    """ return split_fasta_on_sample_ids as {sample_id: [(seq_id, seq), ], }
    
        seqs: (seq_id,seq) pairs, as generated by MinimalFastaParser
    
    """
    result = {}
    for sample_id,seq_id,seq in split_fasta_on_sample_ids(seqs):
        try:
            result[sample_id].append((seq_id,seq))
        except KeyError:
            result[sample_id] = [(seq_id,seq)]
    return result
    
def split_fasta_on_sample_ids_to_files(seqs,output_dir):
    """ output of split_fasta_on_sample_ids to fasta in specified output_dir
    
        seqs: (seq_id,seq) pairs, as generated by MinimalFastaParser
        output_dir: string defining directory where output should be 
         written, will be created if it doesn't exist
    
    """
    create_dir(output_dir)
    file_lookup = {}
    for sample_id,seq_id,seq in split_fasta_on_sample_ids(seqs):
        try:
            file_lookup[sample_id].write('>%s\n%s\n' % (seq_id,seq))
        except KeyError:
            file_lookup[sample_id] = open('%s/%s.fasta' % 
                                          (output_dir,sample_id),'w')
            file_lookup[sample_id].write('>%s\n%s\n' % (seq_id,seq))
    for file_handle in file_lookup.values():
        file_handle.close()
    return None

     
def raise_error_on_parallel_unavailable(qiime_config=None):
    """Raise error if no parallel QIIME bc user hasn't set jobs_to_start
    """
    if qiime_config == None:
        qiime_config = load_qiime_config()
    if 'jobs_to_start' not in qiime_config or \
       int(qiime_config['jobs_to_start']) < 2:
       raise RuntimeError,\
        "Parallel QIIME is not available. (Have you set"+\
        " jobs_to_start to greater than 1 in your qiime_config?"
        
def matrix_stats(headers_list, distmats):
    """does, mean, median, stdev on a series of (dis)similarity matrices
    
    takes a series of parsed matrices (list of headers, list of numpy 2d arrays)
    headers must are either row or colunm headers (those must be identical)
    outputs headers (list), means, medians, stdevs (all numpy 2d arrays)
    """
    
    if len(set(map(tuple,headers_list))) > 1:
        raise ValueError("error, not all input matrices have"+\
          " identical column/row headers")
        
    all_mats = numpy.array(distmats) # 3d numpy array: mtx, row, col
    means = numpy.mean(all_mats, axis=0)
    medians = numpy.median(all_mats, axis=0)
    stdevs = numpy.std(all_mats, axis=0)
    
    return deepcopy(headers_list[0]), means, medians, stdevs

def _flip_vectors(jn_matrix, m_matrix):
    """transforms PCA vectors so that signs are correct"""
    m_matrix_trans = m_matrix.transpose()
    jn_matrix_trans = jn_matrix.transpose()
    new_matrix= zeros(jn_matrix_trans.shape, float)
    for i, m_vector in enumerate(m_matrix_trans):
        jn_vector = jn_matrix_trans[i]
        disT = list(m_vector - jn_vector)
        disT = sum(map(abs, disT))
        jn_flip = jn_vector*[-1]
        disF = list(m_vector - jn_flip)
        disF = sum(map(abs, disF))
        if disT > disF:
            new_matrix[i] = jn_flip
        else:
            new_matrix[i] = jn_vector
    return new_matrix.transpose()

def IQR(x):
    """calculates the interquartile range of x

    x can be a list or an array
    
    returns min_val and  max_val of the IQR"""

    x.sort()
    #split values into lower and upper portions at the median
    odd = len(x) % 2
    midpoint = int(len(x)/2)
    if odd:
        low_vals = x[:midpoint]
        high_vals = x[midpoint+1:]
    else: #if even
        low_vals = x[:midpoint]
        high_vals = x[midpoint:]
    #find the median of the low and high values
    min_val = median(low_vals)
    max_val = median(high_vals)
    return min_val, max_val

def matrix_IQR(x):
    """calculates the IQR for each column in an array
    """
    num_cols = x.shape[1]
    min_vals = zeros(num_cols)
    max_vals = zeros(num_cols)
    for i in range(x.shape[1]):
        col = x[:, i]
        min_vals[i], max_vals[i] = IQR(col)
    return min_vals, max_vals

def idealfourths(data, axis=None):
    """This function returns an estimate of the lower and upper quartiles of the data along
    the given axis, as computed with the ideal fourths. This function was taken
    from scipy.stats.mstat_extra.py (http://projects.scipy.org/scipy/browser/trunk/scipy/stats/mstats_extras.py?rev=6392)
    """
    def _idf(data):
        x = data.compressed()
        n = len(x)
        if n < 3:
            return [numpy.nan,numpy.nan]
        (j,h) = divmod(n/4. + 5/12.,1)
        qlo = (1-h)*x[j-1] + h*x[j]
        k = n - j
        qup = (1-h)*x[k] + h*x[k-1]
        return [qlo, qup]
    data = numpy.sort(data, axis=axis).view(MaskedArray)
    if (axis is None):
        return _idf(data)
    else:
        return apply_along_axis(_idf, axis, data)

def isarray(a):
    """
    This function tests whether an object is an array
    """
    try:
        validity=isinstance(a,ndarray)
    except:
        validity=False

    return validity


def degap_fasta_aln(seqs):
    """degap a Fasta aligment.

    seqs: list of label,seq pairs
    """
    
    for (label,seq) in seqs:
        degapped_seq = Sequence(moltype=DNA_with_more_gaps,
                                seq=seq, name=label).degap()
        degapped_seq.Name = label
        yield degapped_seq

def write_degapped_fasta_to_file(seqs, tmp_dir="/tmp/"):
    """ write degapped seqs to temp fasta file."""

    tmp_filename = get_tmp_filename(tmp_dir=tmp_dir, prefix="degapped_", suffix=".fasta")
    fh = open(tmp_filename,"w")
    
    for seq in degap_fasta_aln(seqs):
        fh.write(seq.toFasta()+"\n")
    fh.close()
    return tmp_filename


def fprintf(file, fmt, *args): 
    file.write(fmt % args)


def printf(fmt, *args): 
    sys.stdout.write(fmt % args)


def eprintf(fmt, *args): 
    sys.stderr.write(fmt % args)


# remove the string "/pathway-tools" to infer the pathway tools dir
def create_pathway_tools_dir_path_From_executable(pathway_tools_executable):
    return( pathway_tools_executable.replace('pathway-tools/pathway-tools', 'pathway-tools'))


#removes an existing pgdb from the  ptools-local/pgdbs/user directory under the 
#pathway tools directory
def remove_existing_pgdb( sample_name, pathway_tools_exec):
   suffix_to_remove = ""
   # crete the pathway tools dir
   pathway_tools_dir = create_pathway_tools_dir_path_From_executable(pathway_tools_exec)

   sample_pgdb_dir = pathway_tools_dir + "/" + "ptools-local/pgdbs/user/" + sample_name + "cyc"
   if os.path.exists(sample_pgdb_dir):
      return rmtree(sample_pgdb_dir)
    
