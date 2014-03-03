#!/usr/bin/python

try:
     import sys, traceback
     import re
     import sys
     from   metapaths_utils  import fprintf, printf, GffFileParser
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     print traceback.print_exc(10)
     sys.exit(3)



def copyList(a, b): 
    [ b.append(x) for x in a ] 

class LCAComputation:
    begin_pattern = re.compile("#")

    # a readable taxon name to numeric string id map as ncbi
    name_to_id={}

    # a readable taxon ncbi tax id to name map
    id_to_name={}

    # this is the tree structure in a id to parent map, you can traverse it to go to the root
    taxid_to_ptaxid = {}

    lca_min_score = 50   # an LCA parameter for min score for a hit to be considered
    lca_top_percent = 10    # an LCA param to confine the hits to within the top hits score upto the top_percent% 
    lca_min_support = 5   # a minimum number of reads in the sample to consider a taxon to be present
    results_dictionary = None

    # initialize with the ncbi tree file 
    def __init__(self, filename):
       taxonomy_file = open(filename, 'r')
       lines = taxonomy_file.readlines()
       taxonomy_file.close()

       for line in lines:
          if self.begin_pattern.search(line):
              continue
          fields =  [ x.strip()  for x in line.rstrip().split('\t')]
          if len(fields) !=3:
              continue
          self.name_to_id[str(fields[0])] = str(fields[1])
          self.id_to_name[str(fields[1])] = str(fields[0])
          # the taxid to ptax map has for each taxid a corresponding 3-tuple
          # the first location is the pid, the second is used as a counter for 
          # lca while a search is traversed up the tree and the third is used for
          # the min support
          self.taxid_to_ptaxid[str(fields[1])] = [ str(fields[2]), 0, 0]


    def setParameters(self, min_score, top_percent, min_support):
       self.lca_min_score = min_score
       self.lca_top_percent =top_percent
       self.lca_min_support = min_support
         
    def sizeTaxnames(self ):
         return len(self.name_to_id)


    def sizeTaxids(self):
         return len(self.taxid_to_ptaxid)
          
    def get_a_Valid_ID(self, name_group):
        for name in name_group:
           if name in self.name_to_id:
               return  self.name_to_id[name]
        return -1

    # given a taxon name it returns the correcponding unique ncbi tax id
    def translateNameToID(self, name):
       if not name in self.name_to_id:
           return None
       return self.name_to_id[name]

    # given a taxon id to taxon name map
    def translateIdToName(self, id):
       if not id in self.id_to_name:
           return None
       return self.id_to_name[id]


    # given a name it returns the parents name
    def getParentName(self, name):
       if not name in  self.name_to_id:  
          return None
       id = self.name_to_id[name]  
       pid = self.getParentTaxId(id)
       return self.translateIdToName( pid )


    # given a ncbi tax id returns the parents tax id
    def getParentTaxId(self, ID):
       if not ID in self.taxid_to_ptaxid:
          return None
       return self.taxid_to_ptaxid[ID][0]


    # given a set of ids it returns the lowest common ancenstor 
    # without caring about min support
    # here LCA for a set of ids are computed as follows
    # first we consider one ID at a time
    #   for each id we traverse up the ncbi tree using the id to parent id map
    #   at the same time increasing the count on the second value of the 3-tuple 
    #   note that at the node where all the of the individual ids ( limit in number)
    #   converges the counter matches the limit for the first time, while climbing up. 
    #   This also this enables us to  make the selection of id arbitrary 
    def get_lca(self, IDs):
        limit = len(IDs)
        for id in IDs:
           tid = id 
           while( tid in self.taxid_to_ptaxid and tid !='1' ):
               self.taxid_to_ptaxid[tid][1]+=1
               if self.taxid_to_ptaxid[tid][1]==limit:
                  return  self.id_to_name[tid]  
               tid = self.taxid_to_ptaxid[tid][0]

        return "root"

    def update_taxon_support_count(self, taxonomy):
         id = self.get_a_Valid_ID( [taxonomy ])
         tid = id 
         while( tid in self.taxid_to_ptaxid and tid !='1' ):
               self.taxid_to_ptaxid[tid][2]+=1
               tid = self.taxid_to_ptaxid[tid][0]

    def get_supported_taxon(self, taxonomy):
         id = self.get_a_Valid_ID( [taxonomy ])
         tid = id 
         while( tid in self.taxid_to_ptaxid and tid !='1' ):
            if self.lca_min_support > self.taxid_to_ptaxid[tid][2] :
                tid = self.taxid_to_ptaxid[tid][0]
            else:
                return self.translateIdToName(tid)

         return  self.translateIdToName(tid)
    
    # need to call this to clear the counts of reads at every node      
    def clear_cells(self, IDs):
        limit = len(IDs)
        for id in IDs:
           tid = id 
           while( tid in self.taxid_to_ptaxid and tid !='1' ):
               #if self.taxid_to_ptaxid[tid][1]==0:
               #   return  self.id_to_name[tid]  
               self.taxid_to_ptaxid[tid][1]=0
               tid = self.taxid_to_ptaxid[tid][0]
        return ""


    #given a set of sets of names it computes an lca 
    # in the format [ [name1, name2], [name3, name4,....namex] ...]
    # here name1 and name2 are synonyms and so are name3 through namex
    def getTaxonomy(self, name_groups):
         IDs = []
         for name_group in name_groups:
            id = self.get_a_Valid_ID(name_group)
            if id!=-1:
              IDs.append(id)
    
         consensus = self.get_lca(IDs)
         self.clear_cells(IDs)
         return consensus


    # extracts taxon names for a refseq annotation
    def get_species(self, hit):
       if not 'product' in hit: 
           return None
       species = []
       try:
           m = re.findall(r'\[([^\[]+)\]', hit['product'])
           if m != None:
             copyList(m,species)
       except:
             return None
   
       if species:
          return species
       else:
          return None
 
    # used for optimization
    def set_results_dictionary(self, results_dictionary):
        self.results_dictionary= results_dictionary

    # this returns the megan taxonomy, i.e., it computes the lca but at the same time
    # takes into consideration the parameters, min score, min support and top percent
    def getMeganTaxonomy(self, orfid):
         #compute the top hit wrt score
         names = []
         species = []
         if 'refseq' in self.results_dictionary:
            if orfid in self.results_dictionary['refseq']:
                 
               top_score = 0 
               for hit in self.results_dictionary['refseq'][orfid]:
                  if hit['bitscore'] >= self.lca_min_score and hit['bitscore'] >= top_score:
                     top_score = hit['bitscore']

               for hit in self.results_dictionary['refseq'][orfid]:
                  if (100-self.lca_top_percent)*top_score/100 < hit['bitscore']:
                     names = self.get_species(hit)
                     #if 'MD_2_95' == orfid:
                     #  for hit in self.results_dictionary['refseq'][orfid]:
                     #     print  orfid  + ':' + str(names)
                     #  else:
                     #     print orfid  + ':' + str([])
                     if names:
                       species.append(names) 

         taxonomy = self.getTaxonomy(species)
         meganTaxonomy = self.get_supported_taxon( taxonomy)
         return meganTaxonomy
 

    # this is use to compute the min support for each taxon in the tree
    # this is called before the  getMeganTaxonomy
    def compute_min_support_tree(self, annotate_gff_file, pickorfs):
        gffreader = GffFileParser(annotate_gff_file)
        try:
           for contig in  gffreader:
              for orf in  gffreader.orf_dictionary[contig]:
                 if not orf['id'] in pickorfs:
                     continue
                 taxonomy = None
                 species = []
                 if 'refseq' in self.results_dictionary:
                   if orf['id'] in self.results_dictionary['refseq']:
                       #compute the top hit wrt score
                       top_score = 0 
                       for hit in self.results_dictionary['refseq'][orf['id']]:
                          if hit['bitscore'] >= self.lca_min_score and hit['bitscore'] >= top_score:
                            top_score = hit['bitscore']
       
                       for hit in self.results_dictionary['refseq'][orf['id']]:
                          if (100-self.lca_top_percent)*top_score/100 < hit['bitscore']:
                             names = self.get_species(hit)
                             if names:
                               species.append(names) 
                 taxonomy=self.getTaxonomy(species)
                 self.update_taxon_support_count(taxonomy)
                 pickorfs[orf['id']] = taxonomy
        except:
           print "ERROR : Cannot read annotated gff file "
          

