#!/usr/bin/python

import re
import sys

class LCAComputation:
    begin_pattern = re.compile("#")

    name_to_id={}
    id_to_name={}
    taxid_to_ptaxid = {}
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
          self.taxid_to_ptaxid[str(fields[1])] = [ str(fields[2]), 0]


    def sizeTaxnames(self ):
         return len(self.name_to_id)

    def sizeTaxids(self):
         return len(self.taxid_to_ptaxid)
          
    def get_a_Valid_ID(self, name_group):
        for name in name_group:
           if name in self.name_to_id:
               return  self.name_to_id[name]
        return -1

    def translateNameToID(self, name):
       if not name in self.name_to_id:
           return None
       return self.name_to_id[name]

    def translateIdToName(self, id):
       if not id in self.id_to_name:
           return None
       return self.id_to_name[id]


    def getParentName(self, name):
       if not name in  self.name_to_id:  
          return None
       id = self.name_to_id[name]  
       pid = self.getParentTaxId(id)
       return self.translateIdToName( pid )


    def getParentTaxId(self, ID):
       if not ID in self.taxid_to_ptaxid:
          return None
       return self.taxid_to_ptaxid[ID][0]


    def get_lca(self, IDs):
        limit = len(IDs)
        for id in IDs:
           tid = id 
           while( tid in self.taxid_to_ptaxid and tid !='1' ):
               self.taxid_to_ptaxid[tid][1]+=1
               if self.taxid_to_ptaxid[tid][1]==limit:
                  return  self.id_to_name[tid]  
               tid = self.taxid_to_ptaxid[tid][0]
        return ""

    def clear_cells(self, IDs):
        limit = len(IDs)
        for id in IDs:
           tid = id 
           while( tid in self.taxid_to_ptaxid and tid !='1' ):
               if self.taxid_to_ptaxid[tid][1]==0:
                  return  self.id_to_name[tid]  
               self.taxid_to_ptaxid[tid][1]=0
               tid = self.taxid_to_ptaxid[tid][0]
        return ""


    def getTaxonomy(self, name_groups):

         IDs = []
         for name_group in name_groups:
            #print name_group
            id = self.get_a_Valid_ID(name_group)
            if id!=-1:
              IDs.append(id)
    
         consensus = self.get_lca(IDs)
         #print "=============="
         self.clear_cells(IDs)
         return consensus

