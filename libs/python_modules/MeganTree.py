#!/usr/bin/python

try:
    #from python_modules.LCAComputation import *
    from LCAComputation import *
    from os import sys
except:
    print """ Could not load some user defined  module functions"""
    sys.exit(3)


class MeganTree:
    begin_pattern = re.compile("#")

    child_to_parent={}
    parent_to_child={}
    id_to_name={}
    taxid_to_ptaxid = {}
    children_at_leaf = {}
    lca = None
    output=None
    def __init__(self, lca):
       self.lca = lca

    def insertTaxon(self, name):
       oid = self.lca.translateNameToID(name)
       #print "insert"
       id = oid
       while id != '1':
          pid = self.lca.getParentTaxId(id)
          #print name + " " + pid
          if  pid is None:
             return
          if not id in self.child_to_parent:
             self.child_to_parent[id] = [pid, 0]
          self.child_to_parent[id][1] += 1
          id = pid

       if not oid in self.children_at_leaf:
         self.children_at_leaf[oid] = 0
       self.children_at_leaf[oid] +=1

    def getChildToParentMap(self):
       return self.child_to_parent

    def getParentToChildrenMap(self):
        if not self.parent_to_child:
          self.computeParentToChildrenMap()

        return self.parent_to_child

    def computeParentToChildrenMap(self):
        for id, pid in self.child_to_parent.iteritems(): 
            if not pid[0] in self.parent_to_child:
                 self.parent_to_child[pid[0]] = [ [],  pid[1] ]
            self.parent_to_child[pid[0]][0].append(id)


    def printTree(self, id):
        if not self.parent_to_child:
          self.computeParentToChildrenMap()

        self.output=""
        self.createTree(id)
        return self.output
        
    def createTree(self, id):
        if not id in self.parent_to_child: 
           if id=='1':
              print self.parent_to_child
              sys.exit(0)
           self.output +=  str(self.lca.translateIdToName(id)) + "}" +  str(self.children_at_leaf[id])
           #self.output +=  id + ":" +  str(0)
           return 

        children =  self.parent_to_child[id][0] 
        count =  self.parent_to_child[id][1] 
        
        self.output += "("
        i =0
        for child in children:
            if i > 0 :
              self.output += "{"
            self.createTree(child)
            i+=1

        self.output += ")" +  str(self.lca.translateIdToName(id)) + "}" +  str(count)
        ##self.output += ")" +  id + ":" +  str(count)

