#!/usr/bin/python -w

import socket, re, sys
from os import _exit 
import re
from libs.python_modules.utils.sysutil import getstatusoutput
from multiprocessing import Process
import time

class PythonCyc:
    
      _organism = 'meta' 
      soc = None
       
      
      def __init__(self):
          pass


      def setOrganism(self, organism):
          self._organism = organism

      def makeSocket(self):
          self.soc = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
          self.soc.connect("/tmp/ptools-socket" )



      def tokenize(self,string):
          LPAREN = '\(';
          RPAREN = '\)';
          WSPACE = '\s+';
          STRING = '\"[^"]*?\"';
          PIPES = '\|[^\|]*?\|';

          regexp = re.compile(r'' + "(" + LPAREN + "|" + RPAREN + "|" + STRING + "|"+PIPES+ ")" + "|" + WSPACE)

          _tokens = [ x.strip() for x in   regexp.split(string) if x ]


          removePipes = re.compile(r'\|([^|]*)\|')
          tokens = [] 
          for _token in _tokens:
              tokens.append(re.sub(removePipes, r'', _token) )  #@tokens;  ## removes outer pipes from the string.

          return tokens;



      def parseLisp(self, string):
          tokens = self.tokenize(string)
          parsed_expr = self.parseExpr(tokens)
          return  parsed_expr



      def parseExpr(self, tokens): 
          if not tokens:
              return []

          if tokens[0]== '(': 
             tokens.pop(0)
             list_elements = []

             while tokens[0] != ')' :
                toAddResult = self.parseExpr(tokens) 
                if toAddResult:
                    list_elements.append(toAddResult)

             tokens.pop(0)
             return  list_elements 
         
          elif not tokens[0] :
             tokens.pop(0)
             return []
          else :
             return tokens.pop(0)


      def send_query(self, query):
          self.makeSocket()
          self.soc.send(query)

      def retrieve_results(self): 
           data = '';
           while True:
               _data = self.soc.recv(1024)
               if not _data:
                   break
               data += _data

           return self.parseLisp(data)

      def retrieve_results_string(self): 
           data = '';
           results = []
           while True:
               _data = self.soc.recv(1024)
               if not _data:
                   break
               _data = _data.strip() 
               if _data != None:
                 results.append(_data)

           return ''.join(results)


      def getOrganismList(self):
         negPatterns = [ re.compile(r'^ECOBASES'), re.compile(r'^[#@]') ]
         posPatterns = [ re.compile(r'BASE$') ]
         
         query= "(mapcar #'object-name (all-orgs :all))"

         self.send_query(query)
         data = self.retrieve_results()
         
         if not data :
            return []

         organisms = []

         for _datum in data:
             if not _datum:
                continue
             datum = _datum.strip()
             for patt in negPatterns:
                 result = patt.search(datum)
                 if result:
                     continue 

             result = posPatterns[0].search(datum)
             if not result:
                continue 

             organisms.append(datum) 
         return organisms


      def stopPathwayTools(self):
         query= "(exit)"
         self.send_query(query)


      def wrap_query(self, function ):
          lisp = "(with-organism (:org-id\'%s) (mapcar #\'object-name(%s)))"\
                   %(self._organism, function)
          return lisp

      def call_func(self, function):

          self.send_query(self.wrap_query(function))

          result = self.retrieve_results()
          return result


      def genes_of_pathway(self, pathway, T):
         function = "genes-of-pathway \'%s" %(pathway)
         result = self.call_func(function)
         return result

      def genes_of_reaction(self, reaction, T):
         function = "genes-of-reaction \'%s" %(reaction)
         result = self.call_func(function)
         return result

      def protectFrameName(self, frame):
         pipePatt = re.compile(r'^\|.*\|$') ## if already pipe protected, don't do anything.
         status  = pipePatt.search(frame)
         if status:
           return frame

         if len(frame.strip())==0:
           return "|" + frame + "|"
        
         return frame


      def get_slot_values(self, frame, slot_name):
        try:
          frame = self.protectFrameName(frame)
          function = "get-slot-values \'%s \'%s" %(frame, slot_name)
          result = self.call_func(function)
          return result
        except:
          print frame, slot_name
          _exit(0)


      def get_slot_value(self, frame, slot_name):
         try:
           frame = self.protectFrameName(frame)
           function =  "get-slot-value \'%s \'%s" %(frame, slot_name)
           result = self.call_func_that_returns_string(function)
           return result
         except:
           print frame, slot_name
           _exit(0)

      def call_func_that_returns_string(self, function):
         # use for functions that will return a string and not a list. 
         # this function doesn't call mapcar and doesn't parse the returned list.
         query = "(with-organism (:org-id\'%s) (object-name (%s)))" %(self._organism, function);
         self.send_query (query)
         result = self.retrieve_results_string()
         return result

      def startPathwayTools(self):
         process = Process(target=startPathwayTools)
         process.start()
         time.sleep(5)


def startPathwayTools():
    cmd = "~/pathway-tools/pathway-tools -api"
    status = getstatusoutput(cmd)

if __name__=="__main__":


    pythonCyc = PythonCyc()

    try:
      pythonCyc.stopPathwayTools()
    except:
      print "nothing to stop"

    pythonCyc.startPathwayTools()


    print 'connecting'
    my_base_pathways = pythonCyc.call_func("all-pathways :all T")


    pwy_count=0
    unique_rxns ={}
    for pathway in my_base_pathways:
       pwy_count +=1
       mygenes = pythonCyc.genes_of_pathway(pathway,'T')
       totalrxns = pythonCyc.get_slot_values(pathway, "REACTION-LIST")

       for rxn in totalrxns:
          unique_rxns[rxn] = 1
     
       pathway_common_name = pythonCyc.get_slot_value(pathway,"common-name") 
       if not  pathway_common_name:
          pathway_common_name = "?"

       num_reactions = len(totalrxns)
       num_predicted_orfs = len(mygenes)
     
       num_covered_rxns =0
       num_genes =0
       orf_strings = {}

       for reaction in totalrxns :

          rngenes = pythonCyc.genes_of_reaction(reaction,"T")
          rxngenes = []
          for rngene in rngenes:
              rxngenes.append(pythonCyc.get_slot_value(rngene,"common-name"))

          if rxngenes: #this reaction is covered
              num_covered_rxns+= 1

          rxn_name = pythonCyc.get_slot_value(reaction,"common-name")
          if not rxn_name:
             rxn_name = '???'

          if reaction in unique_rxns:
             unique_rxns[reaction]= {}
             unique_rxns[reaction]['name'] = rxn_name
             unique_rxns[reaction]['ORFs'] = rxngenes
                  
             if rxngenes:
                unique_rxns[reaction]['covered']= 1 
             else :
                unique_rxns[reaction]['covered']= 0 

             unique_rxns[reaction]['num_pwys']= 1

          else :  
              unique_rxns[reaction]['covered']= 1
              unique_rxns[reaction]['num_pwys']+=1

          num_genes += len(rxngenes);
          for rxngene in rxngenes:
             orf_strings[rxngene] =1 

       # done for the reactions in a  pathway
       outputstr =  "PWY:" + "\t" + pathway + "\t" + pathway_common_name\
                    + "\t" + str(num_reactions) + "\t" + str(num_covered_rxns) + "\t" + str(num_predicted_orfs)

       for orf_string in orf_strings.keys():
           outputstr += "\t" + orf_string
       outputstr += "\n";

    # printing the reactions  MOST important paert of the calculation
    rxnOutputStr="";
    for  reaction in unique_rxns.keys(): 
       rxnOutputStr = reaction + "\t" + unique_rxns[reaction]['name']

       for orf in unique_rxns[reaction]['ORFs']:
          rxnOutputStr += "\t" + orf

       rxnOutputStr +="\n";
       print rxnOutputStr

    pythonCyc.stopPathwayTools()


