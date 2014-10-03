#!/usr/bin/python -w

import socket, re, sys, traceback
from os import _exit, path, remove
import re
from libs.python_modules.utils.sysutil import getstatusoutput
from multiprocessing import Process
from libs.python_modules.utils.metapathways_utils  import fprintf, printf, eprintf
import time

class PythonCyc:

    _organism = 'meta'
    soc = None
    _ptoolsExec = None

    def __init__(self):
        pass

    
    def setPToolsExec(self, ptoolsExec):
        self._ptoolsExec = ptoolsExec

    def setOrganism(self, organism):
        self._organism = organism

    def makeSocket(self):
        self.soc = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.soc.connect("/tmp/ptools-socket" )
        return True

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
        if self.makeSocket():
          self.soc.send(query)
        else:
           printf("ERROR\tCannot create or connect socket\n")

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

        _organisms = []
        for organism in organisms:
            _organisms.append(re.sub(r'BASE$', r'', organism).lower())

        return _organisms


    def sendStopSignal(self):
        query= "(exit)"
        self.send_query(query)
        time.sleep(10)


    def wrap_query(self, function ):
        lisp = "(with-organism (:org-id\'%s) (mapcar #\'object-name(%s)))" \
               %(self._organism, function)
        return lisp

    def call_func(self, function):

        self.send_query(self.wrap_query(function))
        result = self.retrieve_results()

        if result=='NIL':
            return []
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

    TIME = 10

    def sendStartSignal(self):
        #print "Starting up pathway tools"
        process = Process(target=startPathwayTools, args=(self._ptoolsExec,))
        process.start()

    def removeSocketFile(self):
       if remove("/tmp/ptools-socket"):
          return True
       else:
          return False

    def doesSocketExist(self):
       if path.exists("/tmp/ptools-socket"):
          #print "Socket exists"
          return True
       else:
          #print "No socket exist"
          return False

    def startPathwayTools(self):
        try:
           self.stopPathwayTools()
           trial = 0
           while not self.doesSocketExist() :
             self.sendStartSignal()
             time.sleep(self.TIME)
             if trial > 5:
               print "Failed to Start pathway-tools"
               return False
             else:
                trial += 1
        except:
            print traceback.print_exc(10) 
            print "Failed to Start pathway-tools"
            return False
        return True

    def stopPathwayTools(self):
        #print "Stopping the pathway tools"
        try:
           trial = 0
           while  self.doesSocketExist() :
             #print "Stopping the running pathway tools"
             self.sendStopSignal()
             time.sleep(self.TIME)
             if trial > 5:
               print "Failed to stop pathway-tools"
               if self.removeSocketFile():
                  return True
               else:
                  return False
             else:
                trial += 1
        except:
            print traceback.print_exc(10) 
            return False

#                self.getOrganismList()

    def getExpectedTaxonomicRange(self, pwy):
        """
        Given a pathway asks for its expected taxonomic range.
        """
        taxonomic_range = self.get_slot_values(pwy, "taxonomic-range")
        taxonomic_range_names = []
        for taxa in taxonomic_range:
            taxa_common_name = self.get_slot_value(taxa, "common-name")
            taxa = taxa.replace("TAX-","")
            taxonomic_range_names.append([taxa, taxa_common_name])
        return taxonomic_range_names


    def getAllPathways(self):
        """
        Function to get all pathways from the current organism.
        """
        pathway_list = []
        my_base_pathways = self.call_func("all-pathways :all T")
        for pathway in my_base_pathways:
            pathway_list.append(pathway)
        return pathway_list


    def getPathwayORFs(self, pwy):
        """
        Return a list of ORFs for a given pathway.
        """
        pathwayORFs = {}
        mygenes = self.genes_of_pathway(pwy,'T')
        for gene in mygenes:
            orf = self.get_slot_value(gene,"common-name")
            if orf not in pathwayORFs:
                orf = orf.strip("\"")
                pathwayORFs[orf] = 0

        return pathwayORFs.keys()

    def getPathwayReactionInfo(self, pwy):
        """
        Return a list of reactions given an input pathway shortname.
        :param pwy: string pathway name shortname
        :return: list of pathway reactions
        """
        pathwayRxns = {}
        myRxns = self.get_slot_values(pwy, "REACTION-LIST")
        coveredRxns = 0
        for rxn in myRxns:
            pathwayRxns[rxn] = 1
            rxn_genes = self.genes_of_reaction(rxn, "T")
            rxngenes_list = []
            for rxn_gene in rxn_genes:
                rxngenes_list.append(self.get_slot_value(rxn_gene,"common-name"))

            if rxngenes_list: # this reaction is covered
                coveredRxns += 1

        return [len(pathwayRxns.keys()), coveredRxns]



    def getReactionListLines(self):
        my_base_pathways = self.call_func("all-pathways :all T")
        pwy_count=0
        unique_rxns ={}
        #print "Extracting the reaction list"
        for pathway in my_base_pathways:
            # printf(" " + pathway)
            sys.stdout.flush()
            pwy_count +=1
            mygenes = self.genes_of_pathway(pathway,'T')
            totalrxns = self.get_slot_values(pathway, "REACTION-LIST")

            for rxn in totalrxns:
                unique_rxns[rxn] = 1

            pathway_common_name = self.get_slot_value(pathway,"common-name")
            if not  pathway_common_name:
                pathway_common_name = "?"

            num_reactions = len(totalrxns)
            num_predicted_orfs = len(mygenes)

            num_covered_rxns =0
            num_genes =0
            orf_strings = {}

            for reaction in totalrxns :

                rngenes = self.genes_of_reaction(reaction,"T")

                rxngenes = []
                for rngene in rngenes:
                    rxngenes.append(self.get_slot_value(rngene,"common-name"))

                if rxngenes: # this reaction is covered
                    num_covered_rxns+= 1

                rxn_name = self.get_slot_value(reaction,"common-name")
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
            outputstr =  "PWY:" + "\t" + pathway + "\t" + pathway_common_name \
                         + "\t" + str(num_reactions) + "\t" + str(num_covered_rxns) + "\t" + str(num_predicted_orfs)

            for orf_string in orf_strings.keys():
                outputstr += "\t" + orf_string
            outputstr += "\n";

        # printing the reactions MOST important part of the calculation
        rxnOutputStrings=[]

        for reaction in unique_rxns.keys():
            rxnOutputStr = reaction + "\t" + unique_rxns[reaction]['name']

            for orf in unique_rxns[reaction]['ORFs']:
                rxnOutputStr += "\t" + orf

            rxnOutputStrings.append(rxnOutputStr)
        return rxnOutputStrings


def startPathwayTools(ptoolsExec):
    cmd = ptoolsExec + " -api"
    status = getstatusoutput(cmd)


if __name__=="__main__":
    pythonCyc = PythonCyc()
    try:
        pythonCyc.stopPathwayTools()
    except:
        print "nothing to stop"

