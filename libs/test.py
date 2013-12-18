import sys

from  libs.starcluster.cli  import cli
import sys

def teststarcluster():
    print sys.path
    clustername = 'smallcluster'
    
    master = cli( ['listmaster', clustername] )
    if master==None:
       print  "No cluster previously setup"
       start_cluster = cli.cli( ['start', clustername, '-c',  '../amazon_aws.config'  ] )
       master_name = cli.cli( ['listmaster', clustername] )
       if master_name:
          print  "Successfully created a cluster with X nodes"
    else:
       print  "Found a cluster with master name %s" %(master) 



def hi():
    print "hi my path is " , sys.path

