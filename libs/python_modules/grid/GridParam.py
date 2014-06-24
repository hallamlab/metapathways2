import sys


class GridParam:
     userid = None
     serviceAddress = None
     sample_name = None
     walltime= None
     os= None
     isaws=False
     awsparams={}
     batchSize=None
     queueSize=None
     working_directory='~'
     delay=None
     def __init__(self, userid, serviceAddress, walltime="10:00:00", os="linux", isaws=False, awsparams={}, \
           batchSize=500, queueSize=400, working_directory='~', delay=10000):
        self.userid = userid 
        self.serviceAddress = serviceAddress 
        self.walltime = walltime
        self.os = os
        self.isaws = isaws
        self.awsparams = awsparams
        self.batchSize = batchSize
        self.queueSize = queueSize
        self.working_directory = working_directory
        self.delay = delay

