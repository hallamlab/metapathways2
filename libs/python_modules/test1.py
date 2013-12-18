#from libs.test import hi as localtest
from libs.starcluster.test import teststarcluster as sctest
import sys

print sys.path

print "NOW RUNNING FROM SAME DIR"
#localtest()

print "NOW RUNNING FROM SC"
sctest()


