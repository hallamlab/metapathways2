from test import hi as localtest
from libs.starcluster.test import hi as sctest
import sys

print sys.path

print "NOW RUNNING FROM SAME DIR"
localtest()

print "NOW RUNNING FROM SC"
sctest()


