__author__ = 'Rohin Kumar Y'
#from fileios import *
#msg = 'Enter Absolute Path to file: '
#f_name = raw_input(msg).strip()

#path = file_data_and_path(f_name)
#if path != None:
#       print 'Path:',path
#import tpcf
from fileios import *

t1 = checkdatfile('./testfile.dat')
print t1

t2=inputfiles('data')
print t2
