## All this functions were grenerated to dynamically print output.
import sys
import numpy as np
import time

class Printer():
    """
    Print things to stdout on one line dynamically
    """ 
    def __init__(self):
        self.tic=time.clock()
        sys.stdout.flush()
        
    def printtextoneline(self,string):
        sys.stdout.write("\r\x1b[K"+string.__str__())
        sys.stdout.flush()
        
    def timepercentprint(self,minv,maxv,step,i):
        sys.stdout.write("\r 0% [{0}>]{1}% Time Elapsed: {2} s  ".format("="*int(i),round((float(i+1)/maxv)*100.0),round((time.clock()-self.tic))))
        sys.stdout.flush()
            
