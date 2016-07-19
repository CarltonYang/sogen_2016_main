import numpy
import scipy
import pylab as pl
import time
from time import clock
from scipy.stats.stats import pearsonr
import sys
import shared
import struct


def main():
	# check the given arguments
	with open(sys.argv[1]) as f:
            content = f.readlines()
	#f = shared.openFile(, "r")
	#directory = sys.argv[2]
	'''	
	data = [line.split() for line in f]
        max_time = len(data) - 1

        # calculate the tissue size
	cons = numpy.zeros(shape = (max_time, 71))
        
        # create array for row plotting

	# put the concentration values from the file into the matrix
	for i in range(0, len(data)):
		for j in range(0, 71):
			cons[i][j] = shared.toFlo(data[i][j])
                        print cons[i][j]
	'''
	#print content
        f.close()
        wr= shared.openFile("newset.params","w")
        newline=""
        needlist= [0,2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,29,30,31,33,38,39,40,42,47,48,49,51,53,54,56,57,59,60,62,63,65,66,67,68,69,70]
        #print len(needlist)
        for i in range(0, len(content)):
            #print content[i]
            newline=""
            newlist= [x.strip() for x in content[i].split(',')]
            #print len(newlist)
            #print len(needlist)
            for j in range(0, len(newlist)):
                if j in needlist:
                    if newline=="":
                         newline+= newlist[j]
                    else:
                         newline+=","+newlist[j]
            wr.write(newline)
            wr.write("\n") 
        
        wr.close()
        
main()