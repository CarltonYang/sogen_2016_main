
import numpy
import scipy
import pylab as pl
import time
from time import clock
from scipy.stats.stats import pearsonr
import sys
import shared
def main():
	# check the given arguments
	if len(sys.argv) < 5:
		usage()
	else:
		f = shared.openFile(sys.argv[1], "r")
		directory = sys.argv[2]
		image_name = sys.argv[3]
		threshold = shared.toFlo(sys.argv[4])
                
	print 'Extrating all parameter sets with score less than ' + str(threshold) +' from ' + sys.argv[1]
        fw = shared.openFile(directory+image_name+'.txt', 'w')

        first_line = True
	para_count = 0
	current_set = ''
	
	data = [line.split() for line in f]
	for line in range(len(data)):
                
                if 'generation:' in data[line] and 'best' in data[line] and 'fitness:' in data[line]:
                        fitness = float(data[line][-1])
                        if fitness <= threshold:
                                if not current_set== data[line+1][2]:
                                        current_set = data[line+1][2]
                        
                                        if not first_line:
                                                fw.write('\n')
                                        fw.write(current_set)
                                        first_line = False
                                        para_count+=1
        print 'All parameter sets below threshold are saved into '+ directory + image_name+ '.txt'
        f.close()
        fw.close()

        
def usage():
	print "Usage: python extract_para_sets.py <file with parameter sets> <directory to store extracted parameters> <file name> <threshold for fitness score>"
	exit(0)

main()

