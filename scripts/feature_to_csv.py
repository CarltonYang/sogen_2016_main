import numpy
import scipy
import pylab as pl
import time
from time import clock
from scipy.stats.stats import pearsonr
import sys
import shared
import struct
import csv


def main():
	# check the given arguments
	with open(sys.argv[1]) as f, open('features.csv', 'w') as fw:
                fw = csv.DictWriter(fw, ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21'])
                fw.writeheader()
                for line in f:
                        
                        line = line.split()
                        #print len(line)
                        fw.writerow({'1':line[0],
                                    '2':line[1],
                                    '3':line[4],
                                    '4':line[5],
                                    '5':line[6],
                                    '6':line[7],
                                    '7':str(line[2])+"/"+str(line[3]),
                                    '8':line[8],
                                    '9':str(line[9])+"/"+str(line[10]),
                                    '10':str(line[11])+"/"+str(line[12]),
                                    '11':line[13],
                                    '12':str(line[14])+"/"+str(line[15]),
                                    '13':str(line[16])+"/"+str(line[17]),
                                    '14':str(line[18])+"/"+str(line[19]),
                                    '15':line[20],
                                    '16':line[21],
                                    '17':str(line[22])+"/"+str(line[23]),
                                    '18':str(line[24])+"/"+str(line[25]),
                                    '19':line[26],
                                    '20':str(line[27])+"/"+str(line[28]),
                                    '21':str(line[29])+"/"+str(line[30])})
                        #0       1   2      3       4  5     6         7         8        9        10     11  12    13       14      15       16      17
                        #28.2708 1 30.3207 28.2708 0 0.8134 -0.540569 -0.689723 0.105733 104.724 124.619 0 127.607 0.403748 14.2342 104.657 1.22888 38.6712
                        #18       19      20       21       22     23      24         25      26      27      28     29       30  
                        #27.1467 91.2622 0.487785 0.334792 101.69 123.696 0.00204714 118.191 0.18401 27.5067 114.584 19.3371 114.584
main()
