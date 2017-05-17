import scipy
import time
from time import clock
from scipy.stats.stats import pearsonr
import sys
import shared
import struct
import csv
import numpy as np
import matplotlib.pyplot as plt

def main():
	# check the given arguments
	if len(sys.argv) < 3:
		usage()
	else:
		directory = sys.argv[1]
		image_name = sys.argv[2]
                source = sys.argv[3]

        with open(source,'r') as f:
                reader = csv.DictReader(f)
                her1mlist=[]
                her1wlist=[]
                deltamlist=[]
                deltawlist=[]
                mespamlist=[]
                mespawlist=[]
                for line in reader:
                        
                        info1 = line['12'].split('/')
                        info2 = line['13'].split('/')
                        info3 = line['14'].split('/')
                        her1m,her1w = float(info1[0]),float(info1[1])
                        deltam,deltaw = float(info2[0]),float(info2[1])
                        mespam,mespaw = float(info3[0]),float(info3[1])
                        her1mlist.append(her1m)
                        her1wlist.append(her1w)
                        deltamlist.append(deltam)
                        deltawlist.append(deltaw)
                        mespamlist.append(mespam)
                        mespawlist.append(mespaw)
                                
        wx1,wy1= plot_bar(her1wlist)
	wx2,wy2= plot_bar(deltawlist)
	wx3,wy3= plot_bar(mespawlist)
	#wx4,wy4= plot_bar([239.07])
	#wx5,wy5= plot_bar([239.07])
	
	mx1,my1= plot_bar(her1mlist)
	mx2,my2= plot_bar(deltamlist)
	mx3,my3= plot_bar(mespamlist)
	#mx4,my4= plot_bar([76.9071])
	#mx5,my5= plot_bar([100.725])
        N = 3
        men_means = (wx1, wx2,wx3)
        men_std = (wy1,wy2,wy3)

        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, men_means, width, color='#737373',yerr=men_std,error_kw=dict(elinewidth=2,ecolor='black',capsize =5, captick= 2))

        women_means = (mx1, mx2,mx3)
        women_std = (my1, my2,my3)
        rects2 = ax.bar(ind + width, women_means, width, color='#a6a6a6',yerr=women_std,error_kw=dict(elinewidth=2,ecolor='black',capsize =5, captick= 2))

        # add some text for labels, title and axes ticks
        ax.set_ylabel('Amplitude')
        ax.set_title('Wildtype/Her1Over mRNA Amplitude Compare')
        ax.set_xticks(ind)
        ax.set_xticklabels(('Her1,630', 'DeltaC,630','Mespa,630'))

        #ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))
        plt.savefig(directory + "/" + image_name + ".png", format = "png")
	plt.close()
	print 'Done. Your plot is stored in ' + directory + "/" + image_name + ".png"



def usage():
	print "Usage: python plot-cells.py <directory to store plots> <plot name> "
	exit(0)       
		
def plot_bar(sync_list):
        return np.mean(sync_list,axis=0), np.std(sync_list,axis=0)/np.sqrt(5)

main()



