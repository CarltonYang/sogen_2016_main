import scipy
import time
from time import clock
from scipy.stats.stats import pearsonr
import sys
import shared
import struct

import numpy as np
import matplotlib.pyplot as plt

def main():
	# check the given arguments
	if len(sys.argv) < 2:
		usage()
	else:
		directory = sys.argv[1]
		image_name = sys.argv[2]

        wx1,wy1= plot_bar([303.423, 259.566, 220.177, 244.611, 285.232])
	wx2,wy2= plot_bar([264.044, 225.314, 229.264, 277.473, 313.41])
	#wx3,wy3= plot_bar([264.044, 225.314, 229.264, 277.473, 313.41])
		
	mx1,my1= plot_bar([0, 0, 0, 0, 0])
	mx2,my2= plot_bar([0.00619674, 0.00238415, 0.00108412, 0.000448733, 0.000117573])
	#mx3,my3= plot_bar([0.00619674, 0.00238415, 0.00108412, 0.000448733, 0.000117573])
        N = 2
        men_means = (wx1, wx2)
        men_std = (wy1,wy2)

        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, men_means, width, color='#737373', yerr=men_std,error_kw=dict(elinewidth=2,ecolor='black',capsize =5, captick= 2))

        women_means = (mx1, mx2)
        women_std = (my1, my2)
        rects2 = ax.bar(ind + width, women_means, width, color='#a6a6a6', yerr=women_std, error_kw=dict(elinewidth=2,ecolor='black',capsize =5, captick= 2))

        # add some text for labels, title and axes ticks
        ax.set_ylabel('Amplitude')
        ax.set_title('Mespa mRNA Amplitude Compare')
        ax.set_xticks(ind)
        ax.set_xticklabels(('WT/Delta,600', 'WT/DAPT,720'))

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

