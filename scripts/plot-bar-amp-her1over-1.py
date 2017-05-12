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

        wx1,wy1= plot_bar([88.6104])
	wx2,wy2= plot_bar([86.8536])
	wx3,wy3= plot_bar([295.609])
	#wx4,wy4= plot_bar([239.07])
	#wx5,wy5= plot_bar([239.07])
	
	mx1,my1= plot_bar([24.5567])
	mx2,my2= plot_bar([0.92375])
	mx3,my3= plot_bar([51.4278])
	#mx4,my4= plot_bar([76.9071])
	#mx5,my5= plot_bar([100.725])
        N = 3
        men_means = (wx1, wx2,wx3)
        men_std = (wy1,wy2,wy3)

        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, men_means, width, color='#737373',)

        women_means = (mx1, mx2,mx3)
        women_std = (my1, my2,my3)
        rects2 = ax.bar(ind + width, women_means, width, color='#a6a6a6')

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



'''
N = 5
men_means = (20, 35, 30, 35, 27)
men_std = (2, 3, 4, 1, 2)

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, men_means, width, color='r', yerr=men_std)

women_means = (25, 32, 34, 20, 25)
women_std = (3, 5, 2, 3, 3)
rects2 = ax.bar(ind + width, women_means, width, color='y', yerr=women_std)

# add some text for labels, title and axes ticks
ax.set_ylabel('Scores')
ax.set_title('Scores by group and gender')
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(('G1', 'G2', 'G3', 'G4', 'G5'))

ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))


def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)

plt.show()
'''
