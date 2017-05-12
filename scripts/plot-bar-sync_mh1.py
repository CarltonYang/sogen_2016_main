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

		
	x1,y1= plot_bar([0.906244, 0.910485, 0.884772, 0.928256, 0.898356 ])
	x2,y2= plot_bar([-0.142582, 0.0517386, -0.0447916, -0.0607525, 0.0444244])
	x3,y3= plot_bar([0.278456, 0.11652, 0.331941, 0.0593335, 0.312663])

	
        N = 3
        men_means = (x1, x2, x3)
        men_std = (y1,y2,y3)

        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, men_means, width, color='#808080', yerr=men_std,error_kw=dict(elinewidth=2,ecolor='black',capsize =5, captick= 2))

        #women_means = (25, 32, 34, 20, 25)
        #women_std = (3, 5, 2, 3, 3)
        #rects2 = ax.bar(ind + width, women_means, width, color='y', yerr=women_std)

        # add some text for labels, title and axes ticks
        ax.set_ylabel('Syncronization Scores')
        ax.set_title('Her1 mRNA Syncronization Scores')
        ax.set_xticks(ind)
        ax.set_xticklabels(('WT', 'Delta', 'DAPT'))

        #ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))

        '''
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
        #autolabel(rects2)
        '''
        plt.savefig(directory + "/" + image_name + ".png", format = "png")
	plt.close()
	print 'Done. Your plot is stored in ' + directory + "/" + image_name + ".png"



def usage():
	print "Usage: python plot-cells.py <directory to store plots> <plot name> "
	exit(0)       
		
def plot_bar(sync_list):
        return np.mean(sync_list,axis=0), np.std(sync_list,axis=0)/np.sqrt(5)

main()
