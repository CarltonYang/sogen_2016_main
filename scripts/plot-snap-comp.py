"""
Plots cells over time
Copyright (C) 2012 Ahmet Ay, Jack Holland, Adriana Sperlea

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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
	if len(sys.argv) < 5:
		usage()
	else:
		f = shared.openFile(sys.argv[1], "r")
		f2 = shared.openFile(sys.argv[2], "r")
		directory = sys.argv[3]
		image_name = sys.argv[4]
		step_size = shared.toFlo(sys.argv[5])
		start_time = sys.argv[6]
		
	
	print 'Plotting all the cells from ' + sys.argv[1] + '...'
	# split the lines to get data
	data = [line.split() for line in f]
	max_time = len(data) - 1
	data2 = [line.split() for line in f2]
	max_time2 = len(data2) - 1
	
	# calculate the tissue size
	cells_width = shared.toInt(data[0][0])
	cells_height = shared.toInt(data[0][1])
	total_cells = cells_width * cells_height + 1
	cells_width2 = shared.toInt(data2[0][0])
	cells_height2 = shared.toInt(data2[0][1])
	total_cells2 = cells_width2 * cells_height2 + 1
	
	# create a matrix to store the concentration values we obtain from the file
	cons = numpy.zeros(shape = (max_time, total_cells))
	cons2 = numpy.zeros(shape = (max_time2, total_cells2))
	
	# put the concentration values from the file into the matrix
	
	for i in range(1, max_time + 1):
		cons[i - 1][0] = shared.toFlo(data[i][0]) * step_size
		for j in range(1, total_cells):
			cons[i - 1][j] = shared.toFlo(data[i][j])

	for i in range(1, max_time2 + 1):
		cons2[i - 1][0] = shared.toFlo(data2[i][0]) * step_size
		for j in range(1, total_cells2):
			cons2[i - 1][j] = shared.toFlo(data2[i][j])
	# close the file
	f.close()
	f2.close()
	
	# plot colors
	colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
	color = 0
        start_time = float(start_time)/step_size
        pl.ylim((-1,250))
        f, axarr = pl.subplots(5,sharex=True , sharey=True)
        pl.tight_layout()
        col = list(range(50))
        row5=[0]*50
        
        for i in range(5):
                row1 = cons[start_time-1,1:51]
                row2 = cons[start_time-1,51:101]
                row3 = cons[start_time-1,101:151]
                row4 = cons[start_time-1,151:201]
                print row1
                print row2
                print row3
                print row4
                print type(row1)
                print len(row1), len(row2),len(row3),len(row4)
                print row1[0]
                for j in range(len(row1)):
                        row5[j]=(row1[j]+row2[j]+row3[j]+row4[j])/4
                axarr[i].plot(col,row4,'r')

                row1 = cons2[start_time-1,1:51]
                row2 = cons2[start_time-1,51:101]
                row3 = cons2[start_time-1,101:151]
                row4 = cons2[start_time-1,151:201]
                print row1
                print row2
                print row3
                print row4
                print type(row1)
                print len(row1), len(row2),len(row3),len(row4)
                print row1[0]
                for j in range(len(row1)):
                        row5[j]=(row1[j]+row2[j]+row3[j]+row4[j])/4
                axarr[i].plot(col,row4,'b')
                #axarr[i].plot(col,row2,'g')
                #axarr[i].plot(col,row3,'b')
                #axarr[i].plot(col,row4,'c')
                title = "Snapshot at " + str(int(start_time*step_size))
                axarr[i].set_title(title )
                start_time = start_time+ 600
                
	#pl.ylim((-1,200))
	pl.savefig(directory + "/" + image_name + ".png", format = "png")
	pl.close()
	print 'Done. Your plot is stored in ' + directory + "/" + image_name + ".png"

def usage():
	print "Usage: python plot-cells.py <file with concentration levels> <2nd file with concentration levels> <directory to store plots> <plot name> <step size> <starting time(in minutes)>"
	exit(0)

main()

