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




'''
HOW TO USE
python plot-cells.py <file with concentration levels> <directory to store plots> <plot name> <step size> <plot type> <position> <number of column(col)/start time(col_t)>
last one is only used in plotting multiple columns, enter 1 otherwise
plot types are cell, one cell over time , position will be cell#
               col, one or more column over time, position will be starting column, number of column is how many columns following
               row, all rows at some specific time step, position will be time step, 0-#of min *100
               all, entire PSM over time, position and # of column not important, set to 1
               col_t, keep track of concentration level of a column over time(as the PSM grows), position will be the starting position,time step as specified
                    (0-30000 is the same and for our use, only column>=10 can be chosen over 30000)

'''
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
	if len(sys.argv) < 6:
		usage()
	else:
		f = shared.openFile(sys.argv[1], "r")
		directory = sys.argv[2]
		image_name = sys.argv[3]
		step_size = shared.toFlo(sys.argv[4])
		plot_style= sys.argv[5]
		plot_helper1= int(shared.toFlo(sys.argv[6]))
		plot_helper2= int(shared.toFlo(sys.argv[7]))


        print 'Plotting all the cells from ' + sys.argv[1] + '...'
	# split the lines to get data
	data = [line.split() for line in f]
	max_time = len(data) - 1
	
	# calculate the tissue size
	cells_width = shared.toInt(data[0][0])
	cells_height = shared.toInt(data[0][1])
	total_cells = cells_width * cells_height + 1
        #print cells_width
	# create a matrix to store the concentration values we obtain from the file
	cons = numpy.zeros(shape = (max_time, total_cells))
        cons_t = [0]*max_time
        # create array for row plotting
	pos = [0]*50
        for i in range (0, 49):
            pos[i]=i
        time = [0]*max_time
        for i in range (1,max_time+1):
            time[i-1]=i*step_size
	# put the concentration values from the file into the matrix
	for i in range(1, max_time + 1):
		cons[i - 1][0] = shared.toFlo(data[i][0]) * step_size
		for j in range(1, total_cells):
			cons[i - 1][j] = shared.toFlo(data[i][j])
	
	# close the file
	f.close()
	
	# plot colors
	colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
	color = 0
        
        # decide which row/column/cells to plot
        #plot_col = 1
        if (plot_style == "col"):
            startpoint = plot_helper1
            interval = cells_width 
        elif (plot_style == "all"):
            startpoint =1
            interval =1
        elif (plot_style =="cell"):
            startpoint = plot_helper1
            interval = total_cells
        elif (plot_style == "col_t"):
            startpoint = plot_helper1
            interval= cells_width
            start_time = plot_helper2
            
        if (plot_style!= "row" and plot_style!="col_t"):    
            print "not in row"
            for j in range(0,plot_helper2):
	       for i in range(startpoint+j, total_cells, interval):
        
		  start = 0
		
		# Adjust the plotting interval for each cell to account for different columns being staggered
		# as they enter the PSM at intervals of 6 minutes apart frome ach other
	  	  while cons[start][i] == -1: # -1 stands for no data in the output file
		 	    start += 1
		  end = max_time - 1
		  while cons[end][i] == -1:
			 end -= 1;

		  if (i % 4 == 0):
			 pl.plot(cons[start:end, 0], cons[start:end, i], 'r')
		  elif (i % 4 == 1):
			 pl.plot(cons[start:end, 0], cons[start:end, i], 'g')
		  elif (i % 4 == 2):
			 pl.plot(cons[start:end, 0], cons[start:end, i], 'b')
		  else:
			 pl.plot(cons[start:end, 0], cons[start:end, i], 'c')
        
	elif (plot_style == "row"):
            print "in row"
            pl.plot(pos[0:49],cons[plot_helper1, 1:50],'r')
            pl.plot(pos[0:49],cons[plot_helper1, 51:100],'g')
	    pl.plot(pos[0:49],cons[plot_helper1, 101:150],'b')
            pl.plot(pos[0:49],cons[plot_helper1, 151:200],'c')


        elif (plot_style == "col_t"):
            difference=0
            if start_time > 29999:
                difference= start_time-29999
            print "tracking column"
            if (startpoint <=9):
                if start_time<29999:
                    timetodie= 24000
                    timetilappear=0   
                    i=1
                    while i<=29999:
                        #print i
                        #print startpoint
                        #print cons[timetilappear+i,startpoint]
                        cons_t[timetilappear+i]=cons[timetilappear+i,startpoint]
                        i+=1
                    timetilappear=29999
                else:
                    for i in range (start_time):
                        cons_t[i]=0
                    timetilappear=start_time
                    timetodie= (cells_width-1- startpoint)*600-1
            elif (startpoint>9 and startpoint <=49): 
                print ">9"
                timetodie= (cells_width-1- startpoint)*600
                timetilappear= (startpoint +1-10)*600+29999+difference
                #print timetodie
                #print timetilappear
                for i in range (timetilappear):
                    cons_t[i]=0      
            growth=1
            i=1
            while i<=timetodie:
                #print i
                #print growth
                #print startpoint
                #print cons[timetilappear+i,startpoint]
                if (timetilappear+i>= 90000):
                    break
                if (startpoint>= 93000):
                    break
                cons_t[timetilappear+i]=cons[timetilappear+i,startpoint]
                growth+=1
                i+=1
                if (growth %600==0 ):
                    startpoint+=1
                    growth=1
            pl.plot(time[0:max_time], cons_t, 'r')
            
        
        
       
        #new_array=[500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840]
        #pl.grid(True)
        #pl.axis([600, 700,0,1100])
        #pl.xticks(new_array)
	pl.savefig(directory + "/" + image_name + ".png", format = "png")
	pl.close()
	print 'Done. Your plot is stored in ' + directory + "/" + image_name + ".png"

def usage():
	print "Usage: python plot-cells.py <file with concentration levels> <directory to store plots> <plot name> <step size>"
	exit(0)

main()

