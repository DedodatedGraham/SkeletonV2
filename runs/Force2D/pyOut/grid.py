import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np
import sys
import getopt
import os
import csv
import math

def LoadGrid(dim,time,root):
    grid = []
    path = root + "scalarplot-{:1.1f}.txt".format(time)
    if dim == 2 and os.path.exists(path):
        i = 0
        with open(path,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            for row in data:
                grid.append([])
                x = float(row[0])
                y = float(row[1])
                d = float(row[2]) / 2
                grid[i].append(x-d)
                grid[i].append(y-d)
                grid[i].append(x+d)
                grid[i].append(y+d)
                grid[i].append(float(row[3]))
                grid[i].append(float(row[4]))
                grid[i].append(float(row[5]))
                grid[i].append(float(row[6]))
                i += 1
    elif dim == 3 and os.path.exists(path):
        i = 0
        with open(path,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            for row in data:
                grid.append([])
                x = float(row[0])
                y = float(row[1])
                z = float(row[2])
                d = float(row[3]) / 2
                grid[i].append(x-d)
                grid[i].append(y-d)
                grid[i].append(z-d)
                grid[i].append(x+d)
                grid[i].append(y+d)
                grid[i].append(z+d)
                grid[i].append(float(row[4]))
                grid[i].append(float(row[5]))
                grid[i].append(float(row[6]))
                grid[i].append(float(row[7]))
                i += 1
    return grid

def PlotGrid(grid,dim,time,ax,cmap):
    #Plots interface data according to process
    gcol = 7
    maxc = -10000
    minc = 10000
    for i in range(len(grid)):
        maxc = max(maxc,grid[i][gcol]);
        minc = min(minc,grid[i][gcol]);
    norm = plt.Normalize(vmin=minc, vmax=maxc)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # This line is necessary for the color bar to work properly
    cbar = plt.colorbar(sm, ax=ax)
    for i in range(len(grid)):
        color = cmap(norm(grid[i][gcol]))
        #loops through each index of our grid to plot each box
        plotx = [grid[i][0],grid[i][0],grid[i][2],grid[i][2],grid[i][0]]
        ploty = [grid[i][1],grid[i][3],grid[i][3],grid[i][1],grid[i][1]]
        #empty box
        #ax.plot(plotx,ploty,linewidth=0.5,color=color)
        #full box
        ax.fill(plotx, ploty, color=color)
