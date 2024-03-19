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

def LoadGeometry(dim,time,root):
    geometry = []
    path = root + "geometry-{:1.1f}.txt".format(time)
    if dim == 2 and os.path.exists(path):
        i = 0
        with open(path,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            for row in data:
                geometry.append([])
                geometry[i].append(float(row[0]))
                geometry[i].append(float(row[1]))
                #geometry[i].append(float(row[2]))
                #geometry[i].append(float(row[3]))
                #geometry[i].append(float(row[4]))
                #geometry[i].append(float(row[5]))
                #geometry[i].append(float(row[6]))
                #geometry[i].append(float(row[7]))
                i += 1
    elif dim == 3 and os.path.exists(path):
        i = 0
        with open(path,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            for row in data:
                geometry.append([])
                geometry[i].append(float(row[ 0]))
                geometry[i].append(float(row[ 1]))
                geometry[i].append(float(row[ 2]))
                geometry[i].append(float(row[ 3]))
                geometry[i].append(float(row[ 4]))
                geometry[i].append(float(row[ 5]))
                geometry[i].append(float(row[ 6]))
                geometry[i].append(float(row[ 7]))
                geometry[i].append(float(row[ 8]))
                geometry[i].append(float(row[ 9]))
                geometry[i].append(float(row[10]))
                geometry[i].append(float(row[11]))
                i += 1
    return geometry

def PlotGeometry(geometry,dim,time,ax,cmap):
    #Plots interface data according to process
    #scaletop = 1.0
    #scalebot = 0.5
    #smax = 0
    #smin = 100000000
    #ref = 6
    #for i in range(len(geometry)):
    #    smax = max(smax,geometry[i][ref])
    #    smin = min(smin,geometry[i][ref])
    #norm = plt.Normalize(vmin=smin, vmax=smax)
    #sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    #sm.set_array([])  # This line is necessary for the color bar to work properly
    #cbar = plt.colorbar(sm, ax=ax)
    for i in range(len(geometry)):
        #color = cmap(norm(geometry[i][ref]))
        #ax.scatter(geometry[i][0],geometry[i][1],s=4,color=color)
        ax.scatter(geometry[i][0],geometry[i][1],s=1,color="black")


