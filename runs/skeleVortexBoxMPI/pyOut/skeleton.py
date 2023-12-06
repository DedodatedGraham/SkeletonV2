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

def LoadSkeleton(dim,time,np,root):
    #Loads in interface as [p0,p1,...,pnp]
    #each p is [[x],[y],[nx],[ny]]
    #and each list contains p0..pn for said pid :)
    skeleton = []
    for i in range(np):
        skeleton.append([])
        for j in range(dim + 3):
            skeleton[i].append([])
        path = root + "skeletonscatter-{:.3f}-p{:d}.dat".format(time,i)
        if dim == 2 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    if(float(row[0]) != 0.0 and float(row[1]) != 0.0):
                        skeleton[i][0].append(float(row[0]))
                        skeleton[i][1].append(float(row[1]))
                        skeleton[i][2].append(float(row[2]))
                        skeleton[i][3].append(float(row[3]))
                        skeleton[i][4].append(float(row[4]))
        elif dim == 3 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    if(float(row[0]) != 0.0 and float(row[1]) != 0.0):
                        skeleton[i][0].append(float(row[0]))
                        skeleton[i][1].append(float(row[1]))
                        skeleton[i][2].append(float(row[2]))
                        skeleton[i][3].append(float(row[3]))
                        skeleton[i][4].append(float(row[4]))
                        skeleton[i][5].append(float(row[5]))
    return skeleton
def LoadThinSkeleton(dim,time,np,root):
    #Loads in interface as [p0,p1,...,pnp]
    #each p is [[x],[y],[nx],[ny]]
    #and each list contains p0..pn for said pid :)
    skeleton = []
    for i in range(np):
        skeleton.append([])
        for j in range(dim + 3):
            skeleton[i].append([])
        path = root + "skeletonthin-{:.3f}-p{:d}.dat".format(time,i)
        if dim == 2 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    if(float(row[0]) != 0.0 and float(row[1]) != 0.0):
                        skeleton[i][0].append(float(row[0]))
                        skeleton[i][1].append(float(row[1]))
                        skeleton[i][2].append(float(row[2]))
                        skeleton[i][3].append(float(row[3]))
                        skeleton[i][4].append(float(row[4]))
        elif dim == 3 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    if(float(row[0]) != 0.0 and float(row[1]) != 0.0):
                        skeleton[i][0].append(float(row[0]))
                        skeleton[i][1].append(float(row[1]))
                        skeleton[i][2].append(float(row[2]))
                        skeleton[i][3].append(float(row[3]))
                        skeleton[i][4].append(float(row[4]))
                        skeleton[i][5].append(float(row[5]))
    return skeleton

def PlotSkeleton(skeleton,dim,time,np,ax,cmap):
    #Plots interface data according to process
    maxcolor = 0;
    mincolor = 1000000;
    #we track max and min for a colorbar
    for i in range(np):
        if len(skeleton[i][0]) > 0:
            colortrace = []
            #scatter whole skeleton at once
            #ax.scatter(skeleton[i][0],skeleton[i][1],s=4,color=color)
            for q in range(len(skeleton[i][0])):
                ##colortrace.append(skeleton[i][2][q] / skeleton[i][4][q])
                colortrace.append(skeleton[i][2][q])
                color = cmap(colortrace[q])
                #color = cmap(i / np)
                ax.scatter(skeleton[i][0][q],skeleton[i][1][q],s=4,color=color)
            #after the plotting we set max and min if possible
            tmax = max(colortrace)
            tmin = min(colortrace)
            if tmax > maxcolor:
                maxcolor = tmax
            if tmin < mincolor:
                mincolor = tmin
    #finally we plot our colorbar :)
    norm = plt.Normalize(vmin=mincolor, vmax=maxcolor)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # This line is necessary for the color bar to work properly
    cbar = plt.colorbar(sm, ax=ax)

def PlotThinSkeleton(skeleton,dim,time,np,ax,cmap):
    #Plots  data according to process
    for i in range(np):
        if len(skeleton[i][0]) > 0:
            for q in range(len(skeleton[i][0])):
                color = cmap(i / np)
                ax.scatter(skeleton[i][0][q],skeleton[i][1][q],s=4,color=color)












