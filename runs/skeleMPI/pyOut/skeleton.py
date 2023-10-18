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
        for j in range(dim + 1):
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
        elif dim == 3 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    if(float(row[0]) != 0.0 and float(row[1]) != 0.0):
                        skeleton[i][0].append(float(row[0]))
                        skeleton[i][1].append(float(row[1]))
                        skeleton[i][2].append(float(row[2]))
                        skeleton[i][3].append(float(row[3]))
    return skeleton

def PlotSkeleton(skeleton,dim,time,np,ax,cmap):
    #Plots interface data according to process
    for i in range(np):
        color = cmap(i / np )
        if len(skeleton[i][0]) > 0:
            ax.scatter(skeleton[i][0],skeleton[i][1],s=4,color=color)
