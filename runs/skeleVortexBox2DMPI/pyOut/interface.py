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

def LoadInterface(dim,time,np,root):
    #Loads in interface as [p0,p1,...,pnp]
    #each p is [[x],[y],[nx],[ny]]
    #and each list contains p0..pn for said pid :)
    interface = []
    for i in range(np):
        interface.append([])
        for j in range(dim*2 + 1):
            interface[i].append([])
        path = root + "interface-{:.3f}-p{:d}.dat".format(time,i)
        if dim == 2 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    interface[i][0].append(float(row[0]))
                    interface[i][1].append(float(row[1]))
                    interface[i][2].append(float(row[2]))
                    interface[i][3].append(float(row[3]))
                    interface[i][4].append(float(row[4]))
        elif dim == 3 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    interface[i][0].append(float(row[0]))
                    interface[i][1].append(float(row[1]))
                    interface[i][2].append(float(row[2]))
                    interface[i][3].append(float(row[3]))
                    interface[i][4].append(float(row[4]))
                    interface[i][5].append(float(row[5]))
    return interface

def PlotInterface(interface,dim,time,np,ax,cmap,maxk):
    #Plots interface data according to process
    for i in range(np):
        if len(interface[i][0]) > 0:
            #color = cmap(i / np)
            ax.scatter(interface[i][0],interface[i][1],s=1,color='black')
            #ax.quiver(interface[i][0],interface[i][1],interface[i][2],interface[i][3],scale=25.0,color='black')
