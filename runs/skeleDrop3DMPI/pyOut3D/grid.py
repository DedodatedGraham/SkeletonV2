import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
import getopt
import os
import csv
import math

def LoadGrid(dim,time,np,root):
    #Loads in interface as [p0,p1,...,pnp]
    #each processor has its grid as  [[x,y,D]...]
    #we load in corners [[x-D,y-D,x+D,y+D]...]
    grid = []
    for i in range(np):
        grid.append([])
        for j in range(dim * 2):
            grid[i].append([])
        path = root + "vofinfo-{:.3f}-p{:d}.dat".format(time,i)
        if dim == 2 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    x = float(row[0])
                    y = float(row[1])
                    d = float(row[2]) / 2
                    grid[i][0].append(x-d)
                    grid[i][1].append(y-d)
                    grid[i][2].append(x+d)
                    grid[i][3].append(y+d)
        elif dim == 3 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    x = float(row[0])
                    y = float(row[1])
                    z = float(row[2])
                    d = float(row[3]) / 2
                    grid[i][0].append(x-d)
                    grid[i][1].append(y-d)
                    grid[i][2].append(z-d)
                    grid[i][3].append(x+d)
                    grid[i][4].append(y+d)
                    grid[i][5].append(z+d)
    return grid

def PlotGrid(grid,dim,time,np,ax,cmap):
    #Plots interface data according to process
    for i in range(np):
        if len(grid[i][0]) > 0:
            color = cmap(i / np)
            for j in range(len(grid[i][0])):
                print(j,'/',len(grid[i][0]))
                #loops through each index of our grid to plot each box
                plotx = [grid[i][0][j],grid[i][0][j],grid[i][3][j],grid[i][3][j],grid[i][0][j]]
                ploty = [grid[i][1][j],grid[i][4][j],grid[i][4][j],grid[i][1][j],grid[i][1][j]]
                ax.plot(plotx,ploty,grid[i][2][j],linewidth=0.5,color=color)#plot both squares on each z
                ax.plot(plotx,ploty,grid[i][5][j],linewidth=0.5,color=color)
                #finally plot connection of squares
                ax.plot([grid[i][0][j],grid[i][0][j]],[grid[i][1][j],grid[i][1][j]],[grid[i][2][j],grid[i][5][j]],linewidth=0.5,color=color)
                ax.plot([grid[i][0][j],grid[i][0][j]],[grid[i][4][j],grid[i][4][j]],[grid[i][2][j],grid[i][5][j]],linewidth=0.5,color=color)
                ax.plot([grid[i][3][j],grid[i][3][j]],[grid[i][1][j],grid[i][1][j]],[grid[i][2][j],grid[i][5][j]],linewidth=0.5,color=color)
                ax.plot([grid[i][3][j],grid[i][3][j]],[grid[i][4][j],grid[i][4][j]],[grid[i][2][j],grid[i][5][j]],linewidth=0.5,color=color)












