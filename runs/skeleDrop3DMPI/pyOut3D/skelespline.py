import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as nump
import sys
import getopt
import os
import csv
import math

def calcbfunc(n,points,t):
    if n == 1:
        return points[0]
    else:
        i = 0
        npoints = []
        while i < n - 1:
            npoints.append([])
            npoints[i].append((1-t)*points[i][0] + t * points[i + 1][0])
            npoints[i].append((1-t)*points[i][1] + t * points[i + 1][1])
            npoints[i].append((1-t)*points[i][2] + t * points[i + 1][2])
            i += 1
        return calcbfunc(n-1,npoints,t)

def LoadSkeletonSpline(dim,time,np,root):
    #Loads in splines of each important timestep
    spline = []
    for i in range(np):
        spline.append([])
        for j in range(dim + 1):
            spline[i].append([])
        path = root + "splineBranchDat-{:.3f}-P{:d}.dat".format(time,i)
        if dim == 2 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    deg = int(row[0])
                    for q in range(deg):
                        spline[i][0].append(float(row[q * (dim + 1) + 1]))
                        spline[i][1].append(float(row[q * (dim + 1) + 2]))
                        spline[i][2].append(float(row[q * (dim + 1) + 3]))
        elif dim == 3 and os.path.exists(path):
            with open(path,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                for row in data:
                    deg = int(row[0])
                    for q in range(deg):
                        spline[i][0].append(float(row[q * (dim + 1) + 1]))
                        spline[i][1].append(float(row[q * (dim + 1) + 2]))
                        spline[i][2].append(float(row[q * (dim + 1) + 3]))
                        spline[i][3].append(float(row[q * (dim + 1) + 4]))
    return spline

def PlotSkeletonSpline(spline,dim,time,np,ax,cmap):
    #Plots interface data according to process
    for i in range(np):
        color = cmap(i / np)
        if len(spline[i][0]) > 0:
            #gather control points
            control_points = []
            for tt in range(len(spline[i][0])):
                control_points.append([spline[i][0][tt],spline[i][1][tt],spline[i][2][tt]])
            ts = nump.linspace(0,1, num=1000)
            n = len(control_points)
            sx = []
            sy = []
            sr = []
            #gather spline points
            for param in ts:
                spline_point = calcbfunc(n,control_points,param)
                sx.append(spline_point[0])
                sy.append(spline_point[1])
                sr.append(spline_point[2])
            #finally plot each line
            for tt in range(len(sx) - 1):
                px = [sx[tt],sx[tt+1]]
                py = [sy[tt],sy[tt+1]]
                pr = (sr[tt] + sr[tt+1]) / 2
                ax.plot(px,py,linewidth=2,color=color)













