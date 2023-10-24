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
#import skeleton.py
#import skelespline.p
from interface import LoadInterface,PlotInterface
from skeleton import LoadSkeleton,PlotSkeleton
from skelespline import LoadSkeletonSpline,PlotSkeletonSpline
import multiprocessing as mp


#here we run our functions we want
def Load(dim,time,np,root):
    data = []
    datatag = []#data tags; allows for turning on/off different sections
    #0 => interface
    #1 => skeleton(points)
    #2 => skeleton(spline)

    #Load interface
    interface = LoadInterface(dim,time,np,root)
    data.append(interface)
    datatag.append(0)

    #Load Skeleton
    skeleton = LoadSkeleton(dim,time,np,root)
    data.append(skeleton)
    datatag.append(1)

    #Load SkeletonSpline
    skeletonspline = LoadSkeletonSpline(dim,time,np,root)
    data.append(skeletonspline)
    datatag.append(2)

    #finally return data we loaded
    return data,datatag


def Plot(data,datatag,dim,time,np,ax):
    #cmap = matplotlib.colormaps.get_cmap("tab20") #for plotting different pid
    cmap = matplotlib.colormaps.get_cmap("brg") #for plotting smooth
    for i in range(len(data)):
        did = datatag[i]
        if did == 0:
            maxx = 0
            minx = 1000000
            maxy = 0
            miny = 1000000
            maxk = 0
            cx = 0
            cy = 0
            scale = 0
            cnt = 0
            for j in range(len(data[i])):
                #loop through np
                for q in range(len(data[i][j][0])):
                    if(not math.isnan(data[i][j][0][q])):
                        cx += data[i][j][0][q]
                    if(not math.isnan(data[i][j][1][q])):
                        cy += data[i][j][1][q]
                    if(data[i][j][0][q] > maxx):
                        maxx = data[i][j][0][q]
                    if(data[i][j][0][q] < minx):
                        minx = data[i][j][0][q]
                    if(data[i][j][1][q] > maxy):
                        maxy = data[i][j][1][q]
                    if(data[i][j][1][q] < miny):
                        miny = data[i][j][1][q]
                    if(data[i][j][4][q] > maxk):
                        maxk = data[i][j][4][q]
                    cnt += 1
            cx /= cnt
            cy /= cnt
            scale = max(maxx-minx,maxy-miny) / 2
            #we define our center based on our interface aswell
            ax.set_xlim(cx - scale,cx + scale)
            ax.set_ylim(cy        ,cy + scale + 0.15)
            PlotInterface(data[i],dim,time,np,ax,cmap,maxk)
        elif did == 1:
            PlotSkeleton(data[i],dim,time,np,ax,cmap)
        elif did == 2:
            PlotSkeletonSpline(data[i],dim,time,np,ax,cmap)
        elif did == 3:
            print("not implemented yet")
def Save(root,plt,time):
    save = root + r'plot-{:3.5f}.png'.format(time)
    plt.savefig(save,dpi=1000)
def main():
    fig = plt.figure()
    #LoadInformation
    np = int(sys.argv[1])
    dim = int(sys.argv[2])
    start_time = float(sys.argv[3]) / 1000
    max_time = float(sys.argv[4]) / 1000
    time = start_time
    dt = 0.001
    root = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + r'/'
    datalocation = root + r'dat/'
    savelocation = root + r'pyOut/Img/'

    case = True
    mode = 1
    while case:
        if mode == 0:
            try:
                if(time > max_time):
                    case = False
                    break
                ax = fig.add_axes([0.1,0.1,0.85,0.85])
                data,datatag = Load(dim,time,np,datalocation)
                Plot(data,datatag,dim,time,np,ax)
                Save(savelocation,plt,time)
                plt.clf()
                print("finished {:f}".format(time))
                time += dt
            except:
                case = False
                print("no more data")
        else:
            #run for debug
            ax = fig.add_axes([0.1,0.1,0.85,0.85])
            data,datatag = Load(dim,time,np,datalocation)
            Plot(data,datatag,dim,time,np,ax)
            Save(savelocation,plt,time)
            plt.clf()
            print("finished {:f}".format(time))
            time += dt

#run main
if __name__ == "__main__":
    main()
