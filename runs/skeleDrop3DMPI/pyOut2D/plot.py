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
from skeleton import LoadSkeleton,PlotSkeleton,LoadThinSkeleton,PlotThinSkeleton
from skelespline import LoadSkeletonSpline,PlotSkeletonSpline
from grid import LoadGrid,PlotGrid
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
    #skeletonspline = LoadSkeletonSpline(dim,time,np,root)
    #data.append(skeletonspline)
    #datatag.append(2)

    #Load vofGrid
    #grid = LoadGrid(dim,time,np,root)
    #data.append(grid)
    #datatag.append(3)

    #Load ThinSkeleton
    #tskeleton = LoadThinSkeleton(dim,time,np,root)
    #data.append(tskeleton)
    #datatag.append(4)

    #finally return data we loaded
    return data,datatag


def Plot(data,datatag,dim,time,np,ax):
    #cmap = matplotlib.colormaps.get_cmap("tab20") #for plotting different pid
    smoothcmap = matplotlib.colormaps.get_cmap("winter") #for plotting smooth
    jaggcmap = matplotlib.colormaps.get_cmap("tab20") #for plotting pid
    plt.title("time:{:3.3f}".format(time))
    for i in range(len(data)):
        did = datatag[i]
        if did == 0:
            maxk = 0
            ax.set_xlim(2,4)
            ax.set_ylim(-1,1)
            PlotInterface(data[i],dim,time,np,ax,jaggcmap,maxk)
        elif did == 1:
            PlotSkeleton(data[i],dim,time,np,ax,smoothcmap)
            #PlotSkeleton(data[i],dim,time,np,ax,jaggcmap)
        elif did == 2:
            PlotSkeletonSpline(data[i],dim,time,np,ax,smoothcmap)
        elif did == 3:
            PlotGrid(data[i],dim,time,np,ax,jaggcmap)
        elif did == 4:
            PlotThinSkeleton(data[i],dim,time,np,ax,jaggcmap)
        else:
            print("not implemented yet")
def Save(root,plt,time):
    timeform = '{:3.3f}'.format(time)
    save = root + r'plot-' + timeform.zfill(6) + r'.png'.format(time)
    plt.savefig(save,dpi=1000)
def main(fig,np,dim,start_time,max_time,root,datalocation,savelocation):
    time = start_time
    print("starting at",start_time)
    print("ending at",max_time)
    case = True
    dt = 0.001
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
            if(time > max_time):
                case = False
                break
            #run for debug
            ax = fig.add_axes([0.1,0.1,0.85,0.85])
            print("trying {:f}".format(time))
            data,datatag = Load(dim,time,np,datalocation)
            Plot(data,datatag,dim,time,np,ax)
            Save(savelocation,plt,time)
            plt.clf()
            print("finished {:f}".format(time))
            time += dt

#run main
if __name__ == "__main__":
    #LoadInformation
    np = int(sys.argv[1])
    dim = int(sys.argv[2])
    start_time = float(sys.argv[3]) / 1000
    max_time = float(sys.argv[4]) / 1000
    mode = int(sys.argv[5])
    num_processes = int(mp.cpu_count()/2)
    print(num_processes)
    root = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + r'/'
    dt = 0.001
    datalocation = root + r'dat/'
    savelocation = root + r'pyOut2D/Img/'
    if(mode == 0):
        fig = plt.figure()
        main(fig,np,dim,start_time,max_time,root,datalocation,savelocation)
    else:
        nps = []
        dims = []
        starts = []
        maxs = []
        figs = []
        roots = []
        datas = []
        saves = []
        totalt = max_time - start_time
        avgdt = totalt / dt # this gives us the amount of time steps now
        avgdt = math.ceil(avgdt / num_processes)
        for i in range(num_processes):
            figs.append(plt.figure())
            nps.append(np)
            dims.append(dim)
            if i == 0:
                #first case
                starts.append(start_time)
                maxs.append(starts[i] + avgdt * dt)
            elif i == num_processes - 1:
                #last case
                starts.append(maxs[i - 1])
                maxs.append(max_time)
            else:
                #between case
                starts.append(maxs[i - 1])
                maxs.append(starts[i] + avgdt * dt)
            roots.append(root)
            datas.append(datalocation)
            saves.append(savelocation)
        with mp.Pool(processes=num_processes) as pool:
            pool.starmap(main, zip(figs,nps,dims,starts,maxs,roots,datas,saves))
