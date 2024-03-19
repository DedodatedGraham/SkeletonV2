import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
from grid import LoadGrid,PlotGrid
from geometry import LoadGeometry,PlotGeometry
import multiprocessing as mp
from matplotlib.colors import LinearSegmentedColormap

#here we run our functions we want
def Load(dim,time,root):
    data = []
    datatag = []#data tags; allows for turning on/off different sections
    #0 => grid
    #1 => geometry

    #Load Grid
    grid = LoadGrid(dim,time,root)
    data.append(grid)
    datatag.append(0)

    #Load ThinSkeleton
    geometry = LoadGeometry(dim,time,root)
    data.append(geometry)
    datatag.append(1)

    #finally return data we loaded
    return data,datatag

def truncate_colormap(cmap,minv,maxv,n):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minv, b=maxv),cmap(np.linspace(minv, maxv, n)))
    return new_cmap
def Plot(data,datatag,dim,time,ax):
    #cmap = matplotlib.colormaps.get_cmap("tab20") #for plotting different pid
    smoothcmap1 = matplotlib.colormaps.get_cmap("inferno") #for plotting smooth
    smoothcmap2 = matplotlib.colormaps.get_cmap("viridis") #for plotting smooth
    #smoothcmap2 = truncate_colormap(smoothcmap2, 0.5, 0.75, 250)
    jaggcmap = matplotlib.colormaps.get_cmap("tab20") #for plotting pid
    #plt.title("time:{:3.3f}".format(time))
    ax.axis('off')
    for i in range(len(data)):
        did = datatag[i]
        if did == 0:
            PlotGrid(data[i],dim,time,ax,smoothcmap1)
        elif did == 1:
            PlotGeometry(data[i],dim,time,ax,smoothcmap2)
        else:
            print("not implemented yet")
def Save(root,plt,time):
    plt.gca().set_aspect('equal', adjustable='box')
    timeform = '{:1.1f}'.format(time)
    save = root + r'plot-' + timeform.zfill(2) + r'.png'.format(time)
    plt.savefig(save,dpi=1000)
def main(fig,dim,start_time,max_time,root,datalocation,savelocation):
    time = start_time
    print("starting at",start_time)
    print("ending at",max_time)
    case = True
    dt = 0.1
    mode = 1
    while case:
        if mode == 0:
            try:
                if(time > max_time):
                    case = False
                    break
                ax = fig.add_axes([0.1,0.1,0.85,0.85])
                data,datatag = Load(dim,time,datalocation)
                Plot(data,datatag,dim,time,ax)
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
            data,datatag = Load(dim,time,datalocation)
            Plot(data,datatag,dim,time,ax)
            Save(savelocation,plt,time)
            plt.clf()
            print("finished {:f}".format(time))
            time += dt

#run main
if __name__ == "__main__":
    #LoadInformation
    dim = int(sys.argv[2])
    start_time = float(sys.argv[3]) / 10
    max_time = float(sys.argv[4]) / 10
    mode = int(sys.argv[5])
    num_processes = int(mp.cpu_count()/2 - 11)
    print(num_processes)
    root = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + r'/'
    dt = 0.1
    datalocation = root + r'dat/'
    savelocation = root + r'pyOut/Img/'
    if(mode == 0):
        fig = plt.figure()
        main(fig,dim,start_time,max_time,root,datalocation,savelocation)
    else:
        dims = []
        starts = []
        maxs = []
        figs = []
        roots = []
        datas = []
        saves = []
        totalt = max_time - start_time
        dp = math.floor(totalt * 10 / num_processes) #avg dt is amount of time
        for i in range(num_processes):
            figs.append(plt.figure())
            dims.append(dim)
            starts.append(i*dp / 10)
            if i == num_processes - 1:
                maxs.append(max_time)
            else:
                maxs.append((i+1)*dp / 10)
            roots.append(root)
            datas.append(datalocation)
            saves.append(savelocation)
        with mp.Pool(processes=num_processes) as pool:
            pool.starmap(main, zip(figs,dims,starts,maxs,roots,datas,saves))
