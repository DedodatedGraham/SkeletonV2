import matplotlib.pyplot as plt
import sys
import getopt
import os
import csv



def plotmirror(sdat,idat,mirroraxi,mirrorval,save,t):
    #We will mirror results here
    plt.clf()
    plt.xlim(1.5,4)
    plt.ylim(-1.25,1.25)
    vmin = 0
    vmax = 0.6
    plt.rcParams['figure.dpi'] = 1000
    if mirroraxi == 0:
        print('not implemented')
    else:
        plt.scatter(idat[0],idat[1],color='black',s=5,vmin=vmin,vmax=vmax)
        i = 0
        tdat = []
        while i < len(idat[1]):
            tdat.append(-1 * idat[1][i])
            i += 1
        plt.scatter(idat[0],tdat,color='black',s=5,vmin=vmin,vmax=vmax)
        p = plt.scatter(sdat[0],sdat[1],c=sdat[2],cmap='rainbow',s=5,vmin=vmin,vmax=vmax)
        tsdat = []
        i = 0
        while i < len(sdat[1]):
            tsdat.append(-1 * sdat[1][i])
            i += 1
        plt.scatter(sdat[0],tsdat,c=sdat[2],cmap='rainbow',s=5,vmin=vmin,vmax=vmax)
        plt.colorbar(p)
        plt.savefig(save)
        print("plotted:",t)
source = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + r'/'
if __name__ == '__main__':
    fig = plt.figure()
    case = True
    t = 0.00
    dt = 0.01
    j = 0
    while(case):
        try:
            sx = []
            sy = []
            sr = []
            ix = []
            iy = []
            inx = []
            iny = []
            spath = source + r'basiliskRuns/' + "skeleton-{:.3f}.dat".format(t)
            ipath = source + r'basiliskRuns/' + "infc-{:.3f}.dat".format(t)
            with open(spath,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                i = 0
                for row in data:
                    sx.append(float(row[0]))
                    sy.append(float(row[1]))
                    sr.append(float(row[2]))
                    i += 1
            with open(ipath,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                i = 0
                for row in data:
                    ix.append(float(row[0]))
                    iy.append(float(row[1]))
                    inx.append(float(row[2]))
                    iny.append(float(row[3]))
                    i += 1
            print("loaded:",t)
            sdat = [sx,sy,sr]
            idat = [ix,iy,inx,iny]
            save = source + r'pScript/2DEvolve/' + "skeleplt-{:03d}.png".format(j)
            t = round(t,2)
            plotmirror(sdat,idat,1,0,save,j)
            #tries to open time step & plot
            j += 1
            t += dt
        except:
            case = False
    print("done")
