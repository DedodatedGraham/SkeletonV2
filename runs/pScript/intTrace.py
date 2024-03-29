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



def plotmirror(fig,sdat,idat,ndat,bdat,mirroraxi,mirrorval,save,t,focus):
    #We will mirror results here
    #First find focus box
    plt.clf()
    xb = focus[0]
    xt = focus[1]
    yb = focus[2]
    yt = focus[3]
    ax = fig.add_axes([0.1,0.1,0.9,0.9])
    ax.set_xlim(xb,xt)
    ax.set_ylim(yb,yt)

    i = 0
    plotboxx = []
    plotboxy = []
    rplotboxx = []
    rplotboxy = []
    dplotboxx = []
    dplotboxy = []
    dplot = []
    while i < len(bdat[0]):
        plotboxx.append([bdat[0][i],bdat[0][i],bdat[2][i],bdat[2][i],bdat[0][i]])
        plotboxy.append([bdat[1][i],bdat[3][i],bdat[3][i],bdat[1][i],bdat[1][i]])
    #    rplotboxx.append([bdat[4][i],bdat[4][i],bdat[6][i],bdat[6][i],bdat[4][i]])
    #    rplotboxy.append([bdat[5][i],bdat[7][i],bdat[7][i],bdat[5][i],bdat[5][i]])
    #    dplotboxx.append([bdat[0][i],bdat[2][i]])
    #    dplotboxy.append([bdat[1][i],bdat[3][i]])
    #    dplot.append(bdat[8][i])
        i+=1
    #cmap = matplotlib.cm.get_cmap('rainbow')
    #dvmax = max(dplot)
    #img = plt.imshow(np.array([[0,1]]), cmap="rainbow",vmin=0,vmax=dvmax)
    #img.set_visible(False)

    #i = 0
    #patches = []
    #clist = []
    #while i < len(plotboxx):
    #    w = dplotboxx[i][1] - dplotboxx[i][0]
    #    h = dplotboxy[i][1] - dplotboxy[i][0]
    #    rgba = cmap(dplot[i]/dvmax)
    #    r = Rectangle((dplotboxx[i][0],dplotboxy[i][0]),w,h)
    #    patches.append(r)
    #    clist.append(rgba)
    #    i += 1
    #ourcmap = ListedColormap(clist)
    #patches_collection = PatchCollection(patches,cmap=ourcmap)
    #patches_collection.set_array(np.arange(len(patches)))
    #ax.add_collection(patches_collection)

    i = 0
    while i < len(plotboxx):
        ax.plot(plotboxx[i],plotboxy[i],c='black',linewidth=0.6)
        i += 1
    i = 0
    while i < len(plotboxx):
        #plt.plot(rplotboxx[i],rplotboxy[i],c='red')
        i += 1

    title = "Time:{:.2f}".format(t/100)
    ax.set_title(title)
    vmin = 0
    vmax = max(ndat[2])
    plt.rcParams['figure.dpi'] = 1000
    if mirroraxi == 0:
        ax.scatter(idat[0],idat[1],c='black',s=2)
        ax.scatter(sdat[0],sdat[1],c='black',s=2)
        #p = ax.scatter(ndat[0],ndat[1],c=ndat[3],cmap='rainbow',s=5)
        p = ax.scatter(ndat[0],ndat[1],c='red',s=5)
        #plt.scatter(idat[0],idat[1],c=idat[4],cmap='rainbow',s=5)
        #p = plt.scatter(sdat[0],sdat[1],c=sdat[2],cmap='rainbow',s=5)
        #p = plt.scatter(ndat[0],ndat[1],c=ndat[2],cmap='rainbow',s=5,vmin=vmin,vmax=vmax)
        i = 0
        #plt.colorbar(p)
        #plt.show()
        #fig.show()
        plt.colorbar(p)
        plt.savefig(save,dpi=1000)
        print("plotted:",t)



source = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + r'/'
if __name__ == '__main__':
    fig = plt.figure()
    case = True
    t = 0.0
    dt = 0.01
    j = 0
    focus = [0,0,0,0]
    fleng = 0;
    while(case):
        #try:
        sx = []
        sy = []
        sr = []
        ix = []
        iy = []
        inx = []
        iny = []
        col = []
        npx = []
        npy = []
        npr = []
        bnx = []
        bny = []
        bpx = []
        bpy = []
        rbnx = []
        rbny = []
        rbpx = []
        rbpy = []
        bd = []
        hasn = []
        alpha = []
        drow = 0
        dcol = 0
        sm = []
        angle = np.tan(30 * 3.141592 / 180)
        #print("angle:",angle)
        spath = source + r'dbasiliskRuns/' + "reducedskeleton-{:.3f}.dat".format(t)
        with open(spath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                #if (float(row[7])) > angle:
                sx.append(float(row[0]))
                sy.append(float(row[1]))
                sr.append(float(row[2]))
                alpha.append(float(row[3]))
                #else:
                #    sx.append(float(row[0]))
                #    sy.append(float(row[1]))
                #    sr.append(0.0)
                #ix.append(float(row[3]))
                #iy.append(float(row[4]))
                #inx.append(float(row[5]))
                #iny.append(float(row[6]))
                i += 1
        npath = source + r'dbasiliskRuns/' + "nodeDat-{:.3f}.dat".format(t)
        with open(npath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                npx.append(float(row[0]))
                npy.append(float(row[1]))
                npr.append(float(row[2]))
                sm.append(int(row[3]))
                i += 1
        bpath = source + r'dbasiliskRuns/' + "boxDat-{:.3f}.dat".format(t)
        #ipath = source + r'superRuns/' + "infc-{:.3f}.dat".format(t)
        with open(bpath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                bnx.append(float(row[0]))
                bny.append(float(row[1]))
                bpx.append(float(row[2]))
                bpy.append(float(row[3]))
                rbnx.append(float(row[4]))
                rbny.append(float(row[5]))
                rbpx.append(float(row[6]))
                rbpy.append(float(row[7]))
                bd.append(float(row[8]))
                drow = int(row[9])
                dcol = int(row[10])
                hasn.append(int(row[11]))
                i += 1
        ipath = source + r'dbasiliskRuns/' + "intdata-{:.3f}.dat".format(t)
        with open(ipath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                ix.append(float(row[0]))
                iy.append(float(row[1]))
                inx.append(float(row[2]))
                iny.append(float(row[3]))
                i += 1
        #print(alpha)
        print("loaded:",t)
        sdat = [sx,sy,sr,alpha]
        idat = [ix,iy,inx,iny]
        ndat = [npx,npy,npr,sm]
        #ndat = []
        bdat = [bnx,bny,bpx,bpy,rbnx,rbny,rbpx,rbpy,bd,drow,dcol,hasn]
        #bdat = []
        #save = source + r'pScript/2DEvolve/' + "skeleplt-{:03d}.png".format(j)
        save = r'2DEvolve/' + "skele2intplt-{:03d}.png".format(j)
        t = round(t,2)
        #find focus
        focus[0] = min(idat[0])
        focus[1] = max(idat[0])
        focus[2] = min(idat[1])
        focus[3] = max(idat[1])
        fx = focus[1]-focus[0]#lx
        fy = focus[3]-focus[2]#ly
        mx = (focus[1]+focus[0]) / 2 #mid point x
        my = (focus[3]+focus[2]) / 2 #mid point y
        #Make const Square bounds; fleng is the square width
        if j == 0:
            if fx > fy:
                fleng = fx*2
            else:
                fleng = fy*2
        #Adjust square as moves
        if fx < fleng:
            focus[0] = mx - fleng/2
            focus[1] = mx + fleng/2
        if fy < fleng:
            focus[2] = my - fleng/2
            focus[3] = my + fleng/2
        #Plot results
        if len(sdat[0]) != 0:
            plotmirror(fig,sdat,idat,ndat,bdat,0,0,save,j,focus)
        #tries to open time step & plot
        j += 1
        t += dt
        #except:
        #    case = False
    print("done")
