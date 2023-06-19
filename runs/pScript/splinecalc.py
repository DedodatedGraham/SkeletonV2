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



def plotmirror(fig,scdat,idat,ndat,bdat,mirroraxi,mirrorval,save,t,focus):
    #We will mirror results here
    #First find focus box
    plt.clf()
    xb = focus[0]
    xt = focus[1]
    yb = focus[2]
    yt = focus[3]
    ax = fig.add_axes([0.1,0.1,0.85,0.85])
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
    i = 0
    ids = []
    scidx = []
    scidy = []
    scidw = []
    scidh = []
    scidd = []
    while i < len(scdat[0]):
        #We decompose needed
        j = 0
        idpass = False
        while j < len(ids):
            if(ids[j] == scdat[5][i]):
                idpass = True
                scidx[j].append(scdat[0][i])
                scidy[j].append(scdat[1][i])
                scidw[j].append(scdat[2][i])
                scidh[j].append(scdat[3][i])
                scidd[j].append(scdat[4][i])
            j += 1
        if not idpass:
            ids.append(scdat[5][i])
            scidx.append([])
            scidy.append([])
            scidw.append([])
            scidh.append([])
            scidd.append([])
            scidx[len(ids) - 1].append(scdat[0][i])
            scidy[len(ids) - 1].append(scdat[1][i])
            scidw[len(ids) - 1].append(scdat[2][i])
            scidh[len(ids) - 1].append(scdat[3][i])
            scidd[len(ids) - 1].append(scdat[4][i])
        i += 1
    cmap = matplotlib.cm.get_cmap('rainbow')
    #dvmax = max(scdat[4])
    dvmax = len(scidx)
    if(dvmax == 0):
        dvmax = 1
    #img = plt.imshow(np.array([[0,1]]), cmap="rainbow",vmin=0,vmax=dvmax)
    #img.set_visible(False)

    i = 0
    patches = []
    clist = []
    while i < len(scidx):
        j = 0
        rgba = cmap(i/dvmax)
        while j < len(scidx[i]):
            #rgba = cmap(scidd[i][j]/dvmax)
            r = Rectangle((scidx[i][j],scidy[i][j]),scidw[i][j],scidh[i][j])
            patches.append(r)
            clist.append(rgba)
            j += 1
        i += 1
    ourcmap = ListedColormap(clist)
    patches_collection = PatchCollection(patches,cmap=ourcmap)
    patches_collection.set_array(np.arange(len(patches)))
    ax.add_collection(patches_collection)

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
    ##ndat[0].append(10000000)
    ##ndat[1].append(10000000)
    ##ndat[3].append(2)
    ##ndat[0].append(10000000)
    ##ndat[1].append(10000000)
    ##ndat[3].append(0)

    if mirroraxi == 0:
        ax.scatter(idat[0],idat[1],c='black',s=2)
        #ax.scatter(sdat[0],sdat[1],c='black',s=2)
        ax.scatter(ndat[0],ndat[1],c=ndat[3],cmap='Set1',s=5)
        #ax.scatter(ndat[0],ndat[1],c='black',s=5)
        #plt.scatter(idat[0],idat[1],c=idat[4],cmap='rainbow',s=5)
        #p = plt.scatter(sdat[0],sdat[1],c=sdat[2],cmap='rainbow',s=5)
        #p = plt.scatter(ndat[0],ndat[1],c=ndat[2],cmap='rainbow',s=5,vmin=vmin,vmax=vmax)
        i = 0
        #plt.colorbar(p)
        #plt.show()
        #fig.show()
        #plt.colorbar(p)
        plt.savefig(save,dpi=1000,figsize=[12.8,9.6])
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
        drow = 0
        dcol = 0
        sm = []
        #our relavant close data
        scx = []
        scy = []
        scw = []
        sch = []
        scd = []
        scid = []

        angle = np.tan(30 * 3.141592 / 180)
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
        scpath = source + r'dbasiliskRuns/' + "splinecalcDat-{:.3f}.dat".format(t)
        with open(scpath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                scx.append(float(row[0]))
                scy.append(float(row[1]))
                scw.append(float(row[2]))
                sch.append(float(row[3]))
                scd.append(int(row[4]))
                scid.append([int(row[5]),int(row[6])])
                i += 1
        print("loaded:",t)
        scdat = [scx,scy,scw,sch,scd,scid]
        idat = [ix,iy,inx,iny]
        ndat = [npx,npy,npr,sm]
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
        #if len(sdat[0]) != 0:
        plotmirror(fig,scdat,idat,ndat,bdat,0,0,save,j,focus)
        #tries to open time step & plot
        j += 1
        t += dt
        #except:
        #    case = False
    print("done")
