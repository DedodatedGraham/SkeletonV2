import matplotlib.pyplot as plt
import sys
import getopt
import os
import csv



def plotmirror(idat,mirroraxi,mirrorval,save,t,focus):
    #We will mirror results here
    #First find focus box
    plt.clf()
    xb = focus[0]
    xt = focus[1]
    yb = focus[2]
    yt = focus[3]
    plt.xlim(xb,xt)
    plt.ylim(yb,yt)
    title = "Time:{:.2f}".format(t/100)
    plt.title(title)
    vmin = 0
    vmax = 0.5
    plt.rcParams['figure.dpi'] = 1000
    if mirroraxi == 0:
        plt.scatter(idat[0],idat[1],color='Red',s=5,vmin=vmin,vmax=vmax)
        plt.scatter(idat[4],idat[5],color='Blue',s=5,vmin=vmin,vmax=vmax)
        plt.savefig(save,dpi=1000)
        print("plotted:",t)



source = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + r'/'
if __name__ == '__main__':
    fig = plt.figure()
    case = True
    t = 0.00
    dt = 0.01
    j = 0
    focus = [0,0,0,0]
    fleng = 0;
    while(case):
        #try:
        if t == 1.62:
            break
        sx = []
        sy = []
        snx = []
        sny = []
        ix = []
        iy = []
        inx = []
        iny = []
        ipath = source + r'superRuns/' + "intdata-{:.3f}.dat".format(t)
        with open(ipath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                ix.append(float(row[0]))
                iy.append(float(row[1]))
                inx.append(float(row[2]))
                iny.append(float(row[3]))
                sx.append(float(row[4]))
                sy.append(float(row[5]))
                snx.append(float(row[6]))
                sny.append(float(row[7]))
                i += 1
        print("loaded:",t)
        idat = [ix,iy,inx,iny,sx,sy,snx,sny]
        #save = source + r'pScript/2DEvolve/' + "skeleplt-{:03d}.png".format(j)
        save = r'2DEvolve/' + "compplt-{:03d}.png".format(j)
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
        plotmirror(idat,0,0,save,j,focus)
        #tries to open time step & plot
        j += 1
        t += dt
        #except:
        #    case = False
    print("done")
