import matplotlib.pyplot as plt
import sys
import getopt
import os
import csv



def plotmirror(sdat,idat,ndat,bdat,mirroraxi,mirrorval,save,t,focus):
    #We will mirror results here
    #First find focus box
    plt.clf()
    xb = focus[0]
    xt = focus[1]
    yb = focus[2]
    yt = focus[3]
    plt.xlim(xb,xt)
    plt.ylim(yb,yt)

    i = 0
    plotboxx = []
    plotboxy = []
    while i < len(bdat[0]):
        plotboxx.append([bdat[0][i],bdat[0][i],bdat[2][i],bdat[2][i],bdat[0][i]])
        plotboxy.append([bdat[1][i],bdat[3][i],bdat[3][i],bdat[1][i],bdat[1][i]])
        i+=1
    i = 0
    while i < len(plotboxx):
        plt.plot(plotboxx[i],plotboxy[i],c='black')
        i += 1
    title = "Time:{:.2f}".format(t/100)
    plt.title(title)
    vmin = 0
    vmax = max(sdat[2])
    plt.rcParams['figure.dpi'] = 1000
    if mirroraxi == 0:
        plt.scatter(idat[0],idat[1],c='black',s=5)
        plt.scatter(sdat[0],sdat[1],c='black',s=5)
        #p = plt.scatter(sdat[0],sdat[1],c=sdat[2],cmap='rainbow',s=5)
        p = plt.scatter(ndat[0],ndat[1],c=ndat[2],cmap='rainbow',s=5,vmin=vmin,vmax=vmax)
        i = 0
        plt.colorbar(p)
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
        try:
            if t == 1.62:
                break
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
            spath = source + r'dbasiliskRuns/' + "smoothskeleton-{:.3f}.dat".format(t)
            #ipath = source + r'superRuns/' + "infc-{:.3f}.dat".format(t)
            with open(spath,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                i = 0
                for row in data:
                    sx.append(float(row[0]))
                    sy.append(float(row[1]))
                    sr.append(float(row[2]))
                    ix.append(float(row[3]))
                    iy.append(float(row[4]))
                    inx.append(float(row[5]))
                    iny.append(float(row[6]))
                    col.append(float(row[7]))
                    i += 1
            npath = source + r'dbasiliskRuns/' + "nodePoint-{:.3f}.dat".format(t)
            #ipath = source + r'superRuns/' + "infc-{:.3f}.dat".format(t)
            with open(npath,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                i = 0
                for row in data:
                    npx.append(float(row[0]))
                    npy.append(float(row[1]))
                    npr.append(float(row[2]))
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
                    i += 1
            print("loaded:",t)
            sdat = [sx,sy,sr,col]
            idat = [ix,iy,inx,iny,col]
            ndat = [npx,npy,npr]
            bdat = [bnx,bny,bpx,bpy]
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
            plotmirror(sdat,idat,ndat,bdat,0,0,save,j,focus)
            #tries to open time step & plot
            j += 1
            t += dt
        except:
            case = False
    print("done")
