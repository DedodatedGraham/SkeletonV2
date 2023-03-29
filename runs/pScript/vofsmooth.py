import matplotlib.pyplot as plt
import sys
import getopt
import os
import csv



def plotmirror(smoothdat,vofdat,save,t,focussmooth,focusvof,mode):
    #We will mirror results here
    #First find focus box
    fig = plt.figure()
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)
    xb = focussmooth[0]
    xt = focussmooth[1]
    yb = focussmooth[2]
    yt = focussmooth[3]
    xb1 = focusvof[0]
    xt1 = focusvof[1]
    yb1 = focusvof[2]
    yt1 = focusvof[3]
    ax1.set_xlim([xb,xt])
    ax1.set_ylim([yb,yt])
    ax2.set_xlim([xb1,xt1])
    ax2.set_ylim([yb1,yt1])
    ax1.set_box_aspect(1)
    ax2.set_box_aspect(1)
    stitle = "smooth:{:.2f}".format(t/100)
    vtitle = "vof:{:.2f}".format(t/100)
    ax1.set_title(stitle)
    ax2.set_title(vtitle)
    if mode == 0:
        #plot smooth
        ax1.scatter(smoothdat[3],smoothdat[4],color='black',s=5)
        p = ax1.scatter(smoothdat[0],smoothdat[1],c=smoothdat[2],cmap='rainbow',s=5)
        plt.colorbar(p,ax=ax1,fraction=0.046,pad=0.04)
        #plot vof
        ax2.scatter(vofdat[3],vofdat[4],color='black',s=5)
        q = ax2.scatter(vofdat[0],vofdat[1],c=vofdat[2],cmap='rainbow',s=5)
        plt.colorbar(q,ax=ax2,fraction=0.046,pad=0.04)
        plt.savefig(save,dpi=1000)
    elif mode == 1:
        #Make quivers
        soriginx = []
        soriginy = []
        snormx = []
        snormy = []
        voriginx = []
        voriginy = []
        vnormx = []
        vnormy = []
        i = 0
        while(i < len(smoothdat[3])):
            soriginx.append(smoothdat[3])
            soriginy.append(smoothdat[4])
            snormx.append(smoothdat[5])
            snormy.append(smoothdat[6])
            i += 1
        i = 0
        while(i < len(vofdat[3])):
            voriginx.append(vofdat[3])
            voriginy.append(vofdat[4])
            vnormx.append(vofdat[5])
            vnormy.append(vofdat[6])
            i += 1
        #plot smooth
        ax1.scatter(smoothdat[3],smoothdat[4],color='black',s=5)
        ax1.quiver(soriginx,soriginy,snormx,snormy,scale_units = 'x',scale = 20)
        #plot vof
        ax2.scatter(vofdat[3],vofdat[4],color='black',s=5)
        ax2.quiver(voriginx,voriginy,vnormx,vnormy, scale_units = 'x',scale = 20)
        plt.savefig(save,dpi=1000)
    print("plotted:",t)



source = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + r'/'
if __name__ == '__main__':
    case = True
    t = 0.00
    dt = 0.01
    j = 0
    focus1 = [0,0,0,0]
    fleng1 = 0;
    focus2 = [0,0,0,0]
    fleng2 = 0;
    while(case):
        #try:
        smoothx = []
        smoothy = []
        smoothr = []
        smoothix = []
        smoothiy = []
        smoothinx = []
        smoothiny = []
        vofx = []
        vofy = []
        vofr = []
        vofix = []
        vofiy = []
        vofinx = []
        vofiny = []
        smoothpath = source + r'superRuns/' + "smoothskeleton-{:.3f}.dat".format(t)
        vofpath = source + r'superRuns/' + "vofskeleton-{:.3f}.dat".format(t)
        intpath = source + r'superRuns/' + "intdata-{:.3f}.dat".format(t)
        #ipath = source + r'superRuns/' + "infc-{:.3f}.dat".format(t)
        with open(smoothpath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                smoothx.append(float(row[0]))
                smoothy.append(float(row[1]))
                smoothr.append(float(row[2]))
                i += 1
        with open(vofpath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                vofx.append(float(row[0]))
                vofy.append(float(row[1]))
                vofr.append(float(row[2]))
                i += 1
        with open(intpath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                vofix.append(float(row[0]))
                vofiy.append(float(row[1]))
                vofinx.append(float(row[2]))
                vofiny.append(float(row[3]))
                smoothix.append(float(row[4]))
                smoothiy.append(float(row[5]))
                smoothinx.append(float(row[6]))
                smoothiny.append(float(row[7]))
                i += 1
        print("loaded:",t)
        smoothdat = [smoothx,smoothy,smoothr,smoothix,smoothiy,smoothinx,smoothiny]
        vofdat = [vofx,vofy,vofr,vofix,vofiy,vofinx,vofiny]
        #save = source + r'pScript/2DEvolve/' + "skeleplt-{:03d}.png".format(j)
        save = r'2DEvolve/' + "svofplt-{:03d}.png".format(j)
        t = round(t,2)
        #find focus
        focus1[0] = min(smoothdat[3])
        focus1[1] = max(smoothdat[3])
        focus1[2] = min(smoothdat[4])
        focus1[3] = max(smoothdat[4])
        fx1 = focus1[1]-focus1[0]#lx
        fy1 = focus1[3]-focus1[2]#ly
        mx1 = (focus1[1]+focus1[0]) / 2 #mid point x
        my1 = (focus1[3]+focus1[2]) / 2 #mid point y
        #Make const Square bounds; fleng is the square width
        if j == 0:
            if fx1 > fy1:
                fleng1 = fx1*2
            else:
                fleng1 = fy1*2
        #Adjust square as moves
        if fx1 < fleng1:
            focus1[0] = mx1 - fleng1/2
            focus1[1] = mx1 + fleng1/2
        if fy1 < fleng1:
            focus1[2] = my1 - fleng1/2
            focus1[3] = my1 + fleng1/2
        #find focus 2
        focus2[0] = min(vofdat[3])
        focus2[1] = max(vofdat[3])
        focus2[2] = min(vofdat[4])
        focus2[3] = max(vofdat[4])
        fx2 = focus2[1]-focus2[0]#lx
        fy2 = focus2[3]-focus2[2]#ly
        mx2 = (focus2[1]+focus2[0]) / 2 #mid point x
        my2 = (focus2[3]+focus2[2]) / 2 #mid point y
        #Make const Square bounds; fleng is the square width
        if j == 0:
            if fx2 > fy2:
                fleng2 = fx2*2
            else:
                fleng2 = fy2*2
        #Adjust square as moves
        if fx2 < fleng2:
            focus2[0] = mx2 - fleng2/2
            focus2[1] = mx2 + fleng2/2
        if fy2 < fleng2:
            focus2[2] = my2 - fleng2/2
            focus2[3] = my2 + fleng2/2
        #Plot results
        plotmirror(smoothdat,vofdat,save,j,focus1,focus2,1)
        #tries to open time step & plot
        j += 1
        t += dt
        #except:
        #    case = False
    print("done")
