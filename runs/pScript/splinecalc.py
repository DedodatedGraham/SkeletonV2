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


def generate_nonuniform_knot_vector(n, degree):
    num_knots = n + degree + 1
    knot_vector = [0] * num_knots

    for i in range(num_knots):
        if i < degree:
            knot_vector[i] = 0
        elif i <= n:
            knot_vector[i] = i - degree + 1
        else:
            knot_vector[i] = n - degree + 2

    return knot_vector

def basis_function(i, degree, param, knots):
    # Check if the degree is 0
    if degree == 0:
        return 1.0 if knots[i] <= param < knots[i + 1] else 0.0

    # Recursive calculation of the B-spline basis function
    denominator1 = knots[i + degree - 1] - knots[i]  # Updated indexing here
    denominator2 = 0.0

    if i + degree < len(knots):
        denominator2 = knots[i + degree] - knots[i + 1]

    term1 = 0.0
    term2 = 0.0

    if denominator1 != 0:
        term1 = ((param - knots[i]) / denominator1) * basis_function(i, degree - 1, param, knots)
    if denominator2 != 0:
        term2 = ((knots[i + degree] - param) / denominator2) * basis_function(i + 1, degree - 1, param, knots)

    return term1 + term2

def calculate_bspline_point(control_points, param,knots):
    n = len(control_points) - 1  # Number of control points
    degree = n  # Degree of the B-spline curve
    # Initialize the point coordinates
    x = 0.0
    y = 0.0
    r = 0.0

    # Calculate the basis functions and weighted control point contributions
    for i in range(n + 1):
        basis = basis_function(i, degree, param,knots)
        x += control_points[i][0] * basis
        y += control_points[i][1] * basis
        r += control_points[i][2] * basis

    return x, y, r


def plotmirror(fig,scdat,idat,ndat,bdat,condat,sbdat,mirroraxi,mirrorval,save,t,focus):
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
    #while i < len(scdat[0]):
    #    #We decompose needed
    #    j = 0
    #    idpass = False
    #    while j < len(ids):
    #        if(ids[j] == scdat[5][i]):
    #            idpass = True
    #            scidx[j].append(scdat[0][i])
    #            scidy[j].append(scdat[1][i])
    #            scidw[j].append(scdat[2][i])
    #            scidh[j].append(scdat[3][i])
    #            scidd[j].append(scdat[4][i])
    #        j += 1
    #    if not idpass:
    #        ids.append(scdat[5][i])
    #        scidx.append([])
    #        scidy.append([])
    #        scidw.append([])
    #        scidh.append([])
    #        scidd.append([])
    #        scidx[len(ids) - 1].append(scdat[0][i])
    #        scidy[len(ids) - 1].append(scdat[1][i])
    #        scidw[len(ids) - 1].append(scdat[2][i])
    #        scidh[len(ids) - 1].append(scdat[3][i])
    #        scidd[len(ids) - 1].append(scdat[4][i])
    #    i += 1
    #cmap = matplotlib.cm.get_cmap('rainbow')
    ##dvmax = max(scdat[4])
    #dvmax = len(scidx)
    #if(dvmax == 0):
    #    dvmax = 1
    ##img = plt.imshow(np.array([[0,1]]), cmap="rainbow",vmin=0,vmax=dvmax)
    ##img.set_visible(False)

    #i = 0
    #patches = []
    #clist = []
    #while i < len(scidx):
    #    j = 0
    #    rgba = cmap(i/dvmax)
    #    while j < len(scidx[i]):
    #        #rgba = cmap(scidd[i][j]/dvmax)
    #        r = Rectangle((scidx[i][j],scidy[i][j]),scidw[i][j],scidh[i][j])
    #        patches.append(r)
    #        clist.append(rgba)
    #        j += 1
    #    i += 1
    #ourcmap = ListedColormap(clist)
    #patches_collection = PatchCollection(patches,cmap=ourcmap)
    #patches_collection.set_array(np.arange(len(patches)))
    #ax.add_collection(patches_collection)

    #i = 0
    #while i < len(plotboxx):
    #    ax.plot(plotboxx[i],plotboxy[i],c='black',linewidth=0.6)
    #    i += 1
    #i = 0
    #while i < len(plotboxx):
    #    #plt.plot(rplotboxx[i],rplotboxy[i],c='red')
    #    i += 1

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

        if(len(sbdat[0]) > 0):
            #n-1 method
            mr = 0
            i = 0
            while i < len(sbdat[0]):
                if mr < sbdat[2][i]:
                    mr = sbdat[2][i]
                if mr < sbdat[5][i]:
                    mr = sbdat[5][i]
                i += 1
            i = 0
            while i < len(sbdat[0]):
                px = []
                py = []
                pr = []
                #start
                px.append(sbdat[0][i])
                py.append(sbdat[1][i])
                pr.append(sbdat[2][i])
                #calcpoints
                j = 1
                dx = (sbdat[3][i]-sbdat[0][i])/100
                dy = (sbdat[4][i]-sbdat[1][i])/100
                dr = (sbdat[5][i]-sbdat[2][i])/100
                startx = sbdat[0][i]
                starty = sbdat[1][i]
                startr = sbdat[2][i]
                while j < 100:
                    cx = startx+j*dx
                    cy = starty+j*dy
                    cr = startr+j*dr
                    px.append(cx)
                    #if(len(sbdat[7]) == 0):
                    py.append(cy)
                    pr.append(cr)
                    #else:
                    #    py.append(cx*sbdat[6][i] + pow(cx,2)*sbdat[7][i])
                    j += 1
                #end
                px.append(sbdat[3][i])
                py.append(sbdat[4][i])
                pr.append(sbdat[5][i])
                j = 0
                while j < len(px):
                    drawcirc = plt.Circle((px[j],py[j]),pr[j])
                    ax.add_artist(drawcirc)
                    j += 1
                j = 0
                while j < len(px)-1:
                    tpx = []
                    tpy = []
                    tpx.append(px[j])
                    tpy.append(py[j])
                    tpx.append(px[j+1])
                    tpy.append(py[j+1])
                    color = plt.cm.get_cmap('hsv')(((pr[j] + pr[j+1])/2)/mr)
                    ax.plot(tpx,tpy,linewidth=3,color=color)
                    j += 1
                i += 1
            #n-2 method
            #mr = 0
            #i = 0
            #while i < len(sbdat[0]):
            #    j = 0
            #    while j < len(sbdat[0][i]):
            #        if sbdat[2][i][j] > mr:
            #            mr = sbdat[2][i][j]
            #        j += 1
            #    i += 1
            #i = 0
            #while i < len(sbdat[0]):
            #    sx = []
            #    sy = []
            #    sr = []
            #    j = 0
            #    control_points = []
            #    while j < len(sbdat[0][i]):
            #        control_points.append([sbdat[0][i][j],sbdat[1][i][j],sbdat[2][i][j]])
            #        j += 1
            #    degree = len(control_points) - 1
            #    knots = generate_nonuniform_knot_vector(degree, degree)
            #    ts = np.linspace(knots[len(knots)-1], knots[0], num=100)
            #    for param in ts:
            #        spline_point = calculate_bspline_point(control_points, param, knots)
            #        sx.append(spline_point[0])
            #        sy.append(spline_point[1])
            #        sr.append(spline_point[2])
            #    j = 0
            #    print()
            #    print('sx',sx)
            #    print('sy',sy)
            #    print('sr',sr)
            #    print()
            #    while j < len(sx) - 1:
            #        tpx = []
            #        tpy = []
            #        tpx.append(sx[j])
            #        tpy.append(sy[j])
            #        tpx.append(sx[j+1])
            #        tpy.append(sy[j+1])
            #        color = plt.cm.get_cmap('hsv')(((sr[j] + sr[j+1])/2)/mr)
            #        ax.plot(tpx,tpy,linewidth=3,color=color)
            #        j += 1
            #    i += 1
        ax.scatter(idat[0],idat[1],c='black',s=2)
        #print("scattered")
        #if(len(sbdat[0]) > 0):
        #    i = 0
        #    while i < len(sbdat[0]):
        #        ax.scatter(sbdat[0][i],sbdat[1][i],s=2)
        #        i += 1
        ###if(len(sbdat[3]) > 0):
        ###    points = {}
        ###    for i in range(len(sbdat[3])):
        ###        curr_x = sbdat[3][i]
        ###        curr_y = sbdat[4][i]
        ###        angle = sbdat[5][i]

        ###        if (curr_x, curr_y) not in points:
        ###            points[(curr_x, curr_y)] = []

        ###        points[(curr_x, curr_y)].append(angle)

        ###    # Plot lines for each group of points
        ###    i = 0
        ###    for coords, angles in points.items():
        ###        curr_x, curr_y = coords
        ###        color = plt.cm.get_cmap('hsv')(i/len(points))
        ###        px = []
        ###        py = []
        ###        px.append(curr_x)
        ###        py.append(curr_y)
        ###        for angle in angles:
        ###            px.append(curr_x + 10 * math.cos(angle))
        ###            py.append(curr_y + 10 * math.sin(angle))
        ###            px.append(curr_x)
        ###            py.append(curr_y)
        ###        ax.plot(px, py, linewidth=3, color=color)
        ###        i += 1
        #ax.scatter(ndat[0],ndat[1],c=ndat[3],cmap='Set1',s=5)
        #ax.scatter(ndat[0],ndat[1],c='black',s=5)
        #plt.scatter(idat[0],idat[1],c=idat[4],cmap='rainbow',s=5)
        #p = plt.scatter(sdat[0],sdat[1],c=sdat[2],cmap='rainbow',s=5)
        #p = plt.scatter(ndat[0],ndat[1],c=ndat[2],cmap='rainbow',s=5,vmin=vmin,vmax=vmax)
        #Finally we draw the lines of each id
        #i = 0
        #while i < len(condat[0]):
        #    px = [condat[1][condat[0][i][0]][2],condat[1][condat[0][i][1]][2]]
        #    py = [condat[1][condat[0][i][0]][3],condat[1][condat[0][i][1]][3]]
        #    plt.plot(px,py,c='blue',linewidth=3);
        #    i += 1
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

        #spline branch data

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
        #Load in node connections
        con = []
        conpath = source + r'dbasiliskRuns/' + "connectionDat-{:.3f}.dat".format(t)
        with open(conpath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                con.append([int(row[0]),int(row[1])])
                i += 1
        conid = []
        conidpath = source + r'dbasiliskRuns/' + "connectionidDat-{:.3f}.dat".format(t)
        with open(conidpath,'r') as csvfile:
            data  = csv.reader(csvfile,delimiter = ' ')
            i = 0
            for row in data:
                conid.append([int(row[0]),int(row[1]),float(row[2]),float(row[3])])
                i += 1
        sbpath = ""
        #n-1 data
        sn0x = []
        sn0y = []
        sn0r = []
        sn1x = []
        sn1y = []
        sn1r = []
        sa0 = []
        sa1 = []
        try:
            sbpath = source + r'dbasiliskRuns/' + "splineBranchDat-{:.3f}.dat".format(t)
            with open(sbpath,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                i = 0
                for row in data:
                    sn0x.append(float(row[0]))
                    sn0y.append(float(row[1]))
                    sn0r.append(float(row[2]))
                    sn1x.append(float(row[3]))
                    sn1y.append(float(row[4]))
                    sn1r.append(float(row[5]))
                    sa0.append(float(row[6]))
                    if(len(row) > 7):
                        sa1.append(float(row[7]))
                    i += 1
        except:
            print("no existing file:",sbpath)
        #n-2 data input
        #n = 0
        #conx = []
        #cony = []
        #conr = []
        #try:
        #    sbpath = source + r'dbasiliskRuns/' + "splineBranchDat-{:.3f}.dat".format(t)
        #    with open(sbpath,'r') as csvfile:
        #        data  = csv.reader(csvfile,delimiter = ' ')
        #        i = 0
        #        for row in data:
        #            print('ready')
        #            conx.append([])
        #            cony.append([])
        #            conr.append([])
        #            if i == 0:
        #                n = int(row[0])
        #            q = 0
        #            c = 3
        #            while q < n + 1:
        #                conx[i].append(float(row[c*q+1]))
        #                cony[i].append(float(row[c*q+2]))
        #                conr[i].append(float(row[c*q+3]))
        #                q += 1
        #            i += 1
        #except:
        #    print("no existing file:",sbpath)
        angpath = source + r'dbasiliskRuns/' + "angDat-{:.3f}.dat".format(t)
        ax = []
        ay = []
        aa = []
        try:
            with open(angpath,'r') as csvfile:
                data  = csv.reader(csvfile,delimiter = ' ')
                i = 0
                for row in data:
                    ax.append(float(row[0]))
                    ay.append(float(row[1]))
                    aa.append(float(row[2]))
                    i += 1
        except:
            print("couldnt load ang")
        print("loaded:",t)
        scdat = [scx,scy,scw,sch,scd,scid]
        idat = [ix,iy,inx,iny]
        ndat = [npx,npy,npr,sm]
        bdat = [bnx,bny,bpx,bpy,rbnx,rbny,rbpx,rbpy,bd,drow,dcol,hasn]
        condat = [con,conid]
        sbdat = [sn0x,sn0y,sn0r,sn1x,sn1y,sn1r,sa0,sa1]
        #sbdat = [conx,cony,conr]
        print("assigned",t)
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
        print("sent",t)
        plotmirror(fig,scdat,idat,ndat,bdat,condat,sbdat,0,0,save,j,focus)
        #tries to open time step & plot
        j += 1
        t += dt
        #except:
        #    case = False
    print("done")
