import numpy as np
import matplotlib.pyplot as plt
import scipy

def d(a):
    return a[1:]-a[:-1]

case = 'WD1.0LT5.0'
filename = 'Data/{case}/domain_data.csv'.format(case=case)
GENERATE = False
RUN = True

T = 0.01

if GENERATE:
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    data = data[np.where(data[:,2]<0.0)[0]]
    data[:,0] /= T
    data[:,1] /= T
    data[:,1] = data[:,1]-0.5*np.max(data[:,1])
    x = np.sort(np.unique(data[:,0]))
    y = np.sort(np.unique(data[:,1]))

    I,J = len(x), len(y)
    X,Y = np.meshgrid(x, y, indexing='ij')
    U = np.nan*(X)
    V = np.nan*(X)
    P = np.nan*(X)
    for j in range(J):
        print(j)
        for i in range(I):
            x0,y0 = X[i,j], Y[i,j]
            dx,dy = x0-data[:,0], y0-data[:,1]
            r2 = dx*dx+dy*dy
            idx = np.argmin(r2)
            row = data[idx,:]
            U[i,j],V[i,j],P[i,j] = row[3],row[4],row[12]

    np.save('Data/{case}/X'.format(case=case), X)
    np.save('Data/{case}/Y'.format(case=case), Y)
    np.save('Data/{case}/U'.format(case=case), U)
    np.save('Data/{case}/V'.format(case=case), V)
    np.save('Data/{case}/P'.format(case=case), P)

if RUN:
    T = 1.0
    fw = float(case[2:5])
    fL = float(case[7:10])
    w = T*fw
    L = T*fL
    xmin = 0.0
    xmax = 50.0*T
    xLE = 0.2*xmax
    xTE = xLE+L
    yt = 0.5*T
    yb = -0.5*T
    yT = yt + w
    yB = yb - w
    H = T+2.0*w

    xlim = [xTE-0.25*L, xTE+1.75*L]
    ylim = [yB-T/4.0, yT+T/4.0]

    X = np.load('Data/{case}/X.npy'.format(case=case))
    Y = np.load('Data/{case}/Y.npy'.format(case=case))
    U = np.load('Data/{case}/U.npy'.format(case=case))
    V = np.load('Data/{case}/V.npy'.format(case=case))
    P = np.load('Data/{case}/P.npy'.format(case=case))
    UV_mag = np.sqrt(U*U+V*V)
    plate_x = [xLE, xTE, xTE, xLE]
    plate_y = [yb, yb, yt, yt]
    wallT_x = [xmin, xmax, xmax, xmin]
    wallT_y = [yT, yT, yT+10.0, yT+10.0]
    wallB_x = [xmin, xmax, xmax, xmin]
    wallB_y = [yB-10, yB-10, yB, yB]
    plate_fc = (0.95, 0.95, 0.95)
    plate_ec = (0.1, 0.1, 0.1)
    plate_hc = (0.9, 0.9, 0.9)
    plate_h_lw =  0.8
    plate_lw = 1.2

    color_c=(1.0, 0.0, 0.0)
    color_nc=(0.8, 0.8, 0.8)
    lw_c = 0.8
    lw_nc = 0.6

    x = X[:,0]
    y = Y[0,:]


    def run(init_x,init_y,xxA,xxB,yyA,yyB,dt0=3e-2,max_delta=T/2.0,max_k=10000,c0=1.04):
        M = len(init_x)
        for m in range(M):
            x0 = init_x[m]
            y0 = init_y[m]

            xx = [x0]
            yy = [y0]
            delta = 0.0
            k = 0
            flag = False
            frwd = np.random.rand()>=0.5
            
            while (k<=max_k and delta<max_delta and not flag):
                i1 = np.argmin(np.abs(x0-x))
                i2 = int(i1 + np.sign(x0-x[i1]))
                j1 = np.argmin(np.abs(y0-y))
                j2 = int(j1 + np.sign(y0-y[j1]))
                ia,ib = max(0,min(i1,i2)),min(max(i1,i2),len(x)-1)
                ja,jb = max(0,min(j1,j2)),min(max(j1,j2),len(y)-1)
                xa = x[ia]; xb = x[ib]
                ya = y[ja]; yb = y[jb]
                if (xb!=xa): 
                    kx = (x0-xa)/(xb-xa)
                else:
                    kx = 1.0
                if (yb!=ya): 
                    ky = (y0-ya)/(yb-ya)
                else:
                    ky = 1.0
                u = ((1.0-ky)*((1.0-kx)*U[ia,ja]+kx*U[ib,ja])+
                (ky)*((1.0-kx)*U[ia,jb]+kx*U[ib,jb]))
                v = ((1.0-ky)*((1.0-kx)*V[ia,ja]+kx*V[ib,ja])+
                (ky)*((1.0-kx)*V[ia,jb]+kx*V[ib,jb]))
                uv_mag = ((1.0-ky)*((1.0-kx)*UV_mag[ia,ja]+kx*UV_mag[ib,ja])+
                (ky)*((1.0-kx)*UV_mag[ia,jb]+kx*UV_mag[ib,jb]))
                dt = 2.0*(frwd-0.5)*dt0/(1.0+uv_mag)

                x0_new = x0 + dt*u
                y0_new = y0 + dt*v
                delta += np.sqrt((x0_new-x0)**2+(y0_new-y0)**2)
                if x0 == x0_new: flag=True
                if y0 == y0_new: flag=True
                x0 = x0_new
                y0 = y0_new
                xx += [x0]
                yy += [y0]

                k+=1
        
            Dx = np.max(xx)-np.min(xx)
            Dy = np.max(yy)-np.min(yy)
            D = np.sqrt((Dx*Dx+Dy*Dy))
            
            if (delta>0.0):
                c = delta/D>c0
                if c:
                    xxA += [*xx, np.nan]
                    yyA += [*yy, np.nan]

                else:
                    xxB += [*xx, np.nan]
                    yyB += [*yy, np.nan]
        
        return xxA,xxB,yyA,yyB

    M = 4000
    init_x=xlim[0]+(xlim[1]-xlim[0])*(np.random.rand((M)))
    init_y=ylim[0]+(ylim[1]-ylim[0])*(np.random.rand((M)))

    xxA,xxB,yyA,yyB = [],[],[],[]
    xxA,xxB,yyA,yyB = run(init_x,init_y,xxA,xxB,yyA,yyB,max_delta=2.0*T,c0=2.0)



    fig,ax = plt.subplots(1,1)
    
    ax.plot(xxB,yyB,'-',solid_capstyle='round', solid_joinstyle='round',lw=lw_nc,c=color_nc)
    ax.plot(xxA,yyA,'-',solid_capstyle='round', solid_joinstyle='round',lw=lw_c,c=color_c)
    dUdx = (U[1:,1:]-U[:-1,1:])/(X[1:,1:]-X[:-1,1:])
    dUdy = (U[1:,1:]-U[1:,:-1])/(Y[1:,1:]-Y[1:,:-1])
    dVdx = (V[1:,1:]-V[:-1,1:])/(X[1:,1:]-X[:-1,1:])
    dVdy = (V[1:,1:]-V[1:,:-1])/(Y[1:,1:]-Y[1:,:-1])
    O = (dVdx-dUdy)
    
    #m=ax.imshow(O.T, cmap=plt.cm.seismic,
    #extent=[xmin, xmax, yB, yT],
    #vmin=-10.0,vmax=10.0,
    #interpolation='bicubic'
    #)



    
    #x_spc = [0.5*(xlim[0]+xLE), xLE, 0.5*(xLE+xTE), xTE, 0.5*(xlim[1]+xTE)]
    #I_spc = [np.argmin(abs(x-x_spc[i])) for i in range(len(x_spc))]

    #for i in I_spc:
    #    y1 = Y[i,:]
    #    x1 = X[i,:]
    #    x2 = X[i,:]+U[i,:]/5.0
        
    #    ax.plot(x1, y1, c='blue', lw=1.0)
    #    ax.plot(x2, y1, c='blue', lw=1.0)
    #    ax.plot([x1,x2],[y1,y1], c='blue', lw=0.8)

    ax.fill(plate_x,plate_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=100)
    ax.fill(plate_x,plate_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=200)
    
    ax.fill(wallT_x,wallT_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=100)
    ax.fill(wallT_x,wallT_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=200)

    ax.fill(wallB_x,wallB_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=100)
    ax.fill(wallB_x,wallB_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=200)

    ax.set(xlim=xlim, ylim=ylim)
    ax.set(aspect='equal')
    
    dxlim = xlim[1]-xlim[0]
    dylim = ylim[1]-ylim[0]
    dydx = dylim/dxlim
    cm = 1.0/2.54
    mm = cm/10.0
    fig_w = 0.8*20*cm
    fig_h = dydx*fig_w

    #cbar_pos = [
    #    (pad+1.0+padc)/total_w,
    #    pad/total_h,
    #    cw/total_w,
    #   dydx/total_h
    #   ]

    
    


    ax.set(xlim=xlim, ylim=ylim)
    ax.set(aspect='equal')


    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes('right', size=0.2*cm, pad=0.2*cm)
    #fig.colorbar(m,cax=cax,orientation='vertical')
    fig.set_size_inches(fig_w,fig_h)

    
    plt.show()