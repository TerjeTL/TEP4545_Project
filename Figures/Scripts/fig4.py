import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def d(a):
    return a[1:]-a[:-1]

fig_path = 'Figures/Figures/'

cases = [
    'WD1.0LT5.0',
    'WD1.0LT6.0',
#   'WD1.0LT7.0',
    'WD1.0LT8.0',
#   'WD2.0LT3.0',
#   'WD2.0LT4.0',
    'WD2.0LT5.0',
#   'WD2.0LT6.0',
#   'WD2.0LT7.0',
    'WD2.0LT8.0',
#   'WD3.0LT3.0',
#   'WD3.0LT4.0',
    'WD3.0LT5.0',
    'WD3.0LT6.0',
#   'WD3.0LT7.0',
    'WD3.0LT8.0',
]

GENERATE = False
RUN = False

def generate(case):
    filename = 'Data/{case}/domain_data.csv'.format(case=case)
    T = 0.01
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
    print(case)
    for j in range(J):
        print('.',end='')
        for i in range(I):
            x0,y0 = X[i,j], Y[i,j]
            dx,dy = x0-data[:,0], y0-data[:,1]
            r2 = dx*dx+dy*dy
            idx = np.argmin(r2)
            row = data[idx,:]
            U[i,j],V[i,j],P[i,j] = row[3],row[4],row[12]
    print()
    np.save('Data/{case}/X'.format(case=case), X)
    np.save('Data/{case}/Y'.format(case=case), Y)
    np.save('Data/{case}/U'.format(case=case), U)
    np.save('Data/{case}/V'.format(case=case), V)
    np.save('Data/{case}/P'.format(case=case), P)

def run(case):
    
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

    x = X[:,0]
    y = Y[0,:]

    i = np.argmin(abs(x-xTE))
    jt = np.argmin(abs(y-yt))
    jT = np.argmin(abs(y-yT))+1
    u = U[i,jt:jT]
    u = u/np.max(u)
    y = (Y[i,jt:jT]-yt)/(yT-yt)
    plt.plot(u,y, '-x')

def streamplot(x,y,U,V,UV_mag,init_x,init_y,dt0=3e-2,max_delta=1.0,max_k=10000,c0=1.57):
    M = len(init_x)
    xxA,xxB,yyA,yyB = [],[],[],[]
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
            dt = 2.0*(frwd-0.5)*dt0/(1.0+uv_mag*uv_mag)

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

def plot_streamplot(case,figname='stream'):
    
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

    xlim = [xLE-0.5*T, xLE+0.6*L]
    ylim = [yb-0.6*T, 0.5*(yt+yb)]

    color_c=(1.0, 0.0, 0.0)
    color_nc=(0.8, 0.8, 0.8)
    lw_c = 0.6
    lw_nc = 0.4

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
    plate_h_lw =  0.6
    plate_lw = 1.0

    x = X[:,0]
    y = Y[0,:]

    #plt.imshow(UV_mag.T, cmap=plt.cm.jet,extent=[xmin, xmax, yT, yB],interpolation='bicubic')
    #plt.contour(X.T,Y.T,UV_mag.T,colors=['black'], levels=[0.0])
    #plt.quiver(X,Y,U/UV_mag,V/UV_mag,scale=100.0)
    #plt.show()



    M = 1000
    init_x=xlim[0]+(xlim[1]-xlim[0])*(np.random.rand((M)))
    init_y=ylim[0]+(ylim[1]-ylim[0])*(np.random.rand((M)))

    
    xxA,xxB,yyA,yyB = streamplot(x,y,U,V,UV_mag,init_x,init_y,max_delta=4.0*T,c0=np.pi/np.sqrt(2))

    fig,ax = plt.subplots(1,1)
    
    ax.plot(xxB,yyB,'-',solid_capstyle='round', solid_joinstyle='round',lw=lw_nc,c=color_nc)
    ax.plot(xxA,yyA,'-',solid_capstyle='round', solid_joinstyle='round',lw=lw_c,c=color_c)

    ax.fill(plate_x,plate_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=100)
    ax.fill(plate_x,plate_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=200)
    
    ax.fill(wallT_x,wallT_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=100)
    ax.fill(wallT_x,wallT_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=200)

    ax.fill(wallB_x,wallB_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=100)
    ax.fill(wallB_x,wallB_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=200)

    dxlim = xlim[1]-xlim[0]
    dylim = ylim[1]-ylim[0]
    dydx = dylim/dxlim
    cm = 1.0/2.54
    mm = cm/10.0
    fig_w = 20*cm
    fig_h = dydx*fig_w

    fig.set_size_inches(fig_w,fig_h)
    ax.set(xlim=xlim, ylim=ylim,position=[0.0,0.0,1.0,1.0])
    ax.set_axis_off()
    

    fig.savefig(fig_path+'{figname}.svg'.format(figname=figname),format='svg',dpi=1200)

def plot_streamplot2(case,figname='stream'):
    
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

    xlim = [xLE-0.5*T, xLE+0.6*L]
    ylim = [yb-0.6*T, 0.5*(yt+yb)]

    color_c=(1.0, 0.0, 0.0)
    color_nc=(0.8, 0.8, 0.8)
    lw_c = 0.6
    lw_nc = 0.4

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
    plate_h_lw =  0.6
    plate_lw = 1.0

    x = X[:,0]
    y = Y[0,:]

    #plt.imshow(UV_mag.T, cmap=plt.cm.jet,extent=[xmin, xmax, yT, yB],interpolation='bicubic')
    #plt.contour(X.T,Y.T,UV_mag.T,colors=['black'], levels=[0.0])
    #plt.quiver(X,Y,U/UV_mag,V/UV_mag,scale=100.0)
    #plt.show()



    M = 800
    init_x=xlim[0]+(xlim[1]-xlim[0])*(np.random.rand((M)))
    init_y=ylim[0]+(ylim[1]-ylim[0])*(np.random.rand((M)))

    
    xxA,xxB,yyA,yyB = streamplot(x,y,U,V,UV_mag,init_x,init_y,max_delta=4.0*T,c0=np.pi/np.sqrt(2))

    fig,ax = plt.subplots(1,1)
    
    ax.plot(xxB,yyB,'-',solid_capstyle='round', solid_joinstyle='round',lw=lw_nc,c=color_nc,zorder=-4)
    ax.plot(xxA,yyA,'-',solid_capstyle='round', solid_joinstyle='round',lw=lw_c,c=color_c,zorder=-4)

    ax.fill(plate_x,plate_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
    ax.fill(plate_x,plate_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)
    
    ax.fill(wallT_x,wallT_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
    ax.fill(wallT_x,wallT_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)

    ax.fill(wallB_x,wallB_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
    ax.fill(wallB_x,wallB_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)

    dxlim = xlim[1]-xlim[0]
    dylim = ylim[1]-ylim[0]
    dydx = dylim/dxlim
    cm = 1.0/2.54
    mm = cm/10.0
    fig_w = 20*cm
    fig_h = 1.2*dydx*fig_w

    fig.set_size_inches(fig_w,fig_h)
    
    ax.set(xlim=xlim, ylim=ylim, aspect='equal')
    ax.set_yticks([-0.0,-0.5,-1.0])
    ax.set_xticks([10,11,12,13])
    ax.set_xlabel(r'$x/T$')
    ax.set_ylabel(r'$y/T$')
    plt.tight_layout()
    plt.grid(c='black')
    #ax.set_axis_off()
    

    fig.savefig(fig_path+'{figname}.svg'.format(figname=figname),format='svg',dpi=1200)

def plot_streamplot3(case,figname='stream2'):
    
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

    xa,xb = [xLE-0.5*T, xLE+0.6*L]
    ya,ybb = [yb-0.6*T, 0.5*(yt+yb)]
    xlim = [xLE-T, xLE+2.0*L]
    ylim = [yB-0.5*T, yT+0.5*T]

    color_c=(1.0, 0.0, 0.0)
    color_nc=(0.8, 0.8, 0.8)
    lw_c = 0.6
    lw_nc = 0.4

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
    plate_h_lw =  0.6
    plate_lw = 1.0

    x = X[:,0]
    y = Y[0,:]

    #plt.imshow(UV_mag.T, cmap=plt.cm.jet,extent=[xmin, xmax, yT, yB],interpolation='bicubic')
    #plt.contour(X.T,Y.T,UV_mag.T,colors=['black'], levels=[0.0])
    #plt.quiver(X,Y,U/UV_mag,V/UV_mag,scale=100.0)
    #plt.show()



    M = 2500
    init_x=xlim[0]+(xlim[1]-xlim[0])*(np.random.rand((M)))
    init_y=ylim[0]+(ylim[1]-ylim[0])*(np.random.rand((M)))

    
    xxA,xxB,yyA,yyB = streamplot(x,y,U,V,UV_mag,init_x,init_y,max_delta=4.0*T,c0=np.pi/np.sqrt(2))

    fig,ax = plt.subplots(1,1)
    
    ax.plot(xxB,yyB,'-',solid_capstyle='round', solid_joinstyle='round',lw=lw_nc,c=color_nc,zorder=-4)
    ax.plot(xxA,yyA,'-',solid_capstyle='round', solid_joinstyle='round',lw=lw_c,c=color_c,zorder=-4)

    ax.fill(plate_x,plate_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
    ax.fill(plate_x,plate_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)
    
    ax.fill(wallT_x,wallT_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
    ax.fill(wallT_x,wallT_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)

    ax.fill(wallB_x,wallB_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
    ax.fill(wallB_x,wallB_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)

    ax.fill([xa,xb,xb,xa],[ya,ya,ybb,ybb],fill=False,ec='blue',lw=2*plate_lw,zorder=-1)

    dxlim = xlim[1]-xlim[0]
    dylim = ylim[1]-ylim[0]
    dydx = dylim/dxlim
    cm = 1.0/2.54
    mm = cm/10.0
    fig_w = 20*cm
    fig_h = 1.2*dydx*fig_w

    fig.set_size_inches(fig_w,fig_h)
    
    ax.set(xlim=xlim, ylim=ylim, aspect='equal')
    ax.set_yticks([-2,-1,0,1,2])
    ax.set_xticks([10,11,12,13,14,15,16,17,18,19])
    ax.set_xlabel(r'$x/T$')
    ax.set_ylabel(r'$y/T$')
    plt.tight_layout()
    plt.grid(linestyle=':',c='gray')
    #ax.set_axis_off()
    

    fig.savefig(fig_path+'{figname}.svg'.format(figname=figname),format='svg',dpi=1200)

def plot_quantity(cases,figname='quan'):

    color_c=(1.0, 0.0, 0.0)
    color_nc=(0.8, 0.8, 0.8)
    lw_c = 0.6
    lw_nc = 0.4

    num_cases = len(cases)
    fig,axes = plt.subplots(num_cases,2,sharex=True,sharey=True)

    for case_indx,case in enumerate(cases):
        ax = axes[case_indx,:]
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

        xlim = [xLE-0.5*L, xTE+2.0*L]
        ylim = [yB, yT]

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
        plate_h_lw =  0.6
        plate_lw = 1.0

        x = X[:,0]
        y = Y[0,:]

        dUdx = (U[1:,1:]-U[:-1,1:])/(X[1:,1:]-X[:-1,1:])
        dUdy = (U[1:,1:]-U[1:,:-1])/(Y[1:,1:]-Y[1:,:-1])
        dVdx = (V[1:,1:]-V[:-1,1:])/(X[1:,1:]-X[:-1,1:])
        dVdy = (V[1:,1:]-V[1:,:-1])/(Y[1:,1:]-Y[1:,:-1])
        O = (dVdx-dUdy)

            #ax[1].set(
            #    ylabel='W{}L{}'.format(int(fw),int(fL))
            #)
        #if case_indx==0:

        

        
        for i in range(2):
            ax[i].fill(plate_x,plate_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
            ax[i].fill(plate_x,plate_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)
            
            ax[i].fill(wallT_x,wallT_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
            ax[i].fill(wallT_x,wallT_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)

            ax[i].fill(wallB_x,wallB_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
            ax[i].fill(wallB_x,wallB_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)

        maps = []
        u0 = 1.4776
        o0 = np.abs(np.max(O))
        maps += [ax[0].imshow(UV_mag.T/u0,
                cmap=plt.cm.plasma,
                extent=[xmin, xmax, yT, yB],
                interpolation='bicubic',
                zorder=-3,
                vmin=0.0,
                vmax=2.0
                )]
        maps += [ax[1].imshow(O.T/o0,
                cmap=plt.cm.seismic,
                extent=[xmin, xmax, yT, yB],
                interpolation='bicubic',
                zorder=-3,
                vmin=-0.1,
                vmax=0.1
                )]

        dxlim = xlim[1]-xlim[0]
        dylim = ylim[1]-ylim[0]
        dydx = dylim/dxlim
        cm = 1.0/2.54
        #mm = cm/10.0
        #fig_w = 20*cm
        #fig_h = dydx*fig_w
        #fig.set_size_inches(fig_w,2*fig_h)

        pad = 0.04
        w = (1.0-3.0*pad)/2.0
        wspace = pad/w
        cbar_h = pad/4.0

        for i in range(2):
            ax[i].set(xlim=xlim, ylim=ylim)
            if case_indx==num_cases-1:
                fig.subplots_adjust(
                    bottom=2*pad,
                    left=pad,
                    right=1.0-pad,
                    top=1.0-2*pad,
                    wspace=wspace
                )
                if i==0:
                    cax = fig.add_axes([pad,1.0-pad-cbar_h,w,cbar_h])
                    cax.set_title('Velocity magnitude '+r'$\left|\vec{u}\right|/u_0$')
                    colorbar = fig.colorbar(maps[i],cax=cax,orientation='horizontal',ticks=[0,1,2])
                else:
                    cax = fig.add_axes([2*pad+w,1.0-pad-cbar_h,w,cbar_h])
                    cax.set_title('Vorticity '+r'$\omega/\left|\omega\right|_{max}$')
                    colorbar = fig.colorbar(maps[i],cax=cax,orientation='horizontal',ticks=[-0.1,0.0,0.1])
        

        ax[1].tick_params(axis='y',left=True,right=False,labelleft=True,labelright=False)
        ax[0].tick_params(axis='y',left=False,right=True,labelleft=False)
        ax[0].set_ylabel('W{}L{}'.format(int(fw),int(fL)))
        ax[0].set_yticks([-2,0,2])
        ax[1].set_yticks([-2,0,2])


        ax[1].yaxis.set_label_position('right')
        ax[0].yaxis.label.set_size(12)
        ax[1].set_ylabel(r'$y/T$')
        if case_indx==num_cases-1:
            ax[0].set_xlabel(r'$x/T$')
            ax[1].set_xlabel(r'$x/T$')

        fig_s=25*cm
        fig.set_size_inches(fig_s,4.4/8*fig_s)
        fig.savefig(fig_path+'{figname}.svg'.format(figname=figname),format='svg',dpi=1200)

def vel_plot(case,figname='false'):

    color_c=(1.0, 0.0, 0.0)
    color_nc=(0.8, 0.8, 0.8)
    lw_c = 0.6
    lw_nc = 0.4


    fig,ax = plt.subplots(1,1)

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

    xlim = [xLE-0.35*L, xTE+0.5*L]
    ylim = [yB-0.5*T, yT+0.5*T]

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
    plate_h_lw =  0.6
    plate_lw = 1.0



    dUdx = (U[1:,1:]-U[:-1,1:])/(X[1:,1:]-X[:-1,1:])
    dUdy = (U[1:,1:]-U[1:,:-1])/(Y[1:,1:]-Y[1:,:-1])
    dVdx = (V[1:,1:]-V[:-1,1:])/(X[1:,1:]-X[:-1,1:])
    dVdy = (V[1:,1:]-V[1:,:-1])/(Y[1:,1:]-Y[1:,:-1])

    ax.fill(plate_x,plate_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
    ax.fill(plate_x,plate_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)
    
    ax.fill(wallT_x,wallT_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
    ax.fill(wallT_x,wallT_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)

    ax.fill(wallB_x,wallB_y,fc=plate_fc,ec=plate_hc,lw=plate_h_lw,hatch='///',zorder=-2)
    ax.fill(wallB_x,wallB_y,fill=False,ec=plate_ec,lw=plate_lw,zorder=-1)


    u0 = 1.4776
    U = U
    V = V
    x = X[:,0]
    y = Y[0,:]

    dxlim = xlim[1]-xlim[0]
    dylim = ylim[1]-ylim[0]
    dydx = dylim/dxlim
    cm = 1.0/2.54

    #mm = cm/10.0
    #fig_w = 20*cm
    #fig_h = dydx*fig_w
    #fig.set_size_inches(fig_w,2*fig_h)
    jt = np.argmin(abs(yt-y))
    jT = np.argmin(abs(yT-y))+1
    jb = np.argmin(abs(yb-y))+1
    jB = np.argmin(abs(yB-y))
    iLE = np.argmin(abs(xLE-x))
    iTE = np.argmin(abs(xTE-x))+1
    di = (iTE-iLE)//4
    indx = [iLE-di, iLE, iLE+di, iLE+2*di, iLE+3*di, iTE-1, iTE+di]
    for i in indx:
        if i >= iLE and i < iTE:
            x=np.array([*X[i,jB:jb], np.nan, *X[i,jt:jT]])
            y=np.array([*Y[i,jB:jb], np.nan, *Y[i,jt:jT]])
            u=np.array([*U[i,jB:jb], np.nan, *U[i,jt:jT]])
            v=np.array([*V[i,jB:jb], np.nan, *V[i,jt:jT]])
        else:
            x,y,u,v=X[i,:],Y[i,:],U[i,:],V[i,:]
        ax.plot(x,y,c='blue')
        S = 0.3
        ax.plot([x[::4],(x+S*u)[::4]],[y[::4],(y)[::4]],c='blue',lw=0.5)
        ax.plot(x+S*u,y,c='blue')
        if i >= iLE and i < iTE:
            ax.plot([x[0],x[0]+S*u0,x[0]+S*u0,x[0],x[0]],[yB,yB,yb,yb,yB],'--',c='red')
            ax.plot([x[0],x[0]+S*u0,x[0]+S*u0,x[0],x[0]],[yt,yt,yT,yT,yt],'--',c='red')
        else:
            ax.plot([x[0],x[0]+S*u0,x[0]+S*u0,x[0],x[0]],[yB,yB,yT,yT,yB],'--',c='red')
            
    
        #ax.plot([x[0],x[0]],[-0.1,0.1],c='red')
        #ax.plot([x[0]+u0/3,x[0]+u0/3],[-0.1,0.1],c='red')
    
    ax.set(xlim=xlim,ylim=ylim)
    ax.set_aspect('equal')
    x = X[:,0]
    y = Y[0,:]
    xticks = x[indx]
    xticklabels = [
        r'$x_{LE}-\frac{1}{4}L$', r'$x_{LE}$', r'$x_{LE}+\frac{1}{4}L$', r'$x_{LE}+\frac{2}{4}L$', r'$x_{LE}+\frac{3}{4}L$', r'$x_{TE}$', r'$x_{TE}+\frac{1}{4}L$'
    ]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_yticks([])
    plt.grid()
    plt.tight_layout()
    fig.set_size_inches(13*cm,10*cm)


    fig.savefig(fig_path+'{figname}.svg'.format(figname=figname),format='svg',dpi=1200)

#for case in cases[0:5]:
#    run(case)
#plt.show()
#case='WD1.0LT5.0'
'''
plot_quantity([
    'WD1.0LT5.0',
    'WD1.0LT8.0',
    'WD2.0LT5.0',
    'WD2.0LT8.0'
],figname='quan')
plt.show()
'''
case = 'WD2.0LT5.0'
plot_streamplot3(case)
#vel_plot(case)
plt.show()