import numpy as np
import matplotlib.pyplot as plt


def d(a):
    return a[1:]-a[:-1]

cases = [
    'WD1.0LT3.0',
    'WD1.0LT4.0',
    'WD1.0LT5.0',
    'WD1.0LT6.0',
    'WD1.0LT7.0',
    'WD1.0LT8.0',
    'WD2.0LT3.0',
    'WD2.0LT4.0',
    'WD2.0LT5.0',
    'WD2.0LT6.0',
    'WD2.0LT7.0',
    'WD2.0LT8.0',
    'WD3.0LT3.0',
    'WD3.0LT4.0',
    'WD3.0LT5.0',
    'WD3.0LT6.0',
    'WD3.0LT7.0',
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
    for j in range(J):
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

for case in cases:
    run(case)
plt.show()

