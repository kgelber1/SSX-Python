import numpy as np
import scipy as sp
from scipy import interpolate
from scipy import optimize
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import CheckCurve as cc

'''Use this to:
                1. interpolate ratio vs temp data at a differnt density or
                2. Make a plot of the ratio vs temp data at various densitites'''

global data
data = np.array([[  5.00000000e+00,   1.00000000e+01,   1.50000000e+01,
2.00000000e+01,   2.50000000e+01,   3.00000000e+01,
3.50000000e+01,   4.00000000e+01,   4.50000000e+01,
5.00000000e+01,   5.50000000e+01,   6.00000000e+01,
6.50000000e+01,   7.00000000e+01,   7.50000000e+01,
8.00000000e+01,   8.50000000e+01,   9.00000000e+01,
9.50000000e+01,   1.00000000e+02],
[  6.38344000e+01,   1.55291800e-01,   3.82801900e-02,
1.61070800e-02,   8.13544300e-03,   4.95377400e-03,
3.38400000e-03,   2.49453300e-03,   1.93987800e-03,
1.56884400e-03,   1.30707100e-03,   1.11445400e-03,
9.63190400e-04,   8.44969200e-04,   7.50490700e-04,
6.73543000e-04,   6.09871400e-04,   5.56443100e-04,
5.11066200e-04,   4.72117900e-04],
[  6.92516400e+01,   9.09328500e-02,   2.60355900e-02,
1.15609300e-02,   6.04283000e-03,   3.76594000e-03,
2.61975700e-03,   1.95819600e-03,   1.53970700e-03,
1.25666000e-03,   1.05490800e-03,   9.05227600e-04,
7.86833200e-04,   6.93724100e-04,   6.18896500e-04,
5.57642500e-04,   5.06722600e-04,   4.63815700e-04,
4.27235900e-04,   3.95728400e-04],
#Temps at densitites of 2e15
[  9.97906400e+00,   5.10684200e-02,   1.61115100e-02,
7.55031000e-03,   4.07848000e-03,   2.60392200e-03,
1.84527300e-03,   1.39922700e-03,   1.11413000e-03,
9.19338300e-04,   7.78921200e-04,   6.73739400e-04,
5.89920000e-04,   5.23148100e-04,   4.69174800e-04,
4.24882100e-04,   3.87730700e-04,   3.56403900e-04,
3.29381100e-04,   3.06235200e-04],
#temps at density of 5e15
[3.44171051e+00,  4.20421692e-02,  1.38376247e-02,
6.62850620e-03,3.62604004e-03, 2.33592391e-03,
1.66645748e-03, 1.27006203e-03,1.01573079e-03,
8.41309392e-04, 7.15057077e-04, 6.20156184e-04,
5.44328450e-04, 4.83643604e-04, 4.34491807e-04,
3.94122603e-04,3.60156053e-04, 3.31509340e-04,
3.06697594e-04, 2.85488006e-04]])

def getData():
    #just makes sure the program can see all the data.
    global data
    return data

def quad(x, offset, scale):
    return(scale/((x)**2) + offset)

def quadFit(xs, ys, atPoint):
    popt, pcov = sp.optimize.curve_fit(quad, xs,ys, p0 = (0.03, 1.5e28))
    y = quad(atPoint,*popt)
    return y

def logf(x, offset, scale):
    #helper function for logFit
    return(scale*np.log(x) + offset)

def logFit(xs, ys, atPoint):
    #Runs a logarithmic curve fit to y vs x, and evaluates it
    #at a point or an array of points
    popt, pcov = sp.optimize.curve_fit(logf, xs,ys, p0 = (1, 1.5e28))
    y = logf(atPoint,*popt)
    return y

def findHighDen():
    '''This is the funciton that interpolates data at new densities. This is set
    to interpolate to 5e15. '''
    global data
    x = data[0]
    y = data[1]
    z = data[2]
    w = data[3]
    highDenData = w.copy()

    xs = np.array([1e14,5e14,2e15])
    for i in range(0, len(highDenData)):
        if(i == 1):#Code was breaking, so I omitted the 1e14 point at 5eV
            xs = np.array([5e14,2e15])
            ys = np.array([z[i], w[i]])
            highDenData[i] = logFit(xs, ys, 5e15)

        else:#interpolates at higher density
            xs = np.array([1e14,5e14,2e15])
            ys = np.array([y[i], z[i], w[i]])
            highDenData[i] = logFit(xs, ys, 5e15)

    return highDenData

###############Incomplete - bad curve fit.
def plotAll(highDenData):
    #Plots all the ratio vs time curve.
    global data
    x = data[0]
    xr = x[::-1]
    yy = data[1]
    z = data[2]
    w = data[3]
    v = highDenData
    logratio = np.arange(-8.13, 2.5, .01)
    rat = np.exp(logratio)
    data = [yy,z,w,v]
    datar = [dat[::-1] for dat in data]
    splr = [sp.interpolate.splrep(x =np.log(dat), y= xr,  k=1, s = 0) for dat in datar]
    splines = [sp.interpolate.splev(logratio, spl) for spl in splr]
    y4, z4, w4, v4= splines
    plt.semilogy(y4, rat)
    plt.semilogy(z4, rat)
    plt.semilogy(w4, rat)
    plt.semilogy(v4, rat)
    plt.scatter(x,yy,  label = " $\it{n}$ = 1 x $10^{14}$")
    plt.scatter(x,z,  label = " $\it{n}$ = 5 x $10^{14}$")
    plt.scatter(x,w,  label = " $\it{n}$ = 2 x $10^{15}$")
    plt.scatter(x, v,  label = " $\it{n}$ = 5 x $10^{15}$")

    plt.xlabel('$\\bf{T}$ ($\it{e}$V)', size = 14)
    plt.ylabel('$\\bf{Line}$  $\\bf{Ratio}$ ( I$_{97.7}$ /I$_{155}$ )', size = 14)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    # plt.ylabel('Line Ratio (97.7/155)')
    plt.title("Line Ratio vs T")
    plt.legend()
    # plt.show()
    dest  = r'C:\Users\Katie\Documents\SSX\Pictures\AllCurve.png'
    plt.savefig(dest, format='png', dpi=600)
    plt.show()


def plotEach(highDenData):
    '''Plots the 2e15 and 5e15 ratio vs temp curves
    Uses spline curve - it was the best I could find after trying a series of
    other python and scipy curve fits'''
    global data#gets data
    x = data[0]
    xr = x[::-1]
    y = data[1]
    z = data[2]
    w = data[3]
    v = highDenData

    #
    vlogratio = np.arange(-8.25, 1.8, .01)
    vrat = np.exp(vlogratio)
    tck = sp.interpolate.splrep(np.log(v[::-1]), xr, k = 1, s = 0)
    v4 = sp.interpolate.splev(vlogratio, tck)

    wlogratio = np.arange(-8.10, 2.2, .01)
    wrat = np.exp(wlogratio)
    tcks = sp.interpolate.splrep(np.log(w[::-1]), xr, k = 1, s = 0)
    w4 = sp.interpolate.splev(wlogratio, tcks)


    plt.semilogy(v4, vrat)
    plt.semilogy(w4, wrat)


    plt.scatter(x,v, label = " $\it{n}$ = 5 x $10^{15}$")
    plt.scatter(x, w, label = "$\it{n}$ = 2 x $10^{15}$")



    plt.xlabel('$\\bf{Temperature}$ ($\it{e}$V)', size = 14)
    plt.ylabel('$\\bf{Line}$  $\\bf{Ratio}$ ( I$_{97.7}$ /I$_{155}$ )', size = 14)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    # plt.ylabel('Line Ratio (97.7/155)')
    plt.title("Line Ratio vs Te")
    plt.legend()
    plt.show()
    # dest  = r'C:\Users\Katie\Documents\SSX\Pictures\TempCurve.png'
    # plt.savefig(dest, format='png', dpi=600)



xs = findHighDen()
print('xs', xs)

# interps(xs)#plots all
plotEach( xs)#PLOTS Just 2e15, 5e15
