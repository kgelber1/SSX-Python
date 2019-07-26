import numpy as np
import scipy as sp
from scipy import interpolate


__author__ = "Katie Gelber"
__version__ = "1.0"

"""Converts the 97.7/155 ratio to temperature using interpolated PRISMSPECT data.
   --> used by VUV functions
   Really shouldn't need to be edited!"""

global data#I know using global varibales isn't great but it's convient.

#Here's the main array, containing temp, then ratio points for 1e1, 5e14, 2e15, 5e15
data = np.array([[  5.00000000e+00,   1.00000000e+01,   1.50000000e+01,
2.00000000e+01,   2.50000000e+01,   3.00000000e+01,
3.50000000e+01,   4.00000000e+01,   4.50000000e+01,
5.00000000e+01,   5.50000000e+01,   6.00000000e+01,
6.50000000e+01,   7.00000000e+01,   7.50000000e+01,
8.00000000e+01,   8.50000000e+01,   9.00000000e+01,
9.50000000e+01,   1.00000000e+02],#Te
[  6.38344000e+01,   1.55291800e-01,   3.82801900e-02,
1.61070800e-02,   8.13544300e-03,   4.95377400e-03,
3.38400000e-03,   2.49453300e-03,   1.93987800e-03,
1.56884400e-03,   1.30707100e-03,   1.11445400e-03,
9.63190400e-04,   8.44969200e-04,   7.50490700e-04,
6.73543000e-04,   6.09871400e-04,   5.56443100e-04,
5.11066200e-04,   4.72117900e-04],#1e14
[  6.92516400e+01,   9.09328500e-02,   2.60355900e-02,
1.15609300e-02,   6.04283000e-03,   3.76594000e-03,
2.61975700e-03,   1.95819600e-03,   1.53970700e-03,
1.25666000e-03,   1.05490800e-03,   9.05227600e-04,
7.86833200e-04,   6.93724100e-04,   6.18896500e-04,
5.57642500e-04,   5.06722600e-04,   4.63815700e-04,
4.27235900e-04,   3.95728400e-04],#5e14
#Temps at densitites of 2e15
[  9.97906400e+00,   5.10684200e-02,   1.61115100e-02,
7.55031000e-03,   4.07848000e-03,   2.60392200e-03,
1.84527300e-03,   1.39922700e-03,   1.11413000e-03,
9.19338300e-04,   7.78921200e-04,   6.73739400e-04,
5.89920000e-04,   5.23148100e-04,   4.69174800e-04,
4.24882100e-04,   3.87730700e-04,   3.56403900e-04,
3.29381100e-04,   3.06235200e-04],#2e15
#temps at density of 5e15
[6.30777004e+00, 2.47194650e-02, 9.20617892e-03, 4.95206629e-03,
 2.87000382e-03 ,1.91298756e-03 ,1.39789432e-03 ,1.08362991e-03,
 8.77756991e-04, 7.34375597e-04 ,6.29147911e-04 ,5.49204923e-04,
 4.84760452e-04 ,4.32782756e-04, 3.90412181e-04, 3.55418095e-04,
 3.25820761e-04 ,3.00750025e-04 ,2.78944776e-04, 2.60241381e-04]])#5e15


def vuvLookupCurve(ratio, den):
    '''Given the ratio at a value, and the density INDEX of the plasma, this will interpolate the PRISMSPECT
    data points using a scipy spline and return the electron temperature.
    Params:
            ratio - int/float of the ratio of 97.7/155 at one point (don't think it can be an array)
            den - int that is the INDEX of the density. VUV function should do the conversion, otherwise
                    see above

    Return:
            temp - a real number corresponding to the elctron temperature.
                    Normal result is around 5-20 ish'''
    global data
    x = data[0]
    xr = x[::-1]
    y = data[1]
    z = data[2]
    w = data[3]
    v = data[4]
    #Use python's spline interpolation. I tried a few methods, log, quadratic decay curve fits
    #but this worked the best.

    #sets up a different range of interpolation for the two highest densities
    wlogratio = np.arange(-8.10, 2.2, .01)
    wrat = np.exp(wlogratio)
    tcks = sp.interpolate.splrep(np.log(w[::-1]), xr, k = 1, s = 0)
    w4 = sp.interpolate.splev(wlogratio, tcks)

    vlogratio = np.arange(-8.25, 1.8, .01)
    vrat = np.exp(vlogratio)
    tck = sp.interpolate.splrep(np.log(v[::-1]), xr, k = 1, s = 0)
    v4 = sp.interpolate.splev(vlogratio, tck)

    if(den == 2):
        temp = sp.interpolate.splev(np.log(ratio), tcks)
    elif(den == 3):
        temp = sp.interpolate.splev(np.log(ratio), tck)
    else:
        te = np.arange(5,100,.001)
        logratio = np.arange(-8.13, 0, .1)
        rat = np.exp(logratio)
        data = [y,z,w,v]
        datar = [dat[::-1] for dat in data]
        splr = [sp.interpolate.splrep(np.log(dat), xr, k=3) for dat in datar]
        splines = [sp.interpolate.splev(logratio, spl) for spl in splr]
        y4, z4, w4, v4= splines
        splw = splr[den]
        temp = sp.interpolate.splev(np.log(ratio), splw)


    return temp
