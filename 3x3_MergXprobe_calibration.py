# -*- coding: utf-8 -*-
"""
Created on Wed Jun 06 10:27:32 2019

@authors: Manjit, Katie, Adam


Generates calibration matrix for the 4x4 probe array.
Outputs a 48 x 4 array, the first three columns corresponding to the Bx, By, Bz.
These values correct for the sign (some loops are inverted) and "provide
information about how much the B field is being picked up by the rest two components"
The fourth column is a normalization coefficent of sqrt(Bx^2 + By^2 + Bz^2)

Calibration data was taken with several shots of the probe stalks in the x,y,z
direction.
"""

import os
from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pylab as plt

# and from the ssx tools:
import mjtools as mj
import ssx_data_read_basic as sdr
import mjtools as mj
import basic_hiresmag as hdr

def loadcurrent(shot, scope = '3'):
    data = sdr.scope_data(shot, scope)
    time = data.time
    eastcurrent = -data.ch2*170000#East gun 170kA/V
    westcurrent = -data.ch4*190000#West gun 190kA/V
    return time,eastcurrent,westcurrent

def load_smallbdot(shot, scope = '1'):
    data = sdr.scope_data(shot, scope)
    time = data.time
    X_bdot = data.ch1 #*170000#East gun 170kA/V
    Y_bdot = data.ch2 #*190000#West gun 190kA/V
    Z_bdot = data.ch3
    Bx = sp.integrate.cumtrapz(X_bdot,time)
    By = sp.integrate.cumtrapz(Y_bdot,time)
    Bz = sp.integrate.cumtrapz(Z_bdot,time)
    timeb = time[:-1]
    return timeb, Bx, By, Bz

def polyPeak_noPlot(time, data, timerange = [40,80],axis = 'x'):
    """Finds the peak in the specified time range.
    Finds the peak in the data in the specified time range.  You must pass it
    pretrigger and timebase as well as the time and data series.
    - timebase - is an integer equal to (1 us)/(delta t of the digitizer).  For
      example, if capturing data at 10 MHz, timebase = (1e-6/1e-7) = 10.  For
      data captured at 2 MHz, timebase = (1e-6/5e-7) = 2.  In other words, how
      many data points per us.  It's a pretty bad way to code this up, but it
      lets you specify the time range in micro seconds in which to look for the
      peak.
    - pretrigger - how many microseconds of data before t = 0."""

    # Find the indices corresponding to the ends of the time range
    t1 = mj.tindex(time,timerange[0])#+pretrigger)
    t2 = mj.tindex(time,timerange[1])#+pretrigger)
    # print 't1=', t1
    # print 't2=', t2
    # generate an array of indices spanning the range
    ti = np.arange(t1,tn 2)
    # get the time and data points in the range
    t = time[ti]
    d = data[ti]
    # Fit a 2nd degree polynomial and find the min and max.
    p = np.polyfit(t,d,2)
    fit = p[0]*t**2 + p[1]*t + p[2]
    dataMax = fit.max()
    dataMin = fit.min()

    if abs(dataMin) > dataMax:
        dataMax = dataMin
    return dataMax


#compute helmholz coil Bfield
def helmholtz2(r, i = 1.0, coil = 2):
    """Compute B field of Helmholtz coil at (x,y,z)

    i is current in amps and coil is the coil selection.
        - 1 for the old wooden coil
        - 2 for the new delrin coil"""
    r1, r2, r3 = r

    # a is the radius of the coil in meters, and d is the distance between the
    # two coils.
    if coil == 1:
        # coil 1 is the old wooden coil.  It has two turns per side, a radius
        # of 6.1" and a separation of 5.9".
        a = 6.1 * 0.0254
        d = 5.9 * 0.0254
        turns = 2
    elif coil == 2:
        # Coil 2 is the new delrin and copper pipe tube built in 2011.  Each
        # side has one turn.  The ID of the pipe is 6.125" and the OD is 6.625,
        # so we split the difference to calculate the radius.  Same goes for
        # the separation - the coils are 3/4" high, so we go from the midline.
        # These numbers were come from page 26-28 in Tim's 2011 lab notebook
        # and helmholtz.nb.
        a = 0.0809625
        d = 0.085725
        turns = 1

    b = np.zeros(3)

    for j in xrange(360):	# compute line integral
        dth = 1./360 * 2 * np.pi
        th = j * dth

        rho = np.sqrt( (r1 - a * np.cos(th))**2 + (r2 - a * np.sin(th))**2 + (r3 - .5 *
            d)**2 )
        b[0] = b[0] + a * dth * (- r3 * np.cos(th))/rho**3
        b[1] = b[1] + a * dth * (- r3 * np.sin(th))/rho**3
        b[2] = b[2] + a * dth * ( (r1 - a * np.cos(th)) * np.cos(th) + (r2 - a *
            np.sin(th)) * np.sin(th) )/rho**3

        rho = np.sqrt( (r1 - a * np.cos(th))**2 + (r2 - a * np.sin(th))**2 + (r3 + .5 *
            d)**2 )
        b[0] = b[0] + a * dth * (- r3 * np.cos(th))/rho**3
        b[1] = b[1] + a * dth * (- r3 * np.sin(th))/rho**3
        b[2] = b[2] + a * dth * ( (r1 - a * np.cos(th)) * np.cos(th) + (r2 - a *
            np.sin(th)) * np.sin(th) )/rho**3


    # calculate the B field.  The 4pi's cancel out, and the 1e4 is to convert
    # to gauss
    b = b * i * 1.0e-7 * turns * 1e4
    return b

def probe_calib(shot, probe_num, position_vector, dir):
    magdata = hdr.getMagData(shot)

    uncalibB_x=sp.signal.detrend(magdata.iUnCalibData[0,probe_num,:])
    uncalibB_y=sp.signal.detrend(magdata.iUnCalibData[1,probe_num,:])
    uncalibB_z=sp.signal.detrend(magdata.iUnCalibData[2,probe_num,:])

    magmax_x = polyPeak_noPlot(magdata.time,uncalibB_x,[10,100], axis = 'x')
    magmax_y = polyPeak_noPlot(magdata.time,uncalibB_y,[10,100], axis = 'y')
    magmax_z = polyPeak_noPlot(magdata.time,uncalibB_z,[10,100], axis = 'z')

    guntime,east,west=loadcurrent(shot)
    currmax = polyPeak_noPlot(guntime,east,[10,100],axis = 'current')

    helmB=helmholtz2(position_vector,currmax)
    if dir == 2 :#r - z direction of probe is along axis of Helmholtz coil
        # r shots : x,y,z - > r,t,z
        helmB = [ 1*helmB[2], 1*helmB[0], 1*helmB[1] ]
    elif dir == 0:#t-  x direction of probe is along axis of Helmholtz coil
    #### NEED TO VERIFY WHETHER THIS IS CORRECT
        # t shots : x,y,z - > r,t,z
        helmB = [ 1*helmB[0], 1*helmB[2], 1*helmB[1] ]
    elif dir ==1:#z - y direction of probe is along axis of Helmholtz coil
        # same as t, but different transform below
        # ADL 20190501 (coils are in same positions, but field directions rotate
        # z shots : x,y,z - > r,t,z
        helmB = [ 1*helmB[0], 1*helmB[1], 1*helmB[2] ]

    ratio_x = magmax_x/helmB[0]
    ratio_y = magmax_y/helmB[1]
    ratio_z = magmax_z/helmB[2]
    # print 'Mag Max = ', magmax
    # print 'HelmB = ', helmB
    # print 'Ratio = ', ratio
    # return ratio_x,ratio_y, ratio_z, currmax,helmB

    return ratio_x, ratio_y, ratio_z, currmax, helmB

def get_probeLocs(dir, num_probes = 16):
    """
    Provides positions in meters along probe stalks for 4x4 array of probes built by M. Kaur.
    ADL -- 20190501:
    KG  -- 20190604
    -----
    This code generates the 16 by 3 array of positions. Each of the 16 rows
    repersents the position of a single probe locations, and the three
    elements are the x, y and z directions. Also worth noting that this uses
    a LEFT handed -cordinate system. This function returns porbe locations for
    probes in a square configuration
    """
    position_vectors = [[0] * 3 for i in range(num_probes)]

    #every x postion
    # x_pos = [-4.25*1e-3, -1.42*1e-3, 1.41*1e-3, 4.24*1e-3]
    x_pos = [-4.25*1e-3, -4.25*1e-3, 4.24*1e-3, 4.24*1e-3]
    y_pos = [-4.25*1e-3,4.24*1e-3, 4.24*1e-3, -4.25*1e-3]
    z_pos = [-2.25*1e-3*25.4, -0.75*1e-3*25.4, 0.75*1e-3*25.4, 2.25*1e-3*25.4]
    x = 0
    for i in range(num_probes):
        if(i%4 ==0 and i>0):
            x+=1
        position_vectors[i][0] = x_pos[x]
        position_vectors[i][1] = y_pos[x]
        position_vectors[i][2] = z_pos[i%4]
        # print(position_vectors[i][0])

    """ Now  take into account the direction
        r shots : x,y,z - > r,t,z
        t shots : x,y,z - > r,t,z
        z shots : x,y,z - > r,t,z
     """
    if dir ==2 :#r
         # don't need to switch anything
        return position_vectors
    if dir == 0:#t
         # looks like -90 degree rotation about y-axis of probe coil orientation, so switch x and z
        position_vectors[:][0], position_vectors[:][2] = position_vectors[:][2], position_vectors[:][0]
        return position_vectors
    if dir ==1:#z
         # also like -90 degree rotation, switch x and z
         position_vectors[:][0], position_vectors[:][2] = position_vectors[:][2], position_vectors[:][0]
         return position_vectors

    return position_vectors

def getRatio(probe_num, position_vector, shot_range, dir, day ='050119r'):
    """Written by M. Kaur
    --------------------------------
    KG  -- 20190605

    This (I think) goes though every shot, and finds the maximum magetic field
    then averages the maximum signal over several shots
    """
    ratio_x = 0
    ratio_y = 0
    ratio_z = 0
    # helm_B = 0
    divideby  = 0
    for shot in range(shot_range[0], shot_range[1]+1):
        print 'On shot ', day+str(shot), ' for probe ',probe_num
        x,y,z, currmax,helmB = probe_calib(day+str(shot), probe_num, position_vector,dir)
        ratio_x = ratio_x + x
        ratio_y = ratio_y + y
        ratio_z = ratio_z + z
        # helm_B = helm_B + helmB
        divideby = divideby + 1 #averaging over the number of shots
        ratio_Bx = ratio_x/divideby
        ratio_By = ratio_y/divideby
        ratio_Bz = ratio_z/divideby
        # helmB = helm_B/divideby

    # Bx_sqr =[ratio_x[i]**2 for i in range(len(ratio_x))]
    Bx_sqr =ratio_x**2
    By_sqr =ratio_y**2
    Bz_sqr =ratio_z**2
    B = Bx_sqr + By_sqr+ Bz_sqr
    normfactor = np.sqrt(B)
    ratio_Bx, ratio_By, ratio_Bz = [ratio_Bx, ratio_By, ratio_Bz]/normfactor
        # print ratio_Bx, ratio_By, ratio_Bz, helmB
    # print("ratio_Bx %f, ratio_By %f, ratio_Bz %f, helmB%s"%(ratio_Bx, ratio_By, ratio_Bz, helmB))

    #going to need to normalize all the ratios before returning them
    #also going to need to find the magnidue = mu0 I/(2*pi*r)
    return (ratio_Bx, ratio_By, ratio_Bz)

def main():
    """This function is where any user-specified values shoulbe be.
    Given shots in the x, y and z direction, it
    finds position vectors (a 16 by 3 array becuase there are 16 probe locations
    (4 locations on 4 probes) and each probe location has an x, y, and z
    coordinate relative to the center of the Bfield) then uses those positions to
    calculate the magnetic field at that point. Then, """
    num_probes = 16

    # print(position_vectors)
    # Create the (48x4) calibration matrix:
    # calibration_matrix = [[0] * 4 for i in range(num_probes*3)]
    calibration_matrix = [[0] * 9 for i in range(num_probes)]
    counter = 0

    # first populate with x-direction:
    shot_range = [17, 20] #x-direction
    dir = 0 #the direction of the orentation of the probe array
    position_vectors = get_probeLocs(dir)
    for probe_num in range(num_probes):
        position_vector = position_vectors[probe_num]
        ratio_Bx, ratio_By, ratio_Bz = getRatio(probe_num, position_vector, shot_range, dir)
        print("Progress: %d / %d" %(counter+1,num_probes*3 ))
        calibration_matrix[probe_num][(dir*3):(dir*3+1)] = [ratio_Bx, ratio_By, ratio_Bz]
        counter +=1

    # Then populate with y-direction:
    shot_range = [21, 25] #y-direction
    dir = 1 #the direction of the orentation of the probe array
    position_vectors = get_probeLocs(dir)
    for probe_num in range(num_probes):
        position_vector = position_vectors[probe_num]
        ratio_Bx, ratio_By, ratio_Bz = getRatio(probe_num, position_vector, shot_range, dir)
        print("Progress: %d / %d" %(counter+1,num_probes*3 ))
        calibration_matrix[probe_num][(dir*3):(dir*3+1)] = [ratio_Bx, ratio_By, ratio_Bz]
        counter +=1

    # Then populate with z-direction:
    shot_range = [11, 15] #z-direction
    dir = 2
    position_vectors = get_probeLocs(dir)
    for probe_num in range(num_probes):
        position_vector = position_vectors[probe_num]
        ratio_Bx, ratio_By, ratio_Bz = getRatio(probe_num, position_vector, shot_range, dir)
        print("Progress: %d / %d" %(counter+1,num_probes*3 ))
        calibration_matrix[probe_num][(dir*3):(dir*3+1)] = [ratio_Bx, ratio_By, ratio_Bz]
        counter +=1

    pth = os.getcwd()
    date = '050119'
    print("Finished! File saved as calib-%s-4x4_diff.txt in current working directory" %(date))
    savetxt(os.path.join(pth, 'calib-%s-4x4_diff.txt' % (date)) , calibration_matrix)


if __name__ == "__main__":
    main()
