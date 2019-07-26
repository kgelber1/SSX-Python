# -*- coding: utf-8 -*-
"""
Created on Wed Jun 06 10:27:32 2019
@authors: Katie, Adam
Generates calibration matrix for the 4x4 probe array.
Outputs a 48 x 3 array, the first three columns corresponding to the Bx, By, Bz.
These values correct for the sign (some loops are inverted) and "provide
information about how much the B field is being picked up by the rest two components"


Calibration data was taken with several shots of the probe stalks in the x,y,z
direction.

Updated on 07/10/19
"""

import os
from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pylab as plt
from scipy.integrate import cumtrapz
import numpy.ma as ma

# and from the ssx tools:
import mjtools as mj
import ssx_data_read_basic as sdr
import basic_hiresmag as hdr
import mytools as mt


def get_probeLocs_calib_setup(dir, num_probes = 16):
    """
    Provides positions in METERS along probe stalks for 4x4 array of probes built by M. Kaur.
    ADL -- 20190501:
    -----
    This code generates the 16 by 3 array of positions. Each of the 16 rows
    repersents the position of a single probe locations, and the three
    elements are the x, y and z directions. Also worth noting that this uses
    a LEFT handed -cordinate system. This function returns porbe locations for
    probes in a square configuration

    -- KG 06/04/19

    """
    position_vectors = [[0] * 3 for i in range(num_probes)]

    #every x postion
    # x_pos = [-4.25*1e-3, -1.42*1e-3, 1.41*1e-3, 4.24*1e-3]


    #1e-3 converts to meters
    x_pos = [-4.25*1e-3, -4.25*1e-3, 4.24*1e-3, 4.24*1e-3]
    y_pos = [-4.25*1e-3, 4.24*1e-3, 4.24*1e-3, -4.25*1e-3]
    z_pos = [-2.25*25.4*1e-3, -0.75*25.4*1e-3, 0.75*25.4*1e-3, 2.25*25.4*1e-3]
    x = 0
    for i in range(num_probes):
        if(i%4 ==0 and i>0):
            x+=1
        position_vectors[i][0] =x_pos[x]
        position_vectors[i][1] = y_pos[x]
        position_vectors[i][2] =z_pos[i%4]
    #     # print(position_vectors[i][0])
    #
    # """ Now  take into account the direction
    #     r shots : x,y,z - > r,t,z
    #     t shots : x,y,z - > r,t,z
    #     z shots : x,y,z - > r,t,z
    #  """
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
    ti = np.arange(t1,t2)
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
def helmholtz2(r, i = 1.0, dir= 2, coil = 2):
    """Compute B field of Helmholtz coil at (x,y,z)
    i is current in amps and coil is the coil selection.
        - 1 for the old wooden coil
        - 2 for the new delrin coil

    Please, please enter the radial position in METERS.

    This is super important: All of the calculations are done in SI
    and return units of guass, because it gets converted on the last line

    """
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

    for j in range(360):	# compute line integral
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


    # Note mulitply by 1e4 is to convert to gauss

    # calculate the B field.  The 4pi's cancel out, and the

    b = b * i * 1.0e-7 * turns * 1e4

    if dir == 0:
        bb = [ 1*b[2], 1*b[0], 1*b[1] ]
    elif dir == 1:
        bb = [ 1*b[0], 1*b[2], 1*b[1] ]
    else:
        # z shots : x,y,z - > r,t,z
        bb = [ 1*b[0], 1*b[1], 1*b[2] ]

    return bb



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
    # helm_B = [0,0,0]
    divideby  = 0
    for shot in range(shot_range[0], shot_range[1]+1):
        print( 'On shot ', day+str(shot), ' for probe ',probe_num)
        x,y,z, currmax,helmB_new = probe_calib(day+str(shot), probe_num, position_vector,dir)
        ratio_x = ratio_x + x
        ratio_y = ratio_y + y
        ratio_z = ratio_z + z
        # helm_B = [helm_B[i] + helmB_new[i] for i in len(helmB)]
        divideby = divideby + 1 #averaging over the number of shots
        ratio_Bx = ratio_x/divideby
        ratio_By = ratio_y/divideby
        ratio_Bz = ratio_z/divideby
        # helmB = [helm_B]/divideby
        # print ratio_Bx, ratio_By, ratio_Bz, helmB
    # print("ratio_Bx %f, ratio_By %f, ratio_Bz %f, helmB%s"%(ratio_Bx, ratio_By, ratio_Bz, helmB))
    Bx_sqr =ratio_x**2
    By_sqr =ratio_y**2
    Bz_sqr =ratio_z**2
    B = Bx_sqr + By_sqr+ Bz_sqr
    norm_factor = np.sqrt(B)
    ratio_Bx, ratio_By, ratio_Bz = [ratio_Bx, ratio_By, ratio_Bz]/norm_factor

    return (ratio_Bx, ratio_By, ratio_Bz, norm_factor)



def ratio_4_doc(shot, dir, num_probes = 16):
    """ This finds the ratio between the idealized helmholtz field
        and the actual recoreded signal

        This also corrects for inverted signals... however due
        to what I'm assuming is noise, finding the inverted
        ones are a bit tricky - feel free to uncomment
        the plotting lines and see if it needs adjusments, though
        I did get it working reliably for the 050119 calibration

        (SI units, meters, Tesla)

        -- KG 06/24/19
    """
    # data = [[0] *3 for i in range(num_probes)]
    # magdata = hdr.getMagData(shot)
    probe_locs = get_probeLocs_calib_setup(shot)
    data=hdr.getquikData(shot)
    time,eastcurrent,westcurrent = loadcurrent(shot)#using eastcurrent
    ratios = [[0]*3 for i in range(num_probes)]
    for probe in range(num_probes):
        ratio =1
        inverted = False
        # fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
        B=sp.signal.detrend(cumtrapz(data.unCalibData[dir,probe,:], data.time))
        plot_time = data.time[:-1]
        if(np.max(B[2000:6000]) < abs(np.average(B[2000:6000]))):
            # print("\ninverted!")
            inverted = True
            # B =  B* -1
            # ratio = -1

        r = probe_locs[probe]
        max_current = polyPeak_noPlot(time,eastcurrent)
        # if(np.max(eastcurrent) < -1*(np.min(eastcurrent))):
        #     max_current = -1*np.min(eastcurrent)
        helmB = helmholtz2(r,max_current)

        # THis is intentional! I am only using shots where the cmponent is lined
        # up with the z-direction of the helmholz field
        # helmB[2] = helmB[2]*-1
        max_theoretical = np.max(helmB[2])
        max_measured = polyPeak_noPlot(plot_time, B)


        ratio = ratio * max_theoretical/max_measured
        if ratio > 30000 or ratio < -30000:
            ratio = 0


        ratios[probe][dir] = ratio
        # print("\tRatio is: %f" %(ratio))
        # if(inverted and ratio <0):
        #     print("Inverted and ratio reflects that")
        # elif(not inverted and ratio <0):
        if probe ==1:
            print("\n Ratio: %5f \n\t max_measured: %3f, \n\t max_theoretical: %5f"%(ratio,max_measured,max_theoretical ) )

    # Compute the median of the non-zero elements
    # m = np.median(foo[foo > 0])
    # Assign the median to the zero elements
    # foo[foo == 0] = m
    return ratios

def find_max_time_current(shot):
    time,eastcurrent,westcurrent = loadcurrent(shot)
    t_bounds = [10,65]
    # generate an array of indices spanning the range
    indices1 = [i for i, v in enumerate(time >= t_bounds[0]) if v]
    indices2 = [i for i, v in enumerate(time <= t_bounds[-1]) if v]
    indices = np.intersect1d(indices1, indices2)
    eastcurrent = eastcurrent[indices[0]:indices[-1]]
    # print(len(eastcurrent))
    #find index at which max current occurs:
    max_i = np.argmax(eastcurrent)
    min_i = np.argmin(eastcurrent)
    if(eastcurrent[max_i] < np.absolute(eastcurrent[min_i])):
        return min_i, eastcurrent[min_i], time[min_i]
    else:
        return max_i, eastcurrent[max_i], time[max_i]

def get_b(dir, shot, probe):
    position_vectors = get_probeLocs_calib_setup(dir)
    max_i, max_curr, max_time = find_max_time_current(shot)
    return helmholtz2(position_vectors[probe], max_curr, dir), max_time



def get_v(magdata, probe, time):
    # from itertools import izip as zip, count # izip for maximum efficiency
    # print("here")
    i = mt.find_closest(magdata.time, time, False)
    # i = [i for i, j in enumerate(magdata.time) if j == time]
    # print(i)
    win_size = 25
    Vx = cumtrapz(magdata.unCalibData[0,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
    magdata.unCalibData[0,probe,:], win_size), magdata.time)
    Vy = cumtrapz(magdata.unCalibData[1,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
    magdata.unCalibData[1,probe,:], win_size), magdata.time)
    Vz = cumtrapz(magdata.unCalibData[2,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
    magdata.unCalibData[2,probe,:], win_size), magdata.time)

    return Vx[i], Vy[i], Vz[i]




def generateMatrix():
    """This function is where any user-specified values should be.
    Given shots in the x, y and z direction, it
    finds position vectors (a 16 by 3 array becuase there are 16 probe locations
    (4 locations on 4 probes) and each probe location has an x, y, and z
    coordinate relative to the center of the Bfield) then uses those positions to
    calculate the magnetic field at that point. Then, """
    num_probes = 16

    # print(position_vectors)
    # Create the (48x4) calibration matrix:
    calibration_matrix = [[0] * 9 for i in range(num_probes)]
    # counter = 0

    #Not bothering to average:
    day ='050119'
    x_shot = 17
    y_shot = 23
    z_shot = 12

    x_shot = day+'r'+ str(x_shot)
    y_shot = day+'r'+ str(y_shot)
    z_shot = day+'r'+ str(z_shot)


    x_magdata = hdr.getquikData(x_shot)
    y_magdata = hdr.getquikData(y_shot)
    z_magdata = hdr.getquikData(z_shot)


    for probe in range(num_probes):
        #first make the B and V arrays:
        B = np.zeros((3,3))
        V = np.zeros((3,3))

        B[0,:], x_max_time= get_b(0, x_shot, probe)
        B[1,:], y_max_time= get_b(1, y_shot, probe)
        B[2,:], z_max_time= get_b(2, z_shot, probe)


        V[0,:] = get_v(x_magdata, probe, x_max_time)
        V[1,:] = get_v(y_magdata, probe, y_max_time)
        V[2,:] = get_v(z_magdata, probe, z_max_time)


        V_inv = np.linalg.inv(V)
        # take the dot product of the arrays
        C = np.matmul(B, V_inv)
        #store the flattened 3x3 array as a row in the text file
        calibration_matrix[probe] = np.reshape(C, 9)
        # print(calibration_matrix)
        # print(calibration_matrix[probe])

    pth = os.getcwd()
    date = '050119'
    # print(calibration_matrix.shape())
    savetxt(os.path.join(pth, 'calib-%s-4x4_matrix_1.txt' % (date)) , calibration_matrix)
    print("Finished! File saved as 'calib-%s-4x4_matrix.txt_1' in current working directory" %(date))

if __name__ == "__main__":
    generateMatrix()
