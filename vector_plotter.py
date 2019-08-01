from __future__ import division, print_function, absolute_import

# Import libraries:
import numpy as np
import scipy as sp
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
import matplotlib.pylab as pyl
import matplotlib.gridspec as gridspec
import numpy.ma as ma
import os
import time
# import imageio
from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d


#And from the SSX files:
import basic_hiresmag as hdr
import mjtools as mj #shouldn't use this
import ssx_data_read_basic as sdr
import process_mag_data_4x4 as pmd
import mytools as my

"""
This file contains a whole bunch of plotting functions, designed for plotting
the magnetic field in 2 or 3 d.

NOTE: You CANNOT just multiply a list by some value. You will get 100 lists
      Instead of that list with every element 100 times bigger. Hence the
      list comprehensions (one line nested for loops)

-- KG 06/14/19



"""

def get_color_frame(shot, t, num_probes = 16):
    """Not functional, but the goal was to have this be a function that would
       create the data froa frame of an animation used by some external program """
    # data = [[0] *3 for i in range(num_probes)]
    # magdata = hdr.getMagData(shot)
    # # print(len(magdata.time))
    # calibFile = 'calib-050119-4x4_lookup_1.txt'
    # calibFile = os.path.join(os.getcwd(), calibFile)
    # calibFactor = np.loadtxt(calibFile)
    data = pmd.get_data_slice(shot, t)
        # Get position vectors in meters
    position_vectors = get_probeLocs_meters()

    u = [data[i][0] for i in range(len(data))]
    v = [data[i][1] for i in range(len(data))]
    w = [data[i][2] for i in range(len(data))]


    x = [position_vectors[i][0] for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2] for i in range(len(position_vectors))]

    # qq = plt.quiver(y,z, v,w, u)
    # plt.plot()

    qq=plt.quiver(y,z, v,w,u,cmap=plt.cm.jet)
    scale_arrow = plt.quiver(.9,.9, .1, 0, label = ".1 Tesla")
    # plt.clim(-.004,.004)

    # plt.show()
    return plt, qq

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


def get_probeLocs_SSX_setup(num_probes = 16):
    """
    Provides positions in meters along probe stalks for 4x4 array of probes built by M. Kaur.
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
    x_pos = [2.85*1e-3*25.4, 1.4*1e-3*25.4, -1.4*1e-3*25.4, -2.85*1e-3*25.4]
    # y_pos = [4.25*1e-3*25.4, 4.25*1e-3*25.4, -4.24*1e-3*25.4, -4.24*1e-3*25.4]
    # z_pos = [-2.25*1e-3*25.4, -0.75*1e-3*25.4, 0.75*1e-3*25.4, 2.25*1e-3*25.4]
    z_pos = [-1*1e-3*25.4, -2.5*1e-3*25.4, -4.0*1e-3*25.4, -5.5*1e-3*25.4]
    x = 0
    for i in range(num_probes):
        if(i%4 ==0 and i>0):
            x+=1
        position_vectors[i][0] =x_pos[x]
        position_vectors[i][1] = 0
        position_vectors[i][2] =z_pos[i%4]
        # print(position_vectors[i][0])

    return position_vectors


def get_probeLocs_SSX_setup_cm(num_probes = 16):
    position_vectors = get_probeLocs_SSX_setup(num_probes)
    # position_vectors = [position_vectors[i]]
    position_vectors = [[element *100 for element in pos] for pos in position_vectors] #convert to cm
    return position_vectors


def get_probeLocs_meters(dir =2, num_probes = 16):
    """
    Probe locs for SSX Summer 2019 Setup are the same as the wire_tes,
    specified by this function!
    Provides positions in meters along probe stalks for 4x4 array of probes built by M. Kaur.
    ADL -- 20190501:
    -----
    This code generates the 16 by 3 array of positions. Each of the 16 rows
    repersents the position of a single probe locations, and the three
    elements are the x, y and z directions. Also worth noting that this uses
    a LEFT handed -cordinate system. This function returns porbe locations for
    probes in a square configuration
    -- KG  06/04/19
    """
    position_vectors = [[0] * 3 for i in range(num_probes)]

    #every x postion
    # y_pos = [ 4.24*1e-3, 1.41*1e-3, -1.42*1e-3, -4.25*1e-3]
    # x_pos = [-4.25*1e-3, -4.25*1e-3, 4.24*1e-3, 4.24*1e-3]
    # y_pos = [-4.25*1e-3, 4.24*1e-3, 4.24*1e-3, -4.25*1e-3]
    # 1 1/2, 2 1 1/2

    y_pos = [-2.5*1e-2*2.54,-1*1e-2*2.54,1*1e-2*2.54,2.5*1e-2*2.54]#cm
    z_pos = [-2.25*1e-2*2.54, -0.75*1e-2*2.54, 0.75*1e-2*2.54, 2.25*1e-2*2.54]#cm
    #
    # y_pos = [-2.5*2.54,-1*2.54,1*2.54,2.5*2.54]#cm
    # z_pos = [ -2.25*2.54, -0.75*2.54, 0.75*2.54, 2.25*2.54]##cm
    x = 0
    for i in range(num_probes):
        if(i%4 ==0 and i>0):
            x+=1
        position_vectors[i][0] = 0
        position_vectors[i][1] = y_pos[x]
        position_vectors[i][2] =z_pos[i%4]
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


def get_probeLocs_cm(dir =2, num_probes = 16):
    position_vectors = get_probeLocs_meters(num_probes)
    # position_vectors = [position_vectors[i]]
    position_vectors = [[element *100 for element in pos] for pos in position_vectors] #convert to cm
    return position_vectors



#
# def polyPeak_noPlot(time, data, timerange = [40,80],axis = 'x'):
#     """Finds the peak in the specified time range.
#     Finds the peak in the data in the specified time range.  You must pass it
#     pretrigger and timebase as well as the time and data series.
#     - timebase - is an integer equal to (1 us)/(delta t of the digitizer).  For
#       example, if capturing data at 10 MHz, timebase = (1e-6/1e-7) = 10.  For
#       data captured at 2 MHz, timebase = (1e-6/5e-7) = 2.  In other words, how
#       many data points per us.  It's a pretty bad way to code this up, but it
#       lets you specify the time range in micro seconds in which to look for the
#       peak.
#     - pretrigger - how many microseconds of data before t = 0.
#
#     -- M Kuar/T Gray
#
#     """
#
#     # Find the indices corresponding to the ends of the time range
#     t1 = mj.tindex(time,timerange[0])#+pretrigger)
#     t2 = mj.tindex(time,timerange[1])#+pretrigger)
#     # print 't1=', t1
#     # print 't2=', t2
#     # generate an array of indices spanning the range
#     ti = np.arange(t1,t2)
#     # get the time and data points in the range
#     t = time[ti]
#     d = data[ti]
#     # Fit a 2nd degree polynomial and find the min and max.
#     p = np.polyfit(t,d,2)
#     fit = p[0]*t**2 + p[1]*t + p[2]
#     dataMax = fit.max()
#     dataMin = fit.min()
#
#     if abs(dataMin) > dataMax:
#         dataMax = dataMin
#     return dataMax

def bfield_shape(shot, dir = 2, num_probes = 16):
    """
    Plot 16 vectors along the 48ch probe to
    sketch out shape of field at the midplane vs time.

    Default components are 'r' and 'z' in chamber coordinates.

    This will be the base for the movie I imagine.

    --KG 06/14/19
    """
    data = [[0] *3 for i in range(num_probes)]
    magdata = hdr.getMagData(shot)
    for probe_num in range(num_probes):
        Bx=sp.signal.detrend(magdata.fullData[0,probe_num,:])
        By=sp.signal.detrend(magdata.fullData[1,probe_num,:])
        Bz=sp.signal.detrend(magdata.fullData[2,probe_num,:])

        magmax_x = polyPeak_noPlot(magdata.time,Bx,[10,100], axis = 'x')
        magmax_y = polyPeak_noPlot(magdata.time,By,[10,100], axis = 'y')
        magmax_z = polyPeak_noPlot(magdata.time,Bz,[10,100], axis = 'z')

        data[probe_num] = [magmax_x, magmax_y, magmax_z]
    # position_vectors = get_probeLocs_meters(dir)
    position_vectors = get_probeLocs_calib_setup(2)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = [position_vectors[i][0] for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2] for i in range(len(position_vectors))]

    u = [data[i][0] for i in range(len(data))]
    v = [data[i][1] for i in range(len(data))]
    w = [data[i][2] for i in range(len(data))]


    dead_chans = [1, 0,0,0,0,0,0,0,1,0,1,0,1,0,0,0]
    good_u = np.ma.masked_array(u,dead_chans)
    good_v = np.ma.masked_array(v, dead_chans)
    good_w = np.ma.masked_array(w, dead_chans)

    dead_u = np.ma.masked_array(u, np.logical_not(dead_chans))
    dead_v = np.ma.masked_array(v, np.logical_not(dead_chans))
    dead_w = np.ma.masked_array(w, np.logical_not(dead_chans))





    col_map = plt.get_cmap('Blues')
    ax.quiver(x,y,z,good_u,good_v,good_w, length=0.01, normalize=True, cmap = col_map)
    ax.quiver(x,y,z,dead_u,dead_v,dead_w, length=0.01, normalize=True, color = 'r')


    ax.set_xlabel("X direction")
    ax.set_ylabel("Y direction")
    ax.set_zlabel("Z direction")
    ax.set_zlim(-.05, 0.05)
    ax.set_xlim(-.05, 0.05)
    ax.set_ylim(-.05, 0.05)
    plt.axis('equal')
    # plt.xlim()
    # ax.set_title("Stright Wire test")
    plt.show()

def generate_cyl_dat(center_x = 0,center_y= 0 ,radius = 0.08,height_z = 0.095):
    z = np.linspace(-height_z, height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    # return a cylneder rotated to match the coordinate system that everything
    # else is plotted in
    return z_grid, x_grid, y_grid


def plot_color(shot, shot_num, num_probes = 16):
    """ This code plot a 2-d graph, with colored arrows to indicated the third
        dimension. For example, in the straight wire test, this plots in the
        y-z plane since that's where most of the field is, and then has colored
        arrows to indicate how much field is in the x-directionself.
        TODO: add a key to show how the length of the arrow corresponds to
              to field strength

              also just worth noting and maybe worth fixing is that the color
              bar limits are set to be on the same axis as the y and z, but
              isn't actually the strength of the field. It was just to make the
              relative field equal to the other components

       -- KG 06/11/19
    """
    data = [[0] *3 for i in range(num_probes)]
    magdata = hdr.getMagData(shot)
    for probe_num in range(num_probes):
        Bx=sp.signal.detrend(magdata.fullData[0,probe_num,:])
        By=sp.signal.detrend(magdata.fullData[1,probe_num,:])
        Bz=sp.signal.detrend(magdata.fullData[2,probe_num,:])


        magmax_x = polyPeak_noPlot(magdata.time,Bx,[10,100], axis = 'x')
        magmax_y = polyPeak_noPlot(magdata.time,By,[10,100], axis = 'y')
        magmax_z = polyPeak_noPlot(magdata.time,Bz,[10,100], axis = 'z')


        data[probe_num] = [magmax_x, magmax_y, magmax_z]


    position_vectors = get_probeLocs_meters()
    x = [position_vectors[i][0] for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2] for i in range(len(position_vectors))]

    u = [data[i][0] for i in range(len(data))]
    v = [data[i][1] for i in range(len(data))]
    w = [data[i][2] for i in range(len(data))]
    # v[4] = v[4]*-1
    # v[5] = v[5]*-1
    v[6] = 0


    x = [position_vectors[i][0] for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2] for i in range(len(position_vectors))]


    qq=plt.quiver(y,z, v,w,u,cmap=plt.cm.jet)
    plt.clim(-.004,.004)
    plt.colorbar(qq, cmap=plt.cm.jet)
    plt.xlabel("Y direction")
    plt.ylabel("Z direction")
    plt.xlim(-.004,.004)
    plt.ylim(-.004,.004)
    plt.title("Magnetic Data - Shot %s" %(shot_num))
    plt.show()

def plot_2d_offCenter_down(shot,  num_probes = 16, lim = 8.5):
    """ This is mainly used for plotting the straight wire test. In this
        test the field was mainly in the y-z plane, so it plots only the
        y-z components and then overlays calculated magnetic field using
        B = mu_naught I/ 2 pi r. It plots the maximum magentic field, and
        uses the maximum current for plotting. The positions are in meters,
        so the units of B are tesla.The off-center, down refers to the fact
        that the wire was not exactly centered in the probe array, so the
        coordinate sytem was shifted accordingly.

        --KG 07/25/19
    """

    timeB, data = pmd.get_data(shot, num_probes)

    #find the time when the current is at it's max:
    time,eastcurrent,westcurrent = loadcurrent(shot)
    index = np.argmax(eastcurrent)
    maxt = time[index]

    t_index = my.tindex_center(timeB,maxt)
    u = [data[i][0][t_index] for i in range(num_probes)]
    v = [data[i][1][t_index] for i in range(num_probes)]
    w = [data[i][2][t_index] for i in range(num_probes)]


    position_vectors = get_probeLocs_cm()
    x = [position_vectors[i][0] for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2] for i in range(len(position_vectors))]
    z = [i - (0.5*25.4*1e-3) for i in z] #shift everthing back .5 inch


    qq = plt.quiver(y,z, v,w)
    # plt.plot()
    # scale_arrow = plt.quiver(-lim+ .01, -lim+ 0.01, .1, 0, label = ".1 Tesla")
    plt.quiverkey(qq, .92,.92, .1, label = ".1 Tesla")
    plt.plot()

    mu = 2 * 10**-7
    print("The maximum current is: %s" %(eastcurrent[index]))


    # lim = .008

    xlist = np.linspace(-lim-1.3, lim+1.3, 1000)
    ylist = np.linspace(-lim-1.3,lim+1.3, 1000)
    X, Y = np.meshgrid(xlist, ylist)
    Z = np.sqrt(X**2 + Y**2)
    B = mu*maxI * (Z**(-1))

    font_size = 18
    plt.rcParams.update({'font.size': font_size})
    levels = [0, 0.001, 0.05, 0.1, 0.2, .3, .4,0.5, .6, .7, .8, .9,1]
    contour = plt.contour(X, Y, B, levels, alpha = 0.2, colors='k')
    # contour_filled = plt.contourf(X, Y, B, levels, alpha = 0.2, cmap=plt.cm.jet)
    level_colors =('k', 'navy', 'b','dodgerblue', 'c', 'mediumseagreen', 'yellowgreen', 'gold', 'orange', 'orangered', 'tomato', 'crimson', 'r')
    contour_filled = plt.contourf(X, Y, B, levels, alpha = 0.2, colors = level_colors)
    # cbar = plt.colorbar(contour_filled)
    plt.ylabel('T')

    # contour_filled.cmap.set_under('yellow')
    contour_filled.cmap.set_over('r')
    # contour_filled.c

    cbar = plt.colorbar(contour_filled)
    cbar.set_label("Magnetic Field Strength (T)", fontsize = font_size)
    # plt.legen:d()
    plt.xlabel("Y direction (m)", fontsize = font_size)
    plt.ylabel("Z direction (m)", fontsize = font_size)
    # plt.axis('equal')
    plt.xlim(-lim, lim)
    plt.ylim(-lim- 1.3, lim-1.3)
    plt.title("Straight Wire Test, I = 14.3 kA, with Magnetic Field Strength (T)")
    plt.show()


def plot_2d(shot,  num_probes = 16):
    '''Plots a black and white plot of the magentic data for a given shot
    '''
    data = [[0] *3 for i in range(num_probes)]
    magdata = hdr.getMagData(shot)
    for probe_num in range(num_probes):
        Bx=sp.signal.detrend(magdata.fullData[0,probe_num,:])
        By=sp.signal.detrend(magdata.fullData[1,probe_num,:])
        Bz=sp.signal.detrend(magdata.fullData[2,probe_num,:])

        magmax_x = polyPeak_noPlot(magdata.time,Bx,[10,100], axis = 'x')
        magmax_y = polyPeak_noPlot(magdata.time,By,[10,100], axis = 'y')
        magmax_z = polyPeak_noPlot(magdata.time,Bz,[10,100], axis = 'z')

        data[probe_num] = [magmax_x, magmax_y, magmax_z]

    # Get position vectors in meters
    position_vectors = get_probeLocs_meters()

    # x = np. a
    x = [position_vectors[i][0] for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2] for i in range(len(position_vectors))]


    u = [data[i][0] for i in range(len(data))]
    v = [data[i][1] for i in range(len(data))]

    # v[6] = 0
    w = [data[i][2] for i in range(len(data))]


    x = [position_vectors[i][0] for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2] for i in range(len(position_vectors))]

    qq = plt.quiver(y,z, u,w)
    plt.plot()

    # print(position_vectors)
    lim = .08


    plt.xlabel("X direction (m)")
    plt.ylabel("Z direction (m)")

    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    plt.title("Straight Wire Test")
    plt.show()


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


    # calculate the be field.  The 4pi's cancel out, and the 1e4 is to convert
    # to gauss
    b = b * i * 1.0e-7 * turns * 1e4
    return b


def load_smallbdot(shot, scope = '1'):
    """ Little helper function to load in the B-dot """
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


def loadcurrent(shot, scope = '3'):
    """ helper function to read in the current"""
    data = sdr.scope_data(shot, scope)
    time = data.time
    eastcurrent = -data.ch2*170000#East gun 170kA/V
    westcurrent = -data.ch4*190000#West gun 190kA/V
    return time,eastcurrent,westcurrent


def plot_voltage(shot, scope = '3'):
    """ helper function to read in the current"""
    data = sdr.scope_data(shot, scope)
    time = data.time
    eastvolt= data.ch1#East gun 170kA/V
    westvolt = data.ch3#West gun 190kA/V
    plt.plot(time, eastvolt, label = "East Voltage")
    plt.plot(time, westvolt, label = "West Voltage")
    plt.legend()
    plt.show()



def plots_4_doc(shot, num_probes = 16):
    """ This plots the magentic field from every single probe one at a timeself.
        The advantage is that is also calculates what the helmholtz field at
        that point should be and plots it (as a function of time) and also
        current all one one plot, making it easy to find the ratio, and also
        verify it vizually. However, to make this function run in a reasonable
        amount of time, I thinned the data for the calulated helmholtz field,
        so that it only calculates every 10th point, meaning it is less
        accurate at actually producing the ratio (becuase it may not find the
        true peak )


        TODO:
            - check if the units for the third plot are ok???
            - Add y-units


        -- KG 06/13/19"""
    # data = [[0] *3 for i in range(num_probes)]
    magdata = hdr.getMagData(shot)

    probe_locs = get_probeLocs_calib_setup(shot)
    for probe_num in range(num_probes):
        # Bz=sp.signal.detrend(magdata.fullData[2,probe_num,:])
        # Bx=sp.signal.detrend(magdata.fullData[0,probe_num,:])
        # By=sp.signal.detrend(magdata.fullData[1,probe_num,:])
        ratio = 1
        timeb, Bx, By, Bz = load_smallbdot(shot)

        data=hdr.getquikData(shot)
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
        # Bz=cumtrapz(data.unCalibData[2,1,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[2,1,:], 1), data.time)
        Bz=sp.signal.detrend(cumtrapz(data.unCalibData[1,probe_num,:], data.time))
        if(abs(np.max(Bz)) > abs(np.min(Bz))):
            print("Max %d, min %d, abs: %d" %(np.max(Bz), np.min(Bz), abs(np.min(Bz))))
            Bz =  Bz* -1
            ratio = -1

        Bz += np.abs(np.min(Bz))
        btime = data.time[:len(Bz)]
        ax1.plot(btime, Bz)
        ax1.set_title("By - integrated")
        ax1.set_xlim(-25, 175)

        time,eastcurrent,westcurrent = loadcurrent(shot)
        ax2.plot(time[::10], eastcurrent[::10])
        ax2.set_title("Current")

        r = probe_locs[1]
        # helmB = helmholtz_listcomp(r,eastcurrent)
        status = 0
        max = len(eastcurrent)/100
        helmB = np.zeros((3, max))
        i = 0
        while i+100 < len(eastcurrent):
        # currmax = polyPeak_noPlot(time,eastcurrent,[10,100],axis = 'current')
            x,y,z = helmholtz2(r,eastcurrent[i])
            helmB[0][status] = x
            helmB[1][status] = y
            helmB[2][status] = z
            status +=1
            i+=100
        # thin the thing
        thinned_time = time[::100]
        thinned_time = thinned_time[:len(helmB[2])]
        helmB[2] = helmB[2]*-1
        max_measured = np.max(Bz)
        max_theoretical = np.max(helmB[2])
        ratio = ratio * max_theoretical/max_measured


        print("Shot: %s, probe %d, Ratio around: %d" %(shot[:6], probe_num, ratio))
        ax3.plot(thinned_time, helmB[2])
        ax3.set_title("Predicted B based on current and location in helmholtz")
        print("Ratio is: %d" %(ratio))
        plt.show()

def plot_both_currents(shot):
    """ Quick plotting routine that
        plots both currents (east and west)
        so you can double check the timeing delay

        -- KG 06/27/19"""
    time,eastcurrent,westcurrent = loadcurrent(shot)

    # fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    # ax1.plot(time[::10], eastcurrent[::10], label = 'eastcurrent')
    # ax2.plot(time[::10], westcurrent[::10], label = 'westcurrent')

    # ax1.set_title("Current")
    # ax2.set_xlabel("Time $\\mu$s")

    plt.plot(time[::10], eastcurrent[::10], label = 'eastcurrent')
    plt.plot(time[::10], westcurrent[::10], label = 'westcurrent')
    plt.legend()
    plt.title("Shot %s \nCurrents vs Time" %shot)
    plt.show()


def integrated_B(shot, num_probes = 16, probe = 3):
    """
       There are different way that I've seen the data
       integrated: this will plot the different ones, along with
       a variety of different moving averages

        --KG 06/13/19
        """

    # data = [[0] *3 for i in range(num_probes)]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    data=hdr.getquikData(shot)
    btime = data.time
    Bz_nBT_15 =cumtrapz(data.unCalibData[2,probe,:]-mj.get_gaussian_moving_avg(data.time,
            data.unCalibData[2,probe,:], 15), data.time)* -1.103810335732411829e+02

    Bz_nBT_100 =cumtrapz(data.unCalibData[2,probe,:]-mj.get_gaussian_moving_avg(data.time,
            data.unCalibData[2,probe,:], 100), data.time)* -1.103810335732411829e+02

    Bz_nBT_500 =cumtrapz(data.unCalibData[2,probe,:]-mj.get_gaussian_moving_avg(data.time,
            data.unCalibData[2,probe,:], 500), data.time)* -1.103810335732411829e+02
    # num_probes = 1 #only plot the first probe
    # probe_locs = get_probeLocs_calib_setup(shot)
        # Bz=sp.signal.detrend(magdata.fullData[2,probe_num,:])
        # Bx=sp.signal.detrend(magdata.fullData[0,probe_num,:])
    # By=sp.signal.detrend(magdata.fullData[1,probe_num,:])

    # Bz=cumtrapz(data.unCalibData[2,1,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[2,1,:], 1), data.time)
    Bz=sp.signal.detrend(cumtrapz(data.unCalibData[2,probe,:], data.time))
    btime = btime[:len(Bz)]
    ax1.plot(btime, Bz, label= 'Detrend, integrated')
    ax1.plot(btime, Bz_nBT_15, label= '-moving avg, win = 15')
    ax1.plot(btime, Bz_nBT_100, label= '-moving avg, win = 100')
    ax1.plot(btime, Bz_nBT_500, label= '-moving avg, win = 500')
    ax1.legend()

    Bdot25 = np.zeros([3,25,8192])
    timeB=data.time[801:]
    Bdot25[0,3,:]=data.Bdot[0,3,:]
    Bdot25[1,3,:]=data.Bdot[1,3,:]
    Bdot25[2,3,:]=data.Bdot[2,3,:]
    bzero = pyl.detrend_linear(Bdot25[2,3,800:])
    bint = sp.integrate.cumtrapz(bzero,data.time[800:])
    ax2.plot(btime, Bz_nBT_15, label= '-moving avg, win = 15')
    ax2.plot(btime, Bz_nBT_500, label= '-moving avg, win = 500')
    ax2.plot(timeB, bint, label = "detrend LINEAR and integrated time shifted")
    # ax2.set_title("Current")

    plt.legend()
    plt.show()
#
    # fig = plt.figure() ; ax = fig.gca(projection='3d') ; qq=ax.quiver(X, Y, Z, dx, dy, dz, mag,cmap=plt.cm.jet) ; ax.set_title('3D Vector Field') ; ax.view_init(elev=18, azim=30) ; ax.dist=8 ; plt.colorbar(qq, cmap=plt.cm.jet) ; plt.show()

def plot_every_single_UNcalib_chan(shot, dir, lim, save = False, show = True, num_probes = 16):
    """This function will plot a 4x4 (though there is some funcitonality built
        in for a different number of probes) grid of the CALIBRATED,
        INTEGRATED signals. Unless the calibration value was zero
        (when I built this some of the probes were dead during calibration)
        In which case it will plot the UNCALIBRATED (INTEGRATED) DATA IN RED
        for un-dead probes. Give it the shot (date + r + shot number) and
        the direction -
        0 for x
        1 for y
        2 for z
        and it will plot all the signals.

        If you want to mess around with the number of probes, try passing it in
        as the parameter, though you may need to adjust cols or rows or both,
        and also maybe the figure size, (make is smaller if you have more probes)

        --KG 06/13/19
        """


    if(dir == 0):
        direction = 'X'
    elif(dir == 1):
        direction = 'Y'
    else:
        direction = 'Z'

    figsize = (9, 5)
    cols = 4
    rows = (num_probes) // cols
    # define the data for cartesian plots


    def trim_axs(axs, N):
        """helper function to make sure the axs list to have correct length"""
        axs = axs.flat
        for ax in axs[N:]:
            ax.remove()
        return axs[:N]

    fig1, axs = plt.subplots(rows, cols, figsize=figsize)
    axs = trim_axs(axs, num_probes)
    # for ax, case in zip(axs, cases)

    magdata=hdr.getquikData(shot)
    for i in range(len(axs)):
        # y = cumtrapz(magdata.unCalibData[dir,i,:]-mj.get_gaussian_moving_avg(magdata.time, magdata.unCalibData[dir,i,:], 15), magdata.time)
        y = magdata.unCalibData[dir,i,:]
        x = magdata.time
        x = x[:len(y)]
        ax = axs[i]
        ax.set_title('Probe (%s,%s %s)' %(str(i//4+ 1), str(i%4 +1), direction), fontsize = 10)

        if np.average(np.abs(y)) < lim*1.5e-2:
            ax.plot(x, y, color = 'r')
            ax.set_ylim(-lim*1e-2, lim*1e-2)
        else:
            ax.plot(x, y)
            ax.set_ylim(-lim, lim)

        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.tick_params(axis='both', which='minor', labelsize=7)
        if(i%cols == 0 ): #lable all the left cols with y label
            ax.set_ylabel("Raw Voltage (V)", fontsize = 8)
        if(i//rows == rows-1 ): #lable all the left cols with y label
            ax.set_xlabel("Time $\\mu$s", fontsize = 8)
        # plt.title("Shot %s" %(shot))
    if show:
        plt.show()
    if save:
        try:
            figname = os.getcwd() + '\\generated_images\\'+ str(shot) + "-"+ direction + "-chans-uncalib.png"
            # print("Saved fig: %s" %(figname))
            plt.savefig(figname)
        except:
            figname = os.getcwd() + str(shot) + "-"+ direction + "-chans-uncalib.png"
            plt.savefig(figname)
            print("No generated_images folder. \nSaved fig: %s" %(figname))

def plot_every_single_calib_chan(shot, dir, lim, save = False, show = True, num_probes = 16):
    """This function will plot a 4x4 (though there is some funcitonality built
        in for a different number of probes) grid of the CALIBRATED,
        INTEGRATED signals. Unless the calibration value was zero
        (when I built this some of the probes were dead during calibration)
        In which case it will plot the UNCALIBRATED (INTEGRATED) DATA IN RED
        for un-dead probes. Give it the shot (date + r + shot number) and
        the direction -
        0 for x
        1 for y
        2 for z
        and it will plot all the signals.

        If you want to mess around with the number of probes, try passing it in
        as the parameter, though you may need to adjust cols or rows or both,
        and also maybe the figure size, (make is smaller if you have more probes)

        --KG 06/13/19
        """


    if(dir == 0):
        direction = 'X'
    elif(dir == 1):
        direction = 'Y'
    else:
        direction = 'Z'

    figsize = (9, 5)
    cols = 4
    rows = (num_probes) // cols
    # define the data for cartesian plots


    def trim_axs(axs, N):
        """helper function to make sure the axs list to have correct length"""
        axs = axs.flat
        for ax in axs[N:]:
            ax.remove()
        return axs[:N]

    calibFile = 'calib-050119-4x4_lookup_1.txt'
    calibFile = os.path.join(os.getcwd(), calibFile)
    calibFactor = np.loadtxt(calibFile)

    fig1, axs = plt.subplots(rows, cols, figsize=figsize)
    axs = trim_axs(axs, num_probes)
    # for ax, case in zip(axs, cases)

    magdata=hdr.getquikData(shot)
    # Bx,By, Bz = pmd.getMagData(shot)
    for i in range(len(axs)):
        y = cumtrapz(magdata.unCalibData[dir,i,:]-mj.get_gaussian_moving_avg(magdata.time, magdata.unCalibData[dir,i,:], 15), magdata.time)

        # print(calibFactor[i][dir])
        x = magdata.time
        x = x[:len(y)]
        ax = axs[i]
        ax.set_title('Probe (%s,%s %s)' %(str(i//4+ 1), str(i%4 +1), direction), fontsize = 10)
        if(calibFactor[i][dir] == 0):
            # If the probe was dead during calibration
            ax.plot(x, y, color = 'r', linestyle = 'dashed')
        else:
            y = y * calibFactor[i][dir]
            ax.plot(x, y)
            ax.set_ylim(-lim, lim)
        # if np.average(np.abs(y)) < lim*1e-2:
        #     ax.plot(x, y, color = 'r')
        #     ax.set_ylim(-lim*1e-5, lim*1e-5)
            # print(i, np.max(y), np.min(y))



        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.tick_params(axis='both', which='minor', labelsize=7)
        if(i%cols == 0 ): #lable all the left cols with y label
            ax.set_ylabel("Gauss", fontsize = 8)
        if(i//rows == rows-1 ): #lable all the left cols with y label
            ax.set_xlabel("Time $\\mu$s", fontsize = 8)
        # plt.title("Shot %s" %(shot))
    if show:
        plt.show()
    if save:
        try:
            figname = os.getcwd() + '\\generated_images\\'+ str(shot) + "-"+ direction + "-chans.png"
            # print("Saved fig: %s" %(figname))
            plt.savefig(figname)
        except:
            figname = os.getcwd() + str(shot) + "-"+ direction + "-chans.png"
            plt.savefig(figname)
            print("No generated_images folder. \nSaved fig: %s" %(figname))

def bfield_movie(shot, dir = 2, num_probes = 16):
    """
    based on bfield_shape, but this will trace out the magnetic feild though
    time?


    --KG 06/14/19
    """
    data = [[0] *3 for i in range(num_probes)]
    magdata = hdr.getMagData(shot)
    for probe_num in range(num_probes):
        Bx=sp.signal.detrend(magdata.fullData[0,probe_num,:])
        By=sp.signal.detrend(magdata.fullData[1,probe_num,:])
        Bz=sp.signal.detrend(magdata.fullData[2,probe_num,:])

        magmax_x = polyPeak_noPlot(magdata.time,Bx,[10,100], axis = 'x')
        magmax_y = polyPeak_noPlot(magdata.time,By,[10,100], axis = 'y')
        magmax_z = polyPeak_noPlot(magdata.time,Bz,[10,100], axis = 'z')

        data[probe_num] = [magmax_x, magmax_y, magmax_z]
    # position_vectors = get_probeLocs_meters(dir)
    position_vectors = get_probeLocs_calib_setup(2)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = [position_vectors[i][0] for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2] for i in range(len(position_vectors))]

    u = [data[i][0] for i in range(len(data))]
    v = [data[i][1] for i in range(len(data))]
    w = [data[i][2] for i in range(len(data))]


    dead_chans = [1, 0,0,0,0,0,0,0,1,0,1,0,1,0,0,0]
    good_u = np.ma.masked_array(u,dead_chans)
    good_v = np.ma.masked_array(v, dead_chans)
    good_w = np.ma.masked_array(w, dead_chans)

    dead_u = np.ma.masked_array(u, np.logical_not(dead_chans))
    dead_v = np.ma.masked_array(v, np.logical_not(dead_chans))
    dead_w = np.ma.masked_array(w, np.logical_not(dead_chans))





    col_map = plt.get_cmap('Blues')
    ax.quiver(x,y,z,good_u,good_v,good_w, length=0.01, normalize=True, cmap = col_map)
    ax.quiver(x,y,z,dead_u,dead_v,dead_w, length=0.01, normalize=True, color = 'r')


    ax.set_xlabel("X direction")
    ax.set_ylabel("Y direction")
    ax.set_zlabel("Z direction")
    ax.set_zlim(-.05, 0.05)
    ax.set_xlim(-.05, 0.05)
    ax.set_ylim(-.05, 0.05)
    plt.axis('equal')
    # plt.xlim()
    # ax.set_title("Stright Wire test")
    plt.show()


def generate_3d_frame(shot, t, count, num_probes = 16):
    """ saves a frame instead of returning it """
    plt.close("all")
    data = [[0] *3 for i in range(num_probes)]
    magdata = hdr.getMagData(shot)
    # print(len(magdata.time))
    calibFile = 'calib-050119-4x4_lookup_1.txt'
    calibFile = os.path.join(os.getcwd(), calibFile)
    calibFactor = np.loadtxt(calibFile)

    for probe in range(num_probes):
        Bx_all=cumtrapz(magdata.unCalibData[0,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[0,probe,:], num_probes-1), magdata.time)
        Bx = Bx_all[int(t)] * calibFactor[probe][0]

        By_all=cumtrapz(magdata.unCalibData[1,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[1,probe,:], num_probes-1), magdata.time)
        By = By_all[int(t)] * calibFactor[probe][1]

        Bz_all=cumtrapz(magdata.unCalibData[2,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[2,probe,:], num_probes-1), magdata.time)
        Bz = Bz_all[int(t)] * calibFactor[probe][2]



        data[probe] = [Bx, By, Bz]

        # Get position vectors in meters
    position_vectors = get_probeLocs_SSX_setup_cm()

    u = [data[i][0] for i in range(len(data))]
    v = [data[i][1] for i in range(len(data))]
    w = [data[i][2] for i in range(len(data))]

    # too big, so I just zero it out
    # v[6] = 0


    """
    In SSX setup:

        x -> -z
        y -> theta
        z -> -r

    """
    x = [position_vectors[i][0]*-1 for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2]*-1 for i in range(len(position_vectors))]

    # qq = plt.quiver(y,z, v,w, u)
    # plt.plot()

    iffy_chans = [1, 0,0,0,0,0,0,0,1,0,1,0,1,0,0,0]
    good_u = np.ma.masked_array(u,iffy_chans)
    good_v = np.ma.masked_array(v, iffy_chans)
    good_w = np.ma.masked_array(w, iffy_chans)

    dead_chans = [1, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    dead_u = np.ma.masked_array(u, np.logical_not(dead_chans))
    dead_v = np.ma.masked_array(v, np.logical_not(dead_chans))
    dead_w = np.ma.masked_array(w, np.logical_not(dead_chans))

    fixed_chans = np.logical_xor(iffy_chans, dead_chans)
    fixed_u = np.ma.masked_array(u, np.logical_not(fixed_chans))
    fixed_v = np.ma.masked_array(v, np.logical_not(fixed_chans))
    fixed_w = np.ma.masked_array(w, np.logical_not(fixed_chans))



    fig = plt.figure()
    ax = fig.gca(projection='3d')
    col_map = plt.get_cmap('Blues')

    # ax.quiver(x,y,z,good_u,good_v,good_w, length=0.01, normalize=False,cmap = col_map, label = "Calibrated probes" )
    # ax.quiver(x,y,z,dead_u,dead_v,dead_w, length=0.01, normalize=False, color = 'r', label = "Broken components")
    # ax.quiver(x,y,z,fixed_u,fixed_v,fixed_w, length=0.01, normalize=False,color = 'g', label = "interpolated components")

    # Q = ax.quiver(x,y,z,good_u,good_v,good_w, length=0.05,  normalize=False, pivot = 'middle', cmap = col_map)
    # ax.quiver(x,y,z,dead_u,dead_v,dead_w, length=0.05,  normalize=False, pivot = 'mid',  color = 'r')
    # ax.quiver(x,y,z,fixed_u,fixed_v,fixed_w, length=0.05,  normalize=False, pivot = 'mid', color = 'g')

    ax.quiver3D(x,y,z,good_u,good_v,good_w, length=0.02,arrow_length_ratio = .4, pivot = 'middle',   normalize=False, cmap = col_map)
    ax.quiver3D(x,y,z,dead_u,dead_v,dead_w, length=0.02, arrow_length_ratio = .4, pivot = 'middle',  normalize=False,  color = 'r')
    ax.quiver3D(x,y,z,fixed_u,fixed_v,fixed_w, length=0.02, arrow_length_ratio = .4, pivot = 'middle',  normalize=False, color = 'g')

    # qk = Q.quiverkey(Q, .92,.92, .1, label = ".1 Tesla")
    # ax.plot(np.zeros(100),np.zeros(100),'cyan--', alpha = .3)


    #set up the  pivot = 'mid',midline
    N = 100
    lim = 8.5

    #workaround quiverkey:
    ax.quiver(0, 0 , 0, 100, 0, 0,length=0.02,arrow_length_ratio = .4,  pivot = 'middle',color = 'c', label = '100 Gauss')

    # ax.plot(np.linspace(-lim, lim, N), np.zeros(N), 'k--', alpha  = 0.6, label = "midline")
    ax.plot(np.linspace(-lim, lim, N), np.zeros(N), 'k--', alpha  = 0.6)
    plt.legend()


    # add the sca;le
    # plt.legend(loc= 'upper right')


    #add cylinder:
    # Xc,Yc,Zc = generate_cyl_dat()
    # Xc[Xc>lim]= np.nan
    # Xc[Xc<-lim]= np.nan
    # Yc[Yc>lim]= np.nan
    # Yc[Yc<-lim]= np.nan
    # Zc[Zc>lim]= np.nan
    # Zc[Zc<-lim]= np.nan
    # ax.plot_surface(Xc, Yc, Zc,  color = 'tan', alpha=0.3)

    #add edge of SSX
    # z = np.empty(N)
    # z.fill(0.08025)
    # x = np.linspace(-lim, lim, N)

    # X, Z = np.meshgrid(x,z)
    # Y =


    # ax.plot(np.linspace(-lim, lim, N), z, 'y--', alpha  = 0.6, label = "midline")
    # ax.plot_surface(,   z, 'y--', alpha  = 0.6, label = "edge of SSX")


    #add title/labels
    ax.set_xlabel("$\\hat{Z}$, cm")
    ax.set_ylabel("$\\hat{\\theta}, cm$")
    ax.set_zlabel("$\\hat{R}, cm$")
    ax.set_zlim3d(-lim, lim)
    ax.set_xlim3d(-lim, lim)
    ax.set_ylim3d(-lim, lim)
    # plt.zlim(-lim, lim)
    # plt.xlim(-lim, lim)
    # plt.ylim(-lim, lim)
    plt.axis('equal')
    # plt.title("Magnetic Data - Shot %s. \n $\\bf{Time}$: %.2f $\\it{\mu s}$ " %(shot[7:], t))
    plt.title("Magnetic Data - Shot %s. \n $\\bf{Time}$: %.2f $\\it{\mu s}$ " %(shot, t))
    path = os.getcwd()
    path = path + '/mag_images/'

    #quick fix to make sure the ordering makes sense:
    str_count = str(count)
    if count <10:
        str_count = '0' + str_count

    if count <100:
        str_count = '0' + str_count

    if count <1000:
        str_count = '0' + str_count

    path = path + str_count +'-mag'  + '.png'

    plt.savefig(path)
    # plt.show()


def generate_all_3d_frames(shot, thin_every = 1):
    """ Wrapper function to generate a frame for every timestep"""
    magdata = hdr.getMagData(shot)
    # speed things up
    count = 0
    for t in magdata.time[::thin_every]:
    # for t in magdata.time:
        print("Progress: %s / %s" %(str(count), str(len(magdata.time)//thin_every)))
        if(t > 20 and t< 150):
            generate_3d_frame(shot, t, count)

        count += 1

def generate_color_frame(shot, t, count, num_probes = 16):
    plt.close("all")
    """ saves a frame instead of returning it ! """
    data = [[0] *3 for i in range(num_probes)]
    magdata = hdr.getMagData(shot)
    # print(len(magdata.time))
    calibFile = 'calib-050119-4x4_lookup_1.txt'
    calibFile = os.path.join(os.getcwd(), calibFile)
    calibFactor = np.loadtxt(calibFile)

    for probe in range(num_probes):
        Bx_all=cumtrapz(magdata.unCalibData[0,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[0,probe,:], num_probes-1), magdata.time)
        Bx = Bx_all[int(t)] * calibFactor[probe][0]

        By_all=cumtrapz(magdata.unCalibData[1,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[1,probe,:], num_probes-1), magdata.time)
        By = By_all[int(t)] * calibFactor[probe][1]

        Bz_all=cumtrapz(magdata.unCalibData[2,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[2,probe,:], num_probes-1), magdata.time)
        Bz = Bz_all[int(t)] * calibFactor[probe][2]



        data[probe] = [Bx, By, Bz]

        # Get position vectors in meters
    position_vectors = get_probeLocs_SSX_setup()

    u = [data[i][0] for i in range(len(data))]
    v = [data[i][1] for i in range(len(data))]
    w = [data[i][2] for i in range(len(data))]

    # too big, so I just zero it out
    # v[6] = 0


    x = [position_vectors[i][0] for i in range(len(position_vectors))]
    y = [position_vectors[i][1] for i in range(len(position_vectors))]
    z = [position_vectors[i][2] for i in range(len(position_vectors))]

    # qq = plt.quiver(y,z, v,w, u)
    # plt.plot()

    qq=plt.quiver(y,z, v,w,u,cmap=plt.cm.jet)
    scale_arrow = plt.quiver(.9,.9, .1, 0, label = ".1 Tesla")
    cbar = plt.colorbar(qq, cmap=plt.cm.jet)
    cbar.set_label("Magnetic Field Strength (T)")
    # cbar.set_norm((-.085,.085))
    cbar.set_clim(-.085,.085)
    plt.xlabel("Y direction")
    plt.ylabel("Z direction")
    plt.xlim(-.085,.085)
    plt.ylim(-.085,.085)
    plt.title("Magnetic Data - Shot %s. Time: %s" %(shot, str(t)))
    path = os.getcwd()
    path = path + '/mag_images/'
    path = path + 'mag-' + str(count) + '.png'

    plt.savefig(path)


def generate_all_color_frames(shot, thin_every =1):
    """ Wrapper function to generate a frame for every timestep"""
    magdata = hdr.getMagData(shot)
    # speed things up
    count = 0
    for t in magdata.time[::thin_every]:
        generate_color_frame(shot, t, count)
        count +=1

def convert_to_gif(files = " "):
    pth = os.getcwd()
    path = pth + '/mag_images/'
    if files == " ":
        files = os.listdir(path)
    # images = [imageio.imread(filename) for filename in fileLocs ]
    images = []
    for filename in files:
        images.append(imageio.imread(path + filename))

    save_name = 'movie_test.gif'
    imageio.mimsave(save_name, images)
    print("SAVED as %s" %(save_name))


def convert_to_mp4(num_runs = 1, files = " "):
    print("Converting to mp4")
    pth = os.getcwd()
    path = pth + '/mag_images/'

    if files == " ":
        files = os.listdir(path)
        files.sort()
    images = []
    # images = [imageio.imread(path + filename) for filename in files]
    for filename in files:

        for i in range(num_runs):
            images.append(imageio.imread(path + filename))

    save_name = 'movie_test.mp4'
    imageio.mimsave(save_name, images)
    print("Finished! Saved as %s" %(save_name))



def main():
    day = '073019'
    # day  = '062619'
    shot_num =13
    shot = day+'r'+ str(shot_num)
    # thin_every = 1 #for default
    thin_every = 1000
    # integrated_B(shot)
    # create_txt_stalk(shot)
    # generate_all_3d_frames(shqot, thin_every) #can take a very long time to run
    # convert_to_mp4(5)
    # convert_to_gif()
    # plt = get_color_frame(shot, 0)
    # plt.show()
    # shot ='23'#y dir
    # shot 12 in z dir
    # bfield_shape(shot)
    # plot_color(shot, shot_num)
    # plot_2d(shot)
    # shot ='3'
    # shot = day+str(shot)
    # plot_2d_offCenter_down(shot)
    plot_voltage(shot)
    # plots_4_doc(shot)
    # plot_both_currents(shot)
    # integrated_B(shot)

    # plot_every_single_UNcalib_chan(shot, dir = 2, lim = 2, save = False)
    # for i in range(3):
    #     # plot_every_single_UNcalib_chan(shot, dir = i, lim = 2, show = False, save = True)
    #     plot_every_single_UNcalib_chan(shot, dir = i, lim = 2, show = True, save = False)

    # plot_every_single_calib_chan(shot, dir = 2, lim = 50, save = False)


    # for i in range(3):
    #     plot_every_single_calib_chan(shot, dir = i, lim = 3000, save = True, show = False)




# only run main if not being imported
if __name__ == '__main__':
    main()
