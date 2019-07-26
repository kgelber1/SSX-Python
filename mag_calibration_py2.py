#!/usr/bin/env python
"""Calibration routines for the mag probes.
   This is the latest version, written in python 2 for compatibility
   with the rest of Tim's code.

   This is the modified version from DSchaffner combined with the
   edits from ADLIGHT.

   Also combined with SSX_bdot_calibration_routines.py

   New edits are all the way down at the bottom...
"""

# 5/13/08 9:51 PM by Tim Gray
# AECDD892-E7D5-4506-A7EE-428391D254FA

__author__ = "Tim Gray"
__version__ = "1.1"

import os
import sys

from pylab import *
from numpy import *
# from ssx import *
# import ssx_py_utils_david as ssxutil
import ssx_data_read_basic as sdr	# don't need?
import ssx_py_utils_basic as ssxutil
import mjtools as mj
import scipy as sp
from matplotlib.pylab import *

############# only need this? #################################################
def polyPeak(time, data, timerange = [40,80], range = (50,130), timebase = 10, pretrigger = 20):
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
    # t1, t2 = (array(range) + pretrigger) * timebase
    # # generate an array of indices spanning the range
    # ti = arange(t1,t2, dtype='int16')
    # ti = ti[:len(ti)-1]
    # # print(len(ti-))
    # print("ti is: %s"%(ti))
    t1 = mj.tindex(time,timerange[0])#+pretrigger)
    t2 = mj.tindex(time,timerange[1])#+pretrigger)
    # print 't1=', t1
    # print 't2=', t2
    # generate an array of indices spanning the range
    ti = arange(t1,t2)
    # get the time and data points in the range
    t = time[ti]
    d = data[ti]
    # Fit a 2nd degree polynomial and find the min and max.
    p = polyfit(t,d,2)
    fit = p[0]*t**2 + p[1]*t + p[2]
    dataMax = fit.max()
    dataMin = fit.min()
    if abs(dataMin) > dataMax:
        dataMax = dataMin
    return dataMax

def helmholtz2(r, i = 1.0, coil = 1):
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

    b = zeros(3)

    for j in xrange(360):	# compute line integral
        dth = 1./360 * 2 * pi
        th = j * dth

        rho = sqrt( (r1 - a * cos(th))**2 + (r2 - a * sin(th))**2 + (r3 - .5 *
            d)**2 )
        b[0] = b[0] + a * dth * (- r3 * cos(th))/rho**3
        b[1] = b[1] + a * dth * (- r3 * sin(th))/rho**3
        b[2] = b[2] + a * dth * ( (r1 - a * cos(th)) * cos(th) + (r2 - a *
            sin(th)) * sin(th) )/rho**3

        rho = sqrt( (r1 - a * cos(th))**2 + (r2 - a * sin(th))**2 + (r3 + .5 *
            d)**2 )
        b[0] = b[0] + a * dth * (- r3 * cos(th))/rho**3
        b[1] = b[1] + a * dth * (- r3 * sin(th))/rho**3
        b[2] = b[2] + a * dth * ( (r1 - a * cos(th)) * cos(th) + (r2 - a *
            sin(th)) * sin(th) )/rho**3


    # calculate the be field.  The 4pi's cancel out, and the 1e4 is to convert
    # to gauss
    b = b * i * 1.0e-7 * turns * 1e4
    return b


def calibBFields2(dir, calibInfo, position_vectors):
    """Calculates B field of Helmholtz coil for probe tips in
    direction dir.

    probeLocs is in mm.  Could be arbitrarily extended for other probe
    configurations as long as the right locations are entered.

    This code almost always assumes 1 amp of current in the helmholtz coil.
    This is compensated for at a later point in time."""
    probeLocs = calibInfo["probeLocs"]
    coil = calibInfo["coil"]
    numProbes = len(probeLocs)
    bb = zeros((numProbes,3))

    for i in xrange(numProbes):
        # ADL 20190501
        # if only one probe stalk is calibrated at a time, spacing information
        # should be provided in probeLocs as far as I understand

        if dir == 'r': # z direction of probe is along axis of Helmholtz coil
            r = position_vectors[i]
            b = helmholtz2(r, coil = coil)
            # r shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[2], 1*b[0], 1*b[1] ]
        elif dir == 't': # x direction of probe is along axis of Helmholtz coil
    #### NEED TO VERIFY WHETHER THIS IS CORRECT
            # looks like -90 degree rotation about y-axis of probe coil orientation, so switch x and z
            # print("r is: ",  position_vectors[i])
            r  = position_vectors[i]
            r = [r[2], r[1], r[0]]
            # r = position_vectors[i].transpose((2,1,0))
            b = helmholtz2(r, coil = coil)
            # t shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[0], 1*b[2], 1*b[1] ]
        elif dir == 'z': # y direction of probe is along axis of Helmholtz coil
            # same as t, but different transform below
            # ADL 20190501 (coils are in same positions, but field directions rotate
            # r = position_vectors[i].transpose((2,1,0))
            r  = position_vectors[i]
            r = [r[2], r[1], r[0]]
            b = helmholtz2(r, coil = coil)
            # z shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[0], 1*b[1], 1*b[2] ]
    return bb

def get_probeLocs_20190501():
    """
    Provides positions in meters along probe stalks for 4x4 array of probes built by M. Kaur.


    ADL -- 20190501:
    -----
    The way I understand the code to work, it was designed to calibrate a single probe stalk at a time,
    with the spacing between coil sets provided as calibInfo['probeLocs'].

    I think this can be made obsolete by specifying vector coil positions rather than 1D coil positions.
    To maintain backwards compatibility, I have introduced a new keyword to calibBFields2
    (which I understand to be the most current version of the Helmholtz field finder) to
    flag whether the position information is 1D or 3D.
    """
    x1D = np.array([-4.25, 4.24],float)*1e-3  #written down in mm, but convert immediately to meters
    y1D = np.array([-4.25, 4.25],float)*1e-3 #written down in mm, but convert immediately to meters
    z1D = np.array([-2.25, -0.75, 0.75, 2.25],float)*1e-3*25.4 #written down in inches, but convert immediately to meters

    x,y,z = np.meshgrid(x1D,y1D,z1D)

    vector_position = np.array(zip(x.flatten(),y.flatten(),z.flatten()),float)

    return vector_position

def calibrationInfo_4x4(runs, probenum = '1', scope = 2, currentChan = 2, numProbes
    = 4, dir = ['r', 't', 'z'], coil = 2, pretrigger = 100, probeLocs=None):
    """Sets up dictionary with all the pertinent information for calibration.

    runs should be the set of shots taken in the order of the dir list.
    - pretrigger is the number of microseconds before t=0 of data the mag data
      has recorded.
    - coil is the calibration coil used"""

    calibInfo = {}
    calibInfo['probeNum'] = probenum
    calibInfo['runs'] = runs
    calibInfo['scope'] = scope
    calibInfo['numProbes'] = numProbes
    calibInfo['currentChan'] = currentChan
    calibInfo['numChan'] = numProbes * len(dir)
    calibInfo['dir'] = dir
    calibInfo['coil'] = coil
    calibInfo['pretrigger'] = pretrigger
    # So by defualt probe locs is set to be None, even though then it is used in
    # Calib b-fields, so instead use the funciton adam wrote:
    calibInfo['probeLocs'] = get_probeLocs_20190501()


    return calibInfo
#
# def load_smallbdot(shot, scope = '1'):
#     data = sdr.scope_data(shot, scope)
#     time = data.time
#     X_bdot = data.ch1 #*170000#East gun 170kA/V
#     Y_bdot = data.ch2 #*190000#West gun 190kA/V
#     Z_bdot = data.ch3
#     Bx = sp.integrate.cumtrapz(X_bdot,time)
#     By = sp.integrate.cumtrapz(Y_bdot,time)
#     Bz = sp.integrate.cumtrapz(Z_bdot,time)
#     timeb = time[:-1]
#     # return timeb, Bx, By, Bz

def calibrationDTAC_052119(calibInfo, probeinfo, prefunc = None):
    """Modification to modification to calibration() for the dtac data.

    Data from the dtac units needs to be integrated before we calibrate it.
    The integrated data is used for fitting.

    probeinfo has some dtac specific data stashed in it.
    - probe is the probe name (m1).
    - filestrings is the dtac filestring names.
    - dir relates the probe coords (1,2,3) to lab coords (r,t,z).
      Specifically, you pass in ['r', 't', 'z'] in this variable (in the
      appropriate order) if the probe data has *different* names, like ['1',
      '2', '3'].  This probably will never need to be used for the dtac probes,
      but it's here in just in case.  This was originally for the colorado
      probe - look in the functions calibcolorado and calibrationCol for more
      on the implementation."""

    import basic_hiresmag as hr

    nump = calibInfo['numProbes']
    numc = calibInfo['numChan']
    pretrig = calibInfo['pretrigger']
    probeLocs = calibInfo['probeLocs']

    calibData =np.zeros((3, nump, 4))

    # dir = c,b,a
    # dir = a,b,c
    # dir2 = r,t,z
    dir = calibInfo['dir']
    # helmholtz fields - triplet location, shot orientation, components (r,t,z)
    b =np.zeros((nump,3,3))
    # data from 3 calibration shots, triplet position, components (r,z,t)
    voltages =np.zeros((3, nump, 3))
    dir2 = probeinfo['dir']

    for i, k in enumerate(dir):
        # print("i is:", i, "And dir is: ", dir)
        # Not sure if this fixes it but:

        d = dir2[i]
        print ("Calculating b fields - %s" % d)
        # print("calibInfo:", calibInfo)
        # calculate the field from the helmholtz coils as *should* be measured
        # by each probe in the right location
        b0 = calibBFields2(d, calibInfo, probeLocs)
        # print("is it here?")
        # stick it into the big b array
        # print("b.shape: ", b.shape) #4,3,3

        # TODO: should be 16 (locations),3(coils at each location),3(measurements)
        # print("b0.shape: ", b0.shape)#48,3
        # "stored in a 48x4 array, where the first three lines formed a 3x4 array
        # that consisted of 1r,1t,1z data, etc"
        # print(b0)

        b[:,i,:] = b0[,:]

            # IDK WHAT THIS SIZE THING IS:
            # b = b0

    if calibInfo['coil'] == 2:
        peakrange = [25,60]
    elif calibInfo['coil'] == 1:
        peakrange = [50,130]

    # Gets to this point just fine, but then it seems to not actually get the scope_data
    for i,run in enumerate(calibInfo['runs']):
        # print("i is:", i)
        # print("unp.sing run %s" % run)
        # Get the data from the scopes
        print(run, calibInfo['scope'])
        # Think instead of calibInfo[scope], I want prove info??

        shotdata = sdr.scope_data(run, calibInfo['scope'][:1])
        # get the calibrated current for the shot
        current = getattr(shotdata, shotdata.names[calibInfo['currentChan'] -
            1]) * 170000. #(1000./0.00541/1.025)
        # the above is our calibration for the gun current.
        # - 170 kA/V for east gun
        # - 190 kA/V for west gun
        # - 180 kA/V was the historical value

        # get the peak of the current
        currentMax = polyPeak(shotdata.time, current, range = peakrange, timebase
            = 200)
        # read in the dtac mag data and get it ready
        # hiresmag_data should be kosher to use for any dtac setup, but I'm not
        # sure.  We might have to modify this.
        # print("the probe infor thing is:", probeinfo['probe'])
        # magdata = hr.hiresmag_data(run, probeinfo['probe'])
        # magdata.filestrings = probeinfo['filestrings']
        # magdata._getData()
        # magdata._makeFullData()
        # magdata.removeOffsets()

        # or alternaively just pass the filestrings in
        fstrings = probeinfo['filestrings']
        # fstrings.reverse()
        # print(type(fstrings))
        magdata = hr.hiresmag_data(run, probeinfo['probe'], filestrings= fstrings)
        # call our preprocesnp.sing function on the data, if we need it
        if prefunc:
            prefunc(magdata)
        # print("magdata.fullData.shape", magdata.fullData.shape)
        # print("data is", magdata.fullData)
        # Create a new array of the same size, where every element is zero:
        magdata.iFullData = np.ma.zeros(magdata.fullData.shape)
        # for x in magdata.iFullData:
            # magdata.iFullData[x,:]
        # print("zeros/ shape:",  magdata.iFullData, magdata.iFullData.shape)
        magdata.iFullData.mask = True
        # print("Mask:", magdata.iFullData)
        magdata.iFullData[:,1:] = sp.integrate.cumtrapz(magdata.fullData, dx =
            magdata.deltat * 1e6, axis = 1)
        # print("Should have data:", magdata.iFullData)
        # print('fullData ',magdata.fullData.shape)
        # print('iFullData ',magdata.iFullData.shape)

        magmax =np.zeros((nump, 3))

        # Outer loop is for each probe location (i.e. 16 probes)
        # Inner loop is for directions (i.e. r,t,z)
        # So the progression goes: 1r, 1t, 1z, 2r, 2t, 2z...
        for j in range(nump):
            for k,p in enumerate(dir):
                # first figure out the name of our probe channel - ex. m2r1
                probeNum = "%s%s%i" % (probeinfo['probe'], p, j+1)
                # print( "%s: shot_index = %i, position_index = %i, direct_index = %i" % (probeNum, i,
                    # j, k))
                # find the peak in the mag signal
                # chan = magdata.channelNames.index(probeNum)
                chan = (16*k)+j #Dschaffner changed 5/16/2014 to convert [3,16] to [48]
                # chan = k
                # print('Chan = ',chan)
                timebase = int(1e-6 / magdata.deltat)
                # tmpmag = polyPeak(magdata.time, magdata.iFullData[chan],
                #     timebase = timebase, pretrigger = 100)
                # print("should be some data:", magdata.iFullData)
                # print("\n")
                print(peakrange)
                print(timebase,pretrig)
                # tmpmag = polyPeak(magdata.time, magdata.iFullData[chan,:], range
                #     = peakrange, timebase = timebase, pretrigger = pretrig)
                tmpmag = polyPeak(magdata.time, magdata.iFullData[:], range
                    = peakrange, timebase = timebase, pretrigger = pretrig)

                    #DSchaffner changed magdata.iFullData[chan] to magdata.iFullData[k,j]
                # print tmpmag, currentMax, tmpmag/ currentMax
                # normalize to the current
                # No longer normalize to the double turn helmholtz coil.  First
                # of all that was stunp.pid - why not just snp.pit that out in the
                # helmholtz function?  Secondly, we have a new coil which is
                # only one turn, so we need to be able to toggle that on and
                # off, and it seems to make the most sense to do it in the...
                # helmholtz function.
                #
                # We will still normalize the measured field here with the
                # current
                magmax[j, k] =  tmpmag / (currentMax)
        # our normalized B signals - takes into account current
        # magmax is now number of probes x number of axis long.  So for a 16
        # channel, 3 axis probe, it should be 48 channels (16,3).
        # Voltages is the same size as magmax, but with an extra dimension. The
        # extra dimension is number of calibration directions (3
        # usually).  So it would be [3,16,3].
        voltages[i] = magmax

    # Outer loop is for each probe location (i.e. 16 probes)
    # inner loops is r, t, z
    for i in range(nump):
        for j,d in enumerate(dir):
            # so we are going to dot the normalized voltages into the inverse
            # of the expected field.  The indices in voltages makes sure we get
            # a np.single coil (m2r1) for all three calibration orientations
            # (r,t,z).

            u = dot((voltages[:,i,j]), inv(b[i,:,:]))
            # all the magic is in the above line.  All we are doing is dotting
            # the measured voltages (r,t,z) into the inverse of the 3x3 array
            # of the fields we should be measuring ((r,t,z) for each of the 3
            # orientations r, t, and z (though I think the indices of b are as
            # follows:
            #   - first index is probe number
            #   - second index is calibration orientation
            #   - third index is component of field
            # dot sums over the last dimension of the first array and the 2nd
            # to last dimension of the second array.  So we are dotting the
            # calib orientation of voltages and the calib orientation of b.  At
            # this point in time, b should more or less be diagonal (ignoring
            # any error fields in the coil due to non-uniform fields in the
            # helmholtz).
            #
            # this normalizes it all
            modu = np.sqrt(dot(u,u))
            calibData[j,i,:3] = u/modu
            calibData[j,i,3] = modu
    # Reshape before we send it back.  This is getting written to a file so it
    # needs to be 2-D.  Note, I think the data is going to get put down in the
    # file as 1-16r, 1-16t, then 1-16z.  When we read it in later, just run
    # reshape() again, this time like reshape((3, nump, 4)) where nump is 16
    # most likely.
    #
    # This is a change of how this data used to be stored.  It used to be
    # stored in a 48x4 array, where the first three lines formed a 3x4 array
    # that consisted of 1r,1t,1z data, etc.  We also just used to parse the
    # array with a bunch of stunp.pid indices.  Now we reshape to get our 3x16 for
    # axis x probe.  Mind you, we choose 3x16 instead of 16x3 because of the
    # way the dtac data is read in, organized and easily reshaped.  This way we
    # don't have to think about indices (dtac data is 3x16) nor do we have to
    # worry about making the 3D arrays 2D, just use reshape and things go
    # together in the correct way for unreshanp.ping it later.
    return calibData.reshape((nump*3,4))

def calib052119():
    """Calibration from 5/21/19
    written by Katie

    4 4 channel probes in gun region."""
    #date the calibration data was taken
    date = '050119'


    pth = os.getcwd() + '/data/2019/050119/'
    # pth = 'magcalib'
    # pth = ssxutil.ssxPath(pth)
    # print(pth)
    hrinfo = {}
    hrinfo['filestrings'] = ['mag1', 'mag2', 'mag3']
    hrinfo['probe'] = 'm2'
    # not really needed but we'll leave it in.
    # hrinfo['dir'] = ci1['dir']

    # print(pth)
    # probe 1 r,t,z shots and same for probe 2-4
    p1runs, p2runs, p3runs, p4runs = [], [], [], []
    p1 = np.arange(11,16)#shots 11-15 are z -> r dir
    p2 = np.arange(17,21)#shots 17-20 are x -> t dir
    p3 = np.arange(21,26)#shots 21-25 are y -> z dir
    # p2 = ['4', '5', '6']
    # print(p1)
    # p3 = ['9', '8', '7']
    p4 = np.arange(11,26)

    for p in p1:
        p1runs.append(date+'r'+str(p))
    for p in p2:
        p2runs.append(date+'r'+str(p))
    for p in p3:
        p3runs.append(date+'r'+str(p))
    for p in p4:
        p4runs.append(date+'r'+str(p))
    # Calibration info should find this???
    # print("p1 runs, that is being sent into calibration info")
    # print(p1runs)
    ci1 = calibrationInfo_4x4(p1runs, dir = 'r', scope = '2', numProbes = 4)
    # ci2 = calibrationInfo_4x4(p2runs, dir = 't', scope = '2X', numProbes = 4)
    # ci3 = calibrationInfo_4x4(p3runs,dir = 'z', scope = '2X', numProbes = 4)
    ci4 = calibrationInfo_4x4(p4runs, scope = '2', numProbes = 4)

    hrinfo = {}
    hrinfo['filestrings'] = ['mag1', 'mag2', 'mag3']
    hrinfo['probe'] = 'm1'
    # not really needed but we'll leave it in.
    hrinfo['dir'] = ci1['dir']
    # ci4 is all the x, y and z directions
    probe1Calib = calibrationDTAC_052119(ci1, hrinfo)

    print("Finished! File saved as calib-%s-4x4_tims.txt in current working directory" %(date))
    savetxt(os.path.join(pth, 'calib%4x4_tims.txt' % date), probe1Calib)

def main():
    # calib052119()
    wrapper_to_run_davids()

# Only runs main when NOT imported (ie runs when called from terminal)
if __name__ == "__main__":
   main()


def wrapper_to_run_davids():
    # want every shot but 16
    date = '050119'
    shot_nums = np.arange(11,16)
    shot_nums = np.append(np.arrange(17,26))
    shots = []
    for shot in shot_nums:
        shots.append(date+'r'+str(p))


    # pth = os.getcwd() + '/data/2019/050119/'
    pth = os.getcwd()
    print("Finished! File saved as calib-%s-4x4_tims.txt in current working directory" %(date))
    savetxt(os.path.join(pth, 'calib%4x4_tims.txt' % date), probe1Calib)



def calibrationDTAC_david(calibInfo, probeinfo, prefunc = None):
    """Modification to calibration() for the dtac data.

    Data from the dtac units needs to be integrated before we calibrate it.
    The integrated data is used for fitting.

    probeinfo has some dtac specific data stashed in it.
    - probe is the probe name (m1).
    - filestrings is the dtac filestring names.
    - dir relates the probe coords (1,2,3) to lab coords (r,t,z).
      Specifically, you pass in ['r', 't', 'z'] in this variable (in the
      appropriate order) if the probe data has *different* names, like ['1',
      '2', '3'].  This probably will never need to be used for the dtac probes,
      but it's here in just in case.  This was originally for the colorado
      probe - look in the functions calibcolorado and calibrationCol for more
      on the implementation."""

    import hiresmag_david as hr

    nump = calibInfo['numProbes']
    numc = calibInfo['numChan']
    pretrig = calibInfo['pretrigger']

    calibData = zeros((3, nump, 4))

    # dir = c,b,a
    # dir = a,b,c
    # dir2 = r,t,z
    dir = calibInfo['dir']
    # helmholtz fields - triplet location, shot orientation, components (r,t,z)
    b = zeros((nump,3,3))
    # data from 3 calibration shots, triplet position, components (r,z,t)
    voltages = zeros((3, nump, 3))
    dir2 = probeinfo['dir']

    for i, k in enumerate(dir):
        d = dir2[i]
        print "Calculating b fields - %s" % d
        # calculate the field from the helmholtz coils as *should* be measured
        # by each probe in the right location
        b0 = calibBFields2(d, calibInfo)
        # stick it into the big b array
        b[:,i,:] = b0

    if calibInfo['coil'] == 2:
        peakrange = [25,60]
    elif calibInfo['coil'] == 1:
        peakrange = [50,130]

    for i,run in enumerate(calibInfo['runs']):
        print "using run %s" % run
        # Get the data from the scopes
        shotdata = sdr.scope_data(run, calibInfo['scope'])
        # get the calibrated current for the shot
        current = getattr(shotdata, shotdata.names[calibInfo['currentChan'] -
            1]) * 170000. #(1000./0.00541/1.025)
        # the above is our calibration for the gun current.
        # - 170 kA/V for east gun
        # - 190 kA/V for west gun
        # - 180 kA/V was the historical value

        # get the peak of the current
        currentMax = polyPeak(shotdata.time, current, range = peakrange, timebase
            = 200)

        # read in the dtac mag data and get it ready
        # hiresmag_data should be kosher to use for any dtac setup, but I'm not
        # sure.  We might have to modify this.
        magdata = hr.hiresmag_data(run, probeinfo['probe'])
        magdata.filestrings = probeinfo['filestrings']
        magdata._getData()
        magdata._makeFullData()
        magdata.removeOffsets()
        # call our preprocessing function on the data, if we need it
        if prefunc:
            prefunc(magdata)
        magdata.iFullData = np.ma.zeros(magdata.fullData.shape)
        magdata.iFullData.mask = True
        magdata.iFullData[:,1:] = sp.integrate.cumtrapz(magdata.fullData, dx =
            magdata.deltat * 1e6, axis = 1)
        print 'fullData ',magdata.fullData.shape
        print 'iFullData ',magdata.iFullData.shape

        # tmp = ma.masked_outside(magdata.fullData[:,2:100], -.50, .50).mean(1)
        # offsets = np.expand_dims(tmp, 2)
        # dat = magdata.fullData - offsets
        # magdata.iFullData = np.ma.zeros(magdata.fullData.shape)
        # magdata.iFullData.mask = True
        # magdata.iFullData[:,1:] = sp.integrate.cumtrapz(dat, dx =
        #     magdata.deltat * 1e6, axis = 1)

        # magmax = zeros(numc)
        magmax = zeros((nump, 3))

        # Outer loop is for each probe location (i.e. 16 probes)
        # Inner loop is for directions (i.e. r,t,z)
        # So the progression goes: 1r, 1t, 1z, 2r, 2t, 2z...
        for j in xrange(nump):
            for k,p in enumerate(dir):
                # first figure out the name of our probe channel - ex. m2r1
                probeNum = "%s%s%i" % (probeinfo['probe'], p, j+1)
                print "%s: i = %i, j = %i, k = %i" % (probeNum, i,
                    j, k)
                # find the peak in the mag signal
                #chan = magdata.channelNames.index(probeNum)
                chan = (16*k)+j #Dschaffner changed 5/16/2014 to convert [3,16] to [48]
                print 'Chan = ',chan
                timebase = int(1e-6 / magdata.deltat)
                # tmpmag = polyPeak(magdata.time, magdata.iFullData[chan],
                #     timebase = timebase, pretrigger = 100)
                tmpmag = polyPeak(magdata.time, magdata.iFullData[chan,:], range
                    = peakrange, timebase = timebase, pretrigger = pretrig)
                    #DSchaffner changed magdata.iFullData[chan] to magdata.iFullData[k,j]
                # print tmpmag, currentMax, tmpmag/ currentMax
                # normalize to the current
                # No longer normalize to the double turn helmholtz coil.  First
                # of all that was stupid - why not just spit that out in the
                # helmholtz function?  Secondly, we have a new coil which is
                # only one turn, so we need to be able to toggle that on and
                # off, and it seems to make the most sense to do it in the...
                # helmholtz function.
                #
                # We will still normalize the measured field here with the
                # current
                magmax[j, k] =  tmpmag / (currentMax)
        # our normalized B signals - takes into account current
        # magmax is now number of probes x number of axis long.  So for a 16
        # channel, 3 axis probe, it should be 48 channels (16,3).
        # Voltages is the same size as magmax, but with an extra dimension. The
        # extra dimension is number of calibration directions (3
        # usually).  So it would be [3,16,3].
        voltages[i] = magmax

    # Outer loop is for each probe location (i.e. 16 probes)
    # inner loops is r, t, z
    for i in xrange(nump):
        for j,d in enumerate(dir):
            # so we are going to dot the normalized voltages into the inverse
            # of the expected field.  The indices in voltages makes sure we get
            # a single coil (m2r1) for all three calibration orientations
            # (r,t,z).

            u = dot((voltages[:,i,j]), inv(b[i,:,:]))
            # all the magic is in the above line.  All we are doing is dotting
            # the measured voltages (r,t,z) into the inverse of the 3x3 array
            # of the fields we should be measuring ((r,t,z) for each of the 3
            # orientations r, t, and z (though I think the indices of b are as
            # follows:
            #   - first index is probe number
            #   - second index is calibration orientation
            #   - third index is component of field
            # dot sums over the last dimension of the first array and the 2nd
            # to last dimension of the second array.  So we are dotting the
            # calib orientation of voltages and the calib orientation of b.  At
            # this point in time, b should more or less be diagonal (ignoring
            # any error fields in the coil due to non-uniform fields in the
            # helmholtz).
            #
            # this normalizes it all
            modu = sqrt(dot(u,u))
            calibData[j,i,:3] = u/modu
            calibData[j,i,3] = modu
    # Reshape before we send it back.  This is getting written to a file so it
    # needs to be 2-D.  Note, I think the data is going to get put down in the
    # file as 1-16r, 1-16t, then 1-16z.  When we read it in later, just run
    # reshape() again, this time like reshape((3, nump, 4)) where nump is 16
    # most likely.
    #
    # This is a change of how this data used to be stored.  It used to be
    # stored in a 48x4 array, where the first three lines formed a 3x4 array
    # that consisted of 1r,1t,1z data, etc.  We also just used to parse the
    # array with a bunch of stupid indices.  Now we reshape to get our 3x16 for
    # axis x probe.  Mind you, we choose 3x16 instead of 16x3 because of the
    # way the dtac data is read in, organized and easily reshaped.  This way we
    # don't have to think about indices (dtac data is 3x16) nor do we have to
    # worry about making the 3D arrays 2D, just use reshape and things go
    # together in the correct way for unreshaping it later.
    return calibData.reshape((nump*3,4))










###################################
