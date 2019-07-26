#!/usr/bin/env python
"""Calibration routines for the mag probes."""

# 5/13/08 9:51 PM by Tim Gray
# AECDD892-E7D5-4506-A7EE-428391D254FA

__author__ = "Tim Gray"
__version__ = "1.0"

import os
import sys

#import matplotlib
#matplotlib.use('Agg')
#matplotlib.rc('axes', grid=False)

from pylab import *
from numpy import *
# from ssx import *
import ssx_data_read_david as sdr	# don't need?
import ssx_py_utils_david as ssxutil
import scipy as sp
from matplotlib.pylab import *

#ioff()

def helmholtz(r, i = 1.0):
    """Compute B field of Helmholtz coil at (x,y,z)"""
    r1, r2, r3 = r
    a = 6.1 * 0.0254	# radius of hoop in meters
    d = 5.9 * 0.0254	# distance between hoops (not quite helmholtz)

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
    
    b = b * i * 1.0e-3
    # now multiply by the current (in Amps), mu_0/4pi, and convert to Gauss
    return b

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

def calibBFields2(dir, calibInfo):
    """Calculates B field of Helmholtz coil for 8 (or 4) probe tips in
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
        # put the probelocations in meters too
        j = probeLocs[i]/1000.
        # took out the - signs from the original calibBFields.  We used to
        # have the HH field pointing in + r-hat, but in - z-hat and t-hat.  So
        # some minus signs were introduced.
        if dir == 'r':
            r = array([0, 0, j])
            b = helmholtz2(r, coil = coil)
            # r shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[2], 1*b[0], 1*b[1] ]	
        elif dir == 't':
            r = array([j, 0, 0])
            b = helmholtz2(r, coil = coil)
            # t shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[0], 1*b[2], 1*b[1] ]	
        elif dir == 'z':
            # same as t, but different transform below
            r = array([j, 0, 0])
            b = helmholtz2(r, coil = coil)
            # z shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[0], 1*b[1], 1*b[2] ]	
    return bb

def calibrationInfo(runs, probenum = '1', scope = 2, currentChan = 2, numProbes
    = 8, dir = ['r', 't', 'z'], coil = 1, pretrigger = 100):
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

    return calibInfo

def polyPeak(time, data, range = (50,130), timebase = 10, pretrigger = 20):
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
    t1, t2 = (array(range) + pretrigger) * timebase
    # generate an array of indices spanning the range
    ti = arange(t1,t2, dtype='int16')
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

def calibrationDTAC(calibInfo, probeinfo, prefunc = None):
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

def calibrationDTAC_1pt(calibInfo, probeinfo, probeindex,prefunc = None):
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
        peakrange = [25,70]
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
        plot(shotdata.time,current)

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
        #for the z-direction, set all fields to zero to simulate no data taken
        if i == 2:
            magdata.iFullData = 1e-9*np.ma.ones(magdata.iFullData.shape)
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
        
        #for k,p in enumerate(dir):
        if i == 0:
            j = probeindex[1]#vloc
            k = probeindex[0]#vdir
            print "%s: i = %i, j = %i, k = %i" % ('Calculating chan:',i,
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
            magmax[0, 0] =  tmpmag / (currentMax)
            
            j = probeindex[3]#vloc
            k = probeindex[2]#vdir
            print "%s: i = %i, j = %i, k = %i" % ('Calculating chan:',i,
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
            magmax[0, 1] =  tmpmag / (currentMax)
            
            print 'magmax ',magmax
            # our normalized B signals - takes into account current 
            # magmax is now number of probes x number of axis long.  So for a 16
            # channel, 3 axis probe, it should be 48 channels (16,3).
            # Voltages is the same size as magmax, but with an extra dimension. The
            # extra dimension is number of calibration directions (3
            # usually).  So it would be [3,16,3].
        if i == 1:
            j = probeindex[1]#hloc
            k = probeindex[0]#hdir
            print "%s: i = %i, j = %i, k = %i" % ('Calculating chan:',i,
                    j, k)
            chan = (16*k)+j
            print 'Chan = ',chan
            timebase = int(1e-6 / magdata.deltat)
            tmpmag = polyPeak(magdata.time, magdata.iFullData[chan,:], range
                    = peakrange, timebase = timebase, pretrigger = pretrig)
            magmax[0, 0] =  tmpmag / (currentMax)
            
            j = probeindex[3]#hloc
            k = probeindex[2]#hdir
            print "%s: i = %i, j = %i, k = %i" % ('Calculating chan:',i,
                    j, k)
            chan = (16*k)+j
            print 'Chan = ',chan
            timebase = int(1e-6 / magdata.deltat)
            tmpmag = polyPeak(magdata.time, magdata.iFullData[chan,:], range
                    = peakrange, timebase = timebase, pretrigger = pretrig)
            magmax[0, 1] =  tmpmag / (currentMax)

        if i == 2: #Doesn't matter, everything will be zero
            #j = probeindex[1]#vloc
            #k = probeindex[0]#vdir
            #print "%s: i = %i, j = %i, k = %i" % ('Calculating chan:',i,
            #        j, k)
            #chan = (16*k)+j
            #print 'Chan = ',chan
            #timebase = int(1e-6 / magdata.deltat)
            #tmpmag = polyPeak(magdata.time, magdata.iFullData[chan,:], range
            #        = peakrange, timebase = timebase, pretrigger = pretrig)
            magmax[0, 0] = 0#tmpmag / (currentMax)
            magmax[0, 1] = 0
            magmax[0, 2] = 1e-12
        voltages[i] = magmax

    # Outer loop is for each probe location (i.e. 16 probes)
    # inner loops is r, t, z
    print 'voltages :',voltages
    for i in xrange(nump):
        for j,d in enumerate(dir):
            # so we are going to dot the normalized voltages into the inverse
            # of the expected field.  The indices in voltages makes sure we get
            # a single coil (m2r1) for all three calibration orientations
            # (r,t,z).  

            u = dot((voltages[:,i,j]), inv(b[i,:,:]))
            print 'u: ',u
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
            print 'modu ',modu
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
    
def calibHires():
    """Calibration for the hi res probe - June 2011.
    
    This is the actual calibration routine that you call.  It specifies the
    correct shots and what orientation they are.  It also saves the calibration
    file to the calibration directory."""

    date = "061411"
    # find the calib fils in the default ssx location
    pth = 'magcalib'
    pth = ssxutil.ssxPath(pth)

    # probe r,t,z shots
    # for hires probe, it is labelled r, t, z
    # but I think there is some serious crosstalk between t and z
    # the probe was turned a good amount I suspect
    #
    # p1 must be in r, t, z ordering
    # p =  r     t     z
    p1 = ['14', '10', '12']
    p1runs = [date+'r'+p for p in p1]

    # set up some calibration info.  Mainly what the runs were, what the scope
    # is, and the order of the shots.  This order must be the same as the order
    # in p1runs.  Actually, I'm not sure of this - I think it has to be in
    # r,t,z order.  It works that way - don't mess with it.  Just make sure the
    # runs are in r,t,z order as well.
    ci1 = calibrationInfo(p1runs, scope = '3', probenum = 1, dir = ['r', 't',
        'z'], numProbes = 16, coil = 1)

    # Some extra info for the dtac style probes.  This is necessary so the
    # calibrationDTAC routine can load in the data.
    hrinfo = {}
    hrinfo['filestrings'] = ['mag1', 'mag2']
    hrinfo['probe'] = 'm1'
    # not really needed but we'll leave it in.
    hrinfo['dir'] = ci1['dir']

    # probe locs around 0.  2.9" long (0.180" * 16) or 73.7mm.  First probe
    # edge is 0.050" from tip, or is 0.14" from core tip to center of first
    # probe.  Glass adds another 2mm but this can be ignored for calibration.
    ploc = (arange(16) * .180 - .180*15/2) * 25.4
    ci1["probeLocs"] = ploc
    # make the data and write it out
    probe1Calib = calibrationDTAC(ci1, hrinfo)
    savetxt(os.path.join(pth, 'calib-%s-hires1.txt' % date) , probe1Calib)

def hires2pre(hiresdat):
    """Pre processing function for hires probe data.

    This removes some of the initial Bdot spike that is throwing off the
    integration of the single loop hires probe that Ken built."""

    offsets = hiresdat.fullData[:,1:200].mean(1)
    newdat = hiresdat.fullData.transpose() - offsets
    hiresdat.fullData = newdat.transpose()
    hiresdat.fullData[:,0:218] = 0
    oo = sp.integrate.cumtrapz(hiresdat.fullData, dx = hiresdat.deltat * 1.e6,
        axis=1)
    fig = figure()
    fig.clear()
    for o in xrange(oo.shape[0]):
        plot(hiresdat.time[1:], oo[o])

    # fig = figure()
    # fig.clear()
    # for o in xrange(hiresdat.fullData.shape[0]):
    #     plot(hiresdat.time, hiresdat.fullData[o])

def calibFlex(probenum = 1):
    """Calibration for the flex probes - 2012-07-19, with the new helmholtz
    coil.
    
    This is the actual calibration routine that you call.  It specifies the
    correct shots and what orientation they are.  It also saves the calibration
    file to the calibration directory."""

    # find the calib fils in the default ssx location
    pth = 'magcalib'
    pth = ssxutil.ssxPath(pth)

    # probe r,t,z shots
    # for hires probe, it is labelled r, t, z
    # but I think there is some serious crosstalk between t and z
    # the probe was turned a good amount I suspect
    #
    # p1 must be in r, t, z ordering
    # p =  r     t     z
    # The reason this must be in the correct order comes much later in the
    # math, when we dot the recorded voltages into the expected mag fields.
    # The expected mag fields components are ALWAYS in r,t,z order, so this
    # needs to be in the same order.

    # select the correct runs for the right probe
    # date = "071912"
    # if probenum == 1:
    #     #      r     t    z
    #     p1 = ['10', '4', '7']
    #     filestrs = ['mag1', 'mag2']
    #     pname = 'f1'
    # elif probenum == 2:
    #     p1 = ['14', '16', '19']
    #     filestrs = ['mag1', 'mag2']
    #     pname = 'f1'
    # elif probenum == 3:
    #     p1 = ['28', '22', '25']
    #     filestrs = ['mag2', 'mag3']
    #     pname = 'f2'
    # elif probenum == 4:
    #     p1 = ['31', '35', '37']
    #     filestrs = ['mag2', 'mag3']
    #     pname = 'f2'
    date = "072512"
    if probenum == 1:
        #      r     t    z
        p1 = ['46', '50', '51']
        filestrs = ['mag1', 'mag2']
        pname = 'f1'
    elif probenum == 2:
        p1 = ['25', '22', '19']
        filestrs = ['mag1', 'mag2']
        pname = 'f1'
    elif probenum == 3:
        p1 = ['43', '40', '37']
        filestrs = ['mag2', 'mag3']
        pname = 'f2'
    elif probenum == 4:
        date = "072612"
        p1 = ['8', '5', '2']
        filestrs = ['mag2', 'mag3']
        pname = 'f2'

    p1runs = [date+'r'+p for p in p1]
    print p1runs

    # set up some calibration info.  Mainly what the runs were, what the scope
    # is, and the order of the shots.  
    # I think the order needs to be the order of the indices in the array of
    # fullData, as they will be arranged by _makeFullData().  So a,b,c is
    # correct.  r,t,z is correct if the data actually ends up in that order.
    # pretrigger is the number of microseconds before t=0 of data the mag data
    # has recorded.
    ci1 = calibrationInfo(p1runs, scope = '3', probenum = probenum, dir = ['a',
        'b', 'c'], numProbes = 8, coil=2, pretrigger = 2)

    # Some extra info for the dtac style probes.  This is necessary so the
    # calibrationDTAC routine can load in the data.
    hrinfo = {}
    hrinfo['filestrings'] = filestrs
    # pname = 'f' + str(probenum)
    hrinfo['probe'] = pname
    # This is the line that relates our probe axis a,b,c to the calib coil
    # coordinates.  This is absolutely needed if the probe channels are not
    # r,t,z.  The order of this is the same order as the 'dir' argument to
    # calibrationInfo.  In the case of the flex probe, a -> z, b -> t, c -> r,
    # so it should be z,t,r.
    hrinfo['dir'] = ['z', 't', 'r']

    # probe locs around 0.  3" long (0.375" * 8).  Probe spacing is 3/8".
    ploc = (arange(8) * 3/8. - 3/8.*7/2) * 25.4
    ci1["probeLocs"] = ploc
    # make the data and write it out
    # probe1Calib = calibrationDTAC(ci1, hrinfo, prefunc = hires2pre)
    # don't need the prefunc anymore for this probe.
    probe1Calib = calibrationDTAC(ci1, hrinfo)
    # return probe1Calib
    savetxt(os.path.join(pth, 'calib-%s-flex%i.txt' % (date, probenum)) ,
        probe1Calib)
        
def calibHiresLongProbe(channel_number):
    """Calibration for the long probe, with the new helmholtz coil - 5-27-2014.
    
    This is the actual calibration routine that you call.  It specifies the
    correct shots and what orientation they are.  It also saves the calibration
    file to the calibration directory."""

    # find the calib fils in the default ssx location
    pth = 'magcalib'
    pth = ssxutil.ssxPath(pth)

    # probe r,t,z shots
    # for hires probe, it is labelled r, t, z
    # but I think there is some serious crosstalk between t and z
    # the probe was turned a good amount I suspect
    #
    # p1 must be in r, t, z ordering
    # p =  r     t     z

    # Try 1
    # date = "110311"
    # p1 = ['8', '13', '16']
    
    # Try 2
    # date = "121211"
    # p1 = ['9', '5', '2']
    
    # Try 3 - everything is good I hope
    #date = "121611"#(last data run by Tim
    #p1 = ['17', '15', '12']
    #date = "051514"
    #p1 = ['14','5','9']
    #channel_number = 2
    date = "052714"
    #probeindex[vdir,vloc,hdir,hloc]
    vdir = 0 #long probe direction reference vertical (r=0,t=1,z=2)
    hdir = 0 #long probe direction reference vertical (r=0,t=1,z=2)
    vloc = 0 #long probe location vertical
    hloc = 8 #long probe location horizontal
    #example: First probe tip is at v->r1 and h->r9 or v->(0,0) and h->(0,8)
    if channel_number == 1: 
        p1 = ['36','18','1']#shots for channel 1
        probeindex = [0,0,0,8]#indices for channel 1
    if channel_number == 2: 
        p1 = ['35','17','1']#shots for channel 2
        probeindex = [0,1,0,9]#indices for channel 2
    if channel_number == 3: 
        p1 = ['34','16','1']#shots for channel 3
        probeindex = [0,2,0,10]#indices for channel 3
    if channel_number == 4: 
        p1 = ['33','15','1']#shots for channel 4
        probeindex = [0,3,0,11]#indices for channel 4
    if channel_number == 5: 
        p1 = ['32','14','1']#shots for channel 5
        probeindex = [0,4,0,12]#indices for channel 5
    if channel_number == 6: 
        p1 = ['31','14','1']#shots for channel 6 (channel H6 was missed, using shot for H5)
        probeindex = [0,5,0,13]#indices for channel 6
    if channel_number == 7: 
        p1 = ['30','13','1']#shots for channel 7
        probeindex = [0,6,0,14]#indices for channel 7
    if channel_number == 8:
        p1 = ['29','12','1']#shots for channel 8
        probeindex = [0,7,0,15]#indices for channel 8
    if channel_number == 9:
        p1 = ['28','11','1']#shots for channel 9
        probeindex = [1,0,1,8]#indices for channel 9
    if channel_number == 10:
        p1 = ['27','10','1']#shots for channel 10
        probeindex = [1,1,1,10]#indices for channel 10
    if channel_number == 11:
        p1 = ['26','9','1']#shots for channel 11
        probeindex = [1,2,1,11]#indices for channel 11
    if channel_number == 12:
        p1 = ['25','8','1']#shots for channel 12
        probeindex = [1,3,1,14]#indices for channel 12
    if channel_number == 13:
        p1 = ['24','7','1']#shots for channel 13
        probeindex = [1,4,1,15]#indices for channel 13
    if channel_number == 14:
        p1 = ['23','6','1']#shots for channel 14
        probeindex = [1,5,2,9]#indices for channel 14
    if channel_number == 15:
        p1 = ['22','5','1']#shots for channel 15
        probeindex = [1,6,2,10]#indices for channel 15
    if channel_number == 16:
        p1 = ['21','4','1']#shots for channel 16
        probeindex = [1,7,2,11]#indices for channel 16
    if channel_number == 17:
        p1 = ['20','3','1']#shots for channel 17
        probeindex = [2,0,2,12]#indices for channel 17
    if channel_number == 18:
        p1 = ['19','2','1']#shots for channel 18
        probeindex = [2,1,2,13]#indices for channel 18
    p1runs = [date+'r'+p for p in p1]
    channel_label = str(channel_number)
    print p1runs

    # set up some calibration info.  Mainly what the runs were, what the scope
    # is, and the order of the shots.  This order must be the same as the order
    # in p1runs.  Actually, I'm not sure of this - I think it has to be in
    # r,t,z order.  It works that way - don't mess with it.  Just make sure the
    # runs are in r,t,z order as well.
    # pretrigger is the number of microseconds before t=0 of data the mag data
    # has recorded
    ci1 = calibrationInfo(p1runs, scope = '3', probenum = 1, dir = ['r','t','z'], 
        numProbes = 1, coil=2, pretrigger = 2)

    # Some extra info for the dtac style probes.  This is necessary so the
    # calibrationDTAC routine can load in the data.
    hrinfo = {}
    hrinfo['filestrings'] = ['mag1', 'mag2', 'mag3']
    hrinfo['probe'] = 'm2'
    # not really needed but we'll leave it in.
    hrinfo['dir'] = ci1['dir']

    #each probe tip of the long probe was set to be in the center of the coil
    ploc = arange(1)
    ci1["probeLocs"] = ploc
    # make the data and write it out
    # probe1Calib = calibrationDTAC(ci1, hrinfo, prefunc = hires2pre)
    # don't need the prefunc anymore for this probe.
    probe1Calib = calibrationDTAC_1pt(ci1, hrinfo, probeindex)
    # return probe1Calib
    print 'Path = ',pth
    savetxt(os.path.join(pth, 'calib-%s-longprobeDTAC_chan%s.txt' % (date, channel_label)) , probe1Calib)

def calibHiresShortSingle3mm():
    """Calibration for the short single triplet probe with 3mm loops, with the new helmholtz coil - 4-15-2015.
    
    This is the actual calibration routine that you call.  It specifies the
    correct shots and what orientation they are.  It also saves the calibration
    file to the calibration directory."""

    # find the calib fils in the default ssx location
    pth = 'magcalib'
    pth = ssxutil.ssxPath(pth)

    # probe r,t,z shots
    # for hires probe, it is labelled r, t, z
    # but I think there is some serious crosstalk between t and z
    # the probe was turned a good amount I suspect
    #
    # p1 must be in r, t, z ordering
    # p =  r     t     z

    # Try 1
    # date = "110311"
    # p1 = ['8', '13', '16']
    
    # Try 2
    # date = "121211"
    # p1 = ['9', '5', '2']
    
    # Try 3 - everything is good I hope
    #date = "121611"#(last data run by Tim
    #p1 = ['17', '15', '12']
    #date = "051514"
    #p1 = ['14','5','9']
    #channel_number = 2
    date = "041715"
    p1 = ['1','4','7']#shots for channel
    p1runs = [date+'r'+p for p in p1]
    print p1runs

    # set up some calibration info.  Mainly what the runs were, what the scope
    # is, and the order of the shots.  This order must be the same as the order
    # in p1runs.  Actually, I'm not sure of this - I think it has to be in
    # r,t,z order.  It works that way - don't mess with it.  Just make sure the
    # runs are in r,t,z order as well.
    # pretrigger is the number of microseconds before t=0 of data the mag data
    # has recorded
    ci1 = calibrationInfo(p1runs, scope = '3', probenum = 1, dir = ['r','t','z'], 
        numProbes = 1, coil=2, pretrigger = 2)

    # Some extra info for the dtac style probes.  This is necessary so the
    # calibrationDTAC routine can load in the data.
    hrinfo = {}
    hrinfo['filestrings'] = ['mag1', 'mag2', 'mag3']
    hrinfo['probe'] = 'm2'
    # not really needed but we'll leave it in.
    hrinfo['dir'] = ci1['dir']

    #each probe tip of the long probe was set to be in the center of the coil
    ploc = arange(1)
    ci1["probeLocs"] = ploc
    # make the data and write it out
    # probe1Calib = calibrationDTAC(ci1, hrinfo, prefunc = hires2pre)
    # don't need the prefunc anymore for this probe.
    probe1Calib = calibrationDTAC(ci1, hrinfo)
    # return probe1Calib
    print 'Path = ',pth
    savetxt(os.path.join(pth, 'calib-%s-shortsingle3mmDTAC.txt' % date) , probe1Calib)
     
    
def calibHiresKen2():
    """Calibration for the hi res probe #2, the single turn one - November
    2011, with the new helmholtz coil.
    
    This is the actual calibration routine that you call.  It specifies the
    correct shots and what orientation they are.  It also saves the calibration
    file to the calibration directory."""

    # find the calib fils in the default ssx location
    pth = 'magcalib'
    pth = ssxutil.ssxPath(pth)

    # probe r,t,z shots
    # for hires probe, it is labelled r, t, z
    # but I think there is some serious crosstalk between t and z
    # the probe was turned a good amount I suspect
    #
    # p1 must be in r, t, z ordering
    # p =  r     t     z

    # Try 1
    # date = "110311"
    # p1 = ['8', '13', '16']
    
    # Try 2
    # date = "121211"
    # p1 = ['9', '5', '2']
    
    # Try 3 - everything is good I hope
    #date = "121611"#(last data run by Tim
    #p1 = ['17', '15', '12']
    #date = "051514"
    #p1 = ['14','5','9']
    date = "051914"
    p1 = ['11','5','8']
    p1runs = [date+'r'+p for p in p1]
    print p1runs

    # set up some calibration info.  Mainly what the runs were, what the scope
    # is, and the order of the shots.  This order must be the same as the order
    # in p1runs.  Actually, I'm not sure of this - I think it has to be in
    # r,t,z order.  It works that way - don't mess with it.  Just make sure the
    # runs are in r,t,z order as well.
    # pretrigger is the number of microseconds before t=0 of data the mag data
    # has recorded
    ci1 = calibrationInfo(p1runs, scope = '3', probenum = 1, dir = ['r', 't',
        'z'], numProbes = 16, coil=2, pretrigger = 2)

    # Some extra info for the dtac style probes.  This is necessary so the
    # calibrationDTAC routine can load in the data.
    hrinfo = {}
    hrinfo['filestrings'] = ['mag1', 'mag2', 'mag3']
    hrinfo['probe'] = 'm2'
    # not really needed but we'll leave it in.
    hrinfo['dir'] = ci1['dir']

    # probe locs around 0.  2.9" long (0.180" * 16) or 73.7mm.  First probe
    # edge is 0.050" from tip, or is 0.14" from core tip to center of first
    # probe.  Glass adds another 2mm but this can be ignored for calibration.
    ploc = ((arange(16) * .180 - .180*15/2) * 25.4)+20.35
    #Dschaffner 5/16/2014 added 20.35 to ploc to account for position of probe for calibration run on 5/15/14 to 5/19/14   
    #this makes the first probe at -13.94mm from the Helmholtz coil center which is the orientation of the probe based on the calibration pic (100-4381)
    ci1["probeLocs"] = ploc
    # make the data and write it out
    # probe1Calib = calibrationDTAC(ci1, hrinfo, prefunc = hires2pre)
    # don't need the prefunc anymore for this probe.
    probe1Calib = calibrationDTAC(ci1, hrinfo)
    # return probe1Calib
    print 'Path = ',pth
    savetxt(os.path.join(pth, 'calib-%s-hires2.txt' % date) , probe1Calib)


def calibHiresKen():
    """Calibration for the hi res probe #2, the single turn one - July 2011.
    
    This is the actual calibration routine that you call.  It specifies the
    correct shots and what orientation they are.  It also saves the calibration
    file to the calibration directory."""

    date = "072511"
    # find the calib fils in the default ssx location
    pth = 'magcalib'
    pth = ssxutil.ssxPath(pth)

    # probe r,t,z shots
    # for hires probe, it is labelled r, t, z
    # but I think there is some serious crosstalk between t and z
    # the probe was turned a good amount I suspect
    #
    # p1 must be in r, t, z ordering
    # p =  r     t     z
    p1 = ['20', '12', '16']
    p1runs = [date+'r'+p for p in p1]

    # set up some calibration info.  Mainly what the runs were, what the scope
    # is, and the order of the shots.  This order must be the same as the order
    # in p1runs.  Actually, I'm not sure of this - I think it has to be in
    # r,t,z order.  It works that way - don't mess with it.  Just make sure the
    # runs are in r,t,z order as well.
    ci1 = calibrationInfo(p1runs, scope = '3', probenum = 1, dir = ['r', 't',
        'z'], numProbes = 16, coil=1)

    # Some extra info for the dtac style probes.  This is necessary so the
    # calibrationDTAC routine can load in the data.
    hrinfo = {}
    hrinfo['filestrings'] = ['mag2', 'mag3']
    hrinfo['probe'] = 'm2'
    # not really needed but we'll leave it in.
    hrinfo['dir'] = ci1['dir']

    # probe locs around 0.  2.9" long (0.180" * 16) or 73.7mm.  First probe
    # edge is 0.050" from tip, or is 0.14" from core tip to center of first
    # probe.  Glass adds another 2mm but this can be ignored for calibration.
    ploc = ((arange(16) * .180 - .180*15/2) * 25.4)
    ci1["probeLocs"] = ploc
    # make the data and write it out
    probe1Calib = calibrationDTAC(ci1, hrinfo, prefunc = hires2pre)
    savetxt(os.path.join(pth, 'calib-%s-hires2.txt' % date) , probe1Calib)



# Everything below this is older and not as well documented.  It should all
# work but the routines above are the most recent and I tried to comment them
# better.

def calibBFields(dir, numProbes = 8, offset = -3.5):
    """Calculates B field of Helmholtz coil for 8 (or 4) probe tips in
    direction dir.

    Could be arbitrarily extended for other probe configurations as long as
    the right locations are entered."""
    bb = zeros((numProbes,3))

    for i in xrange(numProbes):
        if numProbes == 8:
            j = i*1.0 + offset
        elif numProbes == 4:
# 			j = 0.3125+i*1.0	#  commented out this line for hte 081308 calibrations
            j = i*1.0 + offset
        # We used to have the HH field pointing in + r-hat, but in - z-hat and
        # t-hat. So some minus signs were introduced.  This can be seen in the
        # one - for the t shots and the two - on the z shots below.
        if dir == 'r':
            r = array([0, 0, j]) * 0.0254
            b = helmholtz(r)
            # r shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[2], 1*b[0], 1*b[1] ]	
        elif dir == 't':
            r = array([j, 0, 0]) * 0.0254
            b = helmholtz(r)
            # t shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[0], -1*b[2], 1*b[1] ]	
        elif dir == 'z':
            r = array([j, 0, 0]) * 0.0254	#same as t, but different transform below
            b = helmholtz(r)
            # z shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[0], -1*b[1], -1*b[2] ]	
    return bb

def calibration(calibInfo):

    nump = calibInfo['numProbes']
    numc = calibInfo['numChan']

    calibData = zeros((numc, 4))

    dir = calibInfo['dir']

    b = zeros((nump,3,3))			# helmholtz fields - triplet location, shot orientation, components (r,t,z)
    voltages = zeros((3,numc))	# data from 3 calibration shots, 24 channels each (r,z,t orientation)


    for i,d in enumerate(dir):
        b0 = calibBFields(d, numProbes = nump) # calculate the field from the helmholtz coils as *should* be measured by each probe in the right location
        b[:,i,:] = b0	# stick it into the big b array

    for i,run in enumerate(calibInfo['runs']):
        shotdata = sdr.scope_data(run, calibInfo['scope'])
        current = getattr(shotdata, shotdata.names[calibInfo['currentChan'] - 1]) * (1000./0.00541/1.025)	# get the calibrated current for the shot

        currentMax = polyPeak(shotdata.time, current, timebase = 100)	# get the peak of the current

        magdata = sdr.magnetics_data(run)
        magmax = zeros(numc)

        for j in xrange(nump):
            for k,p in enumerate(dir):
# 				probeNum = magdata.channelGroups[str(j+1)][k]
                probeNum = "m%i%s%i" % (int(calibInfo['probeNum']), p, j+1)
                print probeNum
                tmpmag = polyPeak(magdata.time, getattr(magdata, probeNum)) # find the peak in the mag signal
                magmax[3*j + k] =  tmpmag / (2. * currentMax)	#normalize to the current and double B from the helmholtz because there are 2 turns per coil
        voltages[i] = magmax	# our normalized B signals - takes into account current and two turns of the helmholtz

    for i in xrange(nump):
        for j,d in enumerate(dir):
            u = dot(voltages[:,3*i+j], inv(b[i,:,:]))
            # all the magic is in the above line.  All we are doing is dotting
            # the measured voltages (r,t,z) into the inverse of the 3x3 array
            # of the fields we should be measuring ((r,t,z) for each of the 3
            # orientations r, t, and z
            modu = sqrt(dot(u,u))	# this normalizes it all
            calibData[3*i+j,:3] = u/modu
            calibData[3*i+j,3] = modu
    return calibData

def calibrationCol(calibInfo, dir2 = ['z', 'r', 't']):
    """Modification to calibration() for the colorado probe.

    dir2 relates the probe coords (1,2,3) to lab coords (r,t,z)."""

    nump = calibInfo['numProbes']
    numc = calibInfo['numChan']

    calibData = zeros((numc, 4))

    dir = calibInfo['dir']

    b = zeros((nump,3,3))			# helmholtz fields - triplet location, shot orientation, components (r,t,z)
    voltages = zeros((3,numc))	# data from 3 calibration shots, 24 channels each (r,z,t orientation)

    for i, k in enumerate(dir):
        d = dir2[i]
        print d
        b0 = calibBFields2(d, calibInfo) # calculate the field from the helmholtz coils as *should* be measured by each probe in the right location
        b[:,i,:] = b0	# stick it into the big b array

    for i,run in enumerate(calibInfo['runs']):
        shotdata = sdr.scope_data(run, calibInfo['scope'])
        current = getattr(shotdata, shotdata.names[calibInfo['currentChan'] - 1]) * (1000./0.00541/1.025)	# get the calibrated current for the shot

        currentMax = polyPeak(shotdata.time, current, timebase = 100)	# get the peak of the current

        magdata = sdr.magnetics_data(run)
        magmax = zeros(numc)

        for j in xrange(nump):
            for k,p in enumerate(dir):
# 				probeNum = magdata.channelGroups[str(j+1)][k]
                probeNum = "mc%s_%i" % (p, j+1)
                print probeNum
                tmpmag = polyPeak(magdata.time, getattr(magdata, probeNum)) # find the peak in the mag signal
                magmax[3*j + k] =  tmpmag / (2. * currentMax)	#normalize to the current and double B from the helmholtz because there are 2 turns per coil
        voltages[i] = magmax	# our normalized B signals - takes into account current and two turns of the helmholtz

    for i in xrange(nump):
        for j,d in enumerate(dir):
            u = dot(voltages[:,3*i+j], inv(b[i,:,:]))
            # all the magic is in the above line.  All we are doing is dotting
            # the measured voltages (r,t,z) into the inverse of the 3x3 array
            # of the fields we should be measuring ((r,t,z) for each of the 3
            # orientations r, t, and z
            modu = sqrt(dot(u,u))	# this normalizes it all
            calibData[3*i+j,:3] = u/modu
            calibData[3*i+j,3] = modu
    return calibData

def calib050908():
    """Calibration from 05/09/08.

    Two full 8 channel probes.  
    Probe 1 used cable Y.
    Probe 2 used cable U.

    Saves calibration data using numpy.  Use loadtxt('file') to retrieve."""

    date = '050908'
    pth = '/Users/tgray/Documents/SSX/Python/magcalib'
    # probe 1 r,t,z shots and same for probe 2
    p1 = ['3','1','2']
    p2 = ['6','4','5']

    p1runs, p2runs = [], []
    for p in p1:
        p1runs.append(date+'r'+p)
    for p in p2:
        p2runs.append(date+'r'+p)
    ci1 = calibrationInfo(p1runs)
    ci2 = calibrationInfo(p2runs)

    probe1Calib = calibration(ci1)
    savetxt(os.path.join(pth, 'calib050908p1.txt'), probe1Calib)
    probe2Calib = calibration(ci2)
    savetxt(os.path.join(pth, 'calib050908p2.txt'), probe2Calib)

def calib061208():
    """Calibration for quick quartz probes.

    Probes have 8 z channels and 4 theta."""

    date = '061308'
    pth = '/Users/tgray/Documents/SSX/Python/magcalib'
    # probe 1 r,t,z shots and same for probe 2,3
    p1 = ['2','3']
    p2 = ['4','5']
    p3 = []

    probe1Loc = {}
    probe2Loc = {}
    offset1 = (.434/.0254 - 3.72 - 10.25) * .0254 #434mm is length from bracket to first probe, - 3.72" for the distance between bracket and inner conflat edge while calibrating, - 10.25, which is the distance that the teflon standoff is from the center of the helmholtz coil  -  in mm
    offset2 = (.434/.0254 - 3.32 - 10.25) * .0254 #same as above, but for some reason I set it at 3.32"
    probe1Loc['z'] = (array([8, 23, 38, 53, 68, 83, 121, 160])- 8 )/ 1000. - offset1
    probe1Loc['t'] = (array([8, 38, 68, 121])- 8 )/ 1000. - offset1
    probe2Loc['z'] = (array([8, 23, 38, 53, 68, 83, 121, 160])- 8 )/ 1000. - offset2
    probe2Loc['t'] = (array([8, 38, 68, 121])- 8 )/ 1000. - offset2

    p1runs, p2runs, p3runs = [], [], []
    for p in p1:
        p1runs.append(date+'r'+p)
    for p in p2:
        p2runs.append(date+'r'+p)
# 	for p in p3:
# 		p3runs.append(date+'r'+p)
    ci1 = calibrationInfo(p1runs, numProbes = 8, dir = ['t', 'z'])
    ci1['numChan'] = 12
    ci2 = calibrationInfo(p2runs, numProbes = 8, dir = ['t', 'z'])
    ci2['numChan'] = 12
# 	ci2 = calibrationInfo(p2runs)
# 	ci3 = calibrationInfo(p3runs)
    print p1
    probe1Calib, probe1Calib2 = sloppyCalibration(ci1, probe1Loc)
    probe2Calib, probe2Calib2 = sloppyCalibration(ci2, probe2Loc)
    # merge by channel dicts
    probe1Calib.update(probe2Calib)
    probeCalib = probe1Calib
    print "probe 1 and 2 by channel"
    print probeCalib.__repr__()

    print "probe 1 t, then z"
    print probe1Calib2['t'].__repr__()
    print probe1Calib2['z'].__repr__()
    print "probe 2 t, then z"

    print probe2Calib2['t'].__repr__()
    print probe2Calib2['z'].__repr__()

def calib061708():
    """Calibration for quick quartz probes.  Take 2.

    Probes have 6 r channels and 6 theta or z."""

    date = '061708'
    pth = '/Users/tgray/Documents/SSX/Python/magcalib'
    # probe 1 r,t,z shots and same for probe 2,3
    p1 = ['2','1']
    p2 = ['3','4']
    p3 = []

    probe1Loc = {}
    probe2Loc = {}
    offset1 = (.434/.0254 - 3.72 - 10.25) * .0254 #434mm is length from bracket to first probe, - 3.72" for the distance between bracket and inner conflat edge while calibrating, - 10.25, which is the distance that the teflon standoff is from the center of the helmholtz coil  -  in mm
    offset2 = (.434/.0254 - 3.72 - 10.25) * .0254 #same as above, but for some reason I set it at 3.32"
    probe1Loc['r'] = (array([8, 23, 53, 83, 121, 160])- 8 )/ 1000. - offset1
    probe1Loc['t'] = (array([8, 23, 53, 83, 121, 160])- 8 )/ 1000. - offset1
    probe2Loc['r'] = (array([8, 23, 53, 83, 121, 160])- 8 )/ 1000. - offset2
    probe2Loc['z'] = (array([8, 23, 53, 83, 121, 160])- 8 )/ 1000. - offset2

    p1runs, p2runs, p3runs = [], [], []
    for p in p1:
        p1runs.append(date+'r'+p)
    for p in p2:
        p2runs.append(date+'r'+p)
# 	for p in p3:
# 		p3runs.append(date+'r'+p)
    ci1 = calibrationInfo(p1runs, numProbes = 6, dir = ['r', 't'])
    ci1['numChan'] = 12
    ci2 = calibrationInfo(p2runs, numProbes = 6, dir = ['r', 'z'])
    ci2['numChan'] = 12
# 	ci2 = calibrationInfo(p2runs)
# 	ci3 = calibrationInfo(p3runs)
    print p1
    probe1Calib, probe1Calib2 = sloppyCalibration(ci1, probe1Loc)
    probe2Calib, probe2Calib2 = sloppyCalibration(ci2, probe2Loc)
    # merge by channel dicts
    probe1Calib.update(probe2Calib)
    probeCalib = probe1Calib
    print "# probe 1 and 2 by channel"
    print probeCalib.__repr__()

    print "# probe 1 t, then z"
    print probe1Calib2['r'].__repr__()
    print probe1Calib2['t'].__repr__()
    print "# probe 2 t, then z"

    print probe2Calib2['r'].__repr__()
    print probe2Calib2['z'].__repr__()

def sloppyCalibration(calibInfo, probeLoc):
    """Simple calibration assuming the probes are perfectly aligned.

    Just corrects for magnitude and polarity."""
    fig = figure(1)
    fig.clear()
# 	b = {}
    calibData = {}
    calibData2 = {}
    for i, run in enumerate(calibInfo['runs']):

        # first we need to get the right field values from calculations
        dir = calibInfo['dir'][i]

        print dir, run
        b = sloppyCalibBFields(dir, probeLoc[dir])
        numProbes = len(b)	

        shotdata = readScope(run, calibInfo['scope'])

        current = getattr(shotdata, shotdata.names[calibInfo['currentChan'] - 1]) * (1000./0.00541/1.025)	# get the calibrated current for the shot

        currentMax = polyPeak(shotdata.time, current, timebase = 100)	# get the peak of the current

        magdata = sdr.magnetics_data(run)
        voltages = zeros(numProbes)
        for i in xrange(numProbes):
# 			print "probe: %s" % magdata.axesGroups[dir][i]
            tmpmag = polyPeak(magdata.time, getattr(magdata, magdata.axesGroups[dir][i]))
            voltages[i] = tmpmag / (2. * currentMax) #normalize to the current and double B from the helmholtz because there are 2 turns per coil
            calibData[magdata.axesGroups[dir][i]] =  voltages[i]/b[i]
            plot(getattr(magdata, magdata.axesGroups[dir][i]) / calibData[magdata.axesGroups[dir][i]])
        calibData2[dir] = voltages/b
    return calibData, calibData2

def sloppyCalibBFields(dir, probeLocs):
    """Calculates on axis B field of Helmholtz coil the locations given in the directions given."""
    numProbes = len(probeLocs)
    bb = zeros((numProbes,3))

# 		if numProbes == 8:
# 			j = -3.5+i*1.0
# 		elif numProbes == 4:
# 			j = 0.3125+i*1.0
    for i in xrange(numProbes):
        j = probeLocs[i]
        if dir == 'r':
            r = array([0, 0, j]) * 0.0254
            b = helmholtz(r)
            # r shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[2], 1*b[0], 1*b[1] ]	
        elif dir == 't':
            r = array([j, 0, 0]) * 0.0254
            b = helmholtz(r)
            # t shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[0], -1*b[2], 1*b[1] ]	
        elif dir == 'z':
            r = array([j, 0, 0]) * 0.0254	#same as t, but different transform below
            b = helmholtz(r)
            # z shots : x,y,z - > r,t,z
            bb[i] = [ 1*b[0], -1*b[1], -1*b[2] ]	
    if dir == 'r':
        bb = bb[:,0]
    elif dir == 't':
        bb = bb[:,1]
    elif dir == 'z':
        bb = bb[:,2]
    return bb


def calib080508():
    """Calibration from 8/5/08.

    Two full 8 channel quartz probes.  
    Probe 1 used cable Y.
    Probe 2 used cable U.

    Saves calibration data using numpy.  Use loadtxt('file') to retrieve."""

    date = '080508'
    pth = '/Users/tgray/Documents/SSX/Python/magcalib'
    # probe 1 r,t,z shots and same for probe 2
    p1 = ['4','5','6']
    p2 = ['3','2','1']

    p1runs, p2runs = [], []
    for p in p1:
        p1runs.append(date+'r'+p)
    for p in p2:
        p2runs.append(date+'r'+p)
    ci1 = calibrationInfo(p1runs, scope = '2X')
    ci2 = calibrationInfo(p2runs, scope = '2X')

    probe1Calib = calibration(ci1)
    savetxt(os.path.join(pth, 'calib%sp1.txt' % date), probe1Calib)
    probe2Calib = calibration(ci2)
    savetxt(os.path.join(pth, 'calib%sp2.txt' % date), probe2Calib)

def calib081308():
    """Calibration from 8/13/08.

    4 4 channel probes in gun region."""

    date = "081308"
    pth = '/Users/tgray/Documents/SSX/Python/magcalib'
    # probe 1 r,t,z shots and same for probe 2-4
    p1 = ['3', '1', '2']
    p2 = ['4', '5', '6']
    p3 = ['9', '8', '7']
    p4 = ['10', '11', '12']

    p1runs, p2runs, p3runs, p4runs = [], [], [], []
    for p in p1:
        p1runs.append(date+'r'+p)
    for p in p2:
        p2runs.append(date+'r'+p)
    for p in p3:
        p3runs.append(date+'r'+p)
    for p in p4:
        p4runs.append(date+'r'+p)
    ci1 = calibrationInfo(p1runs, scope = '2X', numProbes = 4)
    ci2 = calibrationInfo(p2runs, scope = '2X', numProbes = 4)
    ci3 = calibrationInfo(p3runs, scope = '2X', numProbes = 4)
    ci4 = calibrationInfo(p4runs, scope = '2X', numProbes = 4)

    # swap lines 60 and 61 for this calibration
    probe1Calib = calibration(ci1)
    savetxt(os.path.join(pth, 'calib%sp1.txt' % date), probe1Calib)
    probe2Calib = calibration(ci2)
    savetxt(os.path.join(pth, 'calib%sp2.txt' % date), probe2Calib)
    probe3Calib = calibration(ci3)
    savetxt(os.path.join(pth, 'calib%sp3.txt' % date), probe3Calib)
    probe4Calib = calibration(ci4)
    savetxt(os.path.join(pth, 'calib%sp4.txt' % date), probe4Calib)

def calib111109():
    """Calibration from 11/11/09.

    4 8 channel probes with new Jeff Santner integrators."""

    date = "111109"
    pth = '/Users/tgray/Documents/SSX/Python/magcalib'
    # probe 1 r,t,z shots and same for probe 2-4
    p1 = ['10', '11', '12']
    p2 = ['3', '2', '1']
    p3 = ['17', '16', '14']
    p4 = ['18', '19', '21']

    p1runs, p2runs, p3runs, p4runs = [], [], [], []
    for p in p1:
        p1runs.append(date+'r'+p)
    for p in p2:
        p2runs.append(date+'r'+p)
    for p in p3:
        p3runs.append(date+'r'+p)
    for p in p4:
        p4runs.append(date+'r'+p)
    ci1 = calibrationInfo(p1runs, scope = '3', probenum=1)
    ci2 = calibrationInfo(p2runs, scope = '3', probenum=2)
    ci3 = calibrationInfo(p3runs, scope = '3', probenum=3)
    ci4 = calibrationInfo(p4runs, scope = '3', probenum=4)

    # swap lines 60 and 61 for this calibration
    probe1Calib = calibration(ci1)
    savetxt(os.path.join(pth, 'calib%sp1.txt' % date), probe1Calib)
    probe2Calib = calibration(ci2)
    savetxt(os.path.join(pth, 'calib%sp2.txt' % date), probe2Calib)
    probe3Calib = calibration(ci3)
    savetxt(os.path.join(pth, 'calib%sp3.txt' % date), probe3Calib)
    probe4Calib = calibration(ci4)
    savetxt(os.path.join(pth, 'calib%sp4.txt' % date), probe4Calib)

def caliblongprobe():
    """Calibration for long probe - 1/27/10"""

    date = "012710"
    pth = '/Users/tgray/Documents/SSX/Python/magcalib'

    px = [2, 3, 4, 5, 6, 7, 8, 24, 10, 11, 12, 13, 14, 25, 16, 18, 20, 21, 22]
    py = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
        44, 45]
    px = [str(i) for i in px]
    py = [str(i) for i in py]

    pxruns, pyruns = [], []
    for p in px:
        pxruns.append(date + 'r' + p)
    for p in py:
        pyruns.append(date + 'r' + p)
    ci1 = {}
    ci1['scope'] = 3
    ci1['dir'] = 'x'
    ci1['runs'] = pxruns
    ci1['currentChan'] = 2
    ci2 = {}
    ci2['scope'] = 3
    ci2['dir'] = 'y'
    ci2['runs'] = pyruns
    ci2['currentChan'] = 2
    cf1 = quicklongcalib(ci1)
    cf2 = quicklongcalib(ci2)
    cf = hstack((cf1, cf2))
    savetxt(os.path.join(pth, 'calib%slong.txt' % date), cf)
    return cf

def quicklongcalib(calibInfo):
    dir = calibInfo['dir']
    scope = calibInfo['scope']
    chan = calibInfo['currentChan']
    if dir == 'x':
        bdir = 't'
    elif dir == 'y':
        bdir = 'r'
    B = sloppyCalibBFields(bdir, [0.0,])
    calibData = {}
    calibData2 = {}
    pnum = arange(19) + 1
    voltages = zeros(19)
    for i, run in enumerate(calibInfo['runs']):
        magchanname = 'm1' + dir + str(i + 1)
        print "Calibrating probe %s" % magchanname

        shotdata = sdr.scope_data(run, scope)
        current = getattr(shotdata, shotdata.names[chan-1]) * (1000./0.00541/1.025)	# get the calibrated current for the shot
        currentMax = polyPeak(shotdata.time, current, timebase = 100)	# get the peak of the current
        magdata = sdr.magnetics_data(run)
        magchan = getattr(magdata, magchanname)
        tmpmag = polyPeak(magdata.time, magchan)
        voltages[i] = tmpmag / (2. * currentMax) #normalize to the current and double B from the helmholtz because there are 2 turns per coil
    calibfactors = voltages / B
    return calibfactors

def calibcolorado():
    """Calibration for the colorado probe - July 2010"""

    date = "083110"
    pth = '/Users/tgray/Documents/SSX/Python/magcalib'

    # probe r,t,z shots
    # for colorado probe, it is labeled 1, 2, 3
    # 1 -> z
    # 2 -> r
    # 3 -> t
    # p    2     3     1
    p1 = ['19', '14', '10']
    p1runs = []
    for p in p1:
        p1runs.append(date+'r'+p)

    ci1 = calibrationInfo(p1runs, scope = '3', probenum=1, dir = ['1', '2', '3'], numProbes = 16)

    # probe locs around 0.  80mm long (16 * 5mm) - 2.5mm*2 for the edge to
    # center distance.
    ploc = 2.5+arange(16)*5-40
    ci1["probeLocs"] = ploc

    probe1Calib = calibrationCol(ci1, ['z', 'r', 't'])
    savetxt(os.path.join(pth, 'calib%scolorado.txt' % date), probe1Calib)
    return ci1
