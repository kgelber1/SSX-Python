import ssx_data_read_basic as sdr
import numpy as np
import mytools

from pylab import figure, axes, show
from pylab import ylabel, xlabel, title

import os
from matplotlib import pyplot as plt
import ssxdefaults_basic as ssxdef
import ssx_py_utils_basic as ssxutil

import iter_smooth as ism


def dens_calib(calib_shot, cutoff = 1, scope = '1'):

    """
    day - str, mmddyy, with date of the shot taken
    shot - int, shot number
    cutoff - float, cutoff frequency in MHz for lowpass filter used to improve data
    scope - str, scope


    returns:
        env - envelope, range of data for signal 1 and 2
        offset - minima of both signals
        phasediff - phase error between cos and sin, in radians

    """

    #load the data
    data = sdr.scope_data(calib_shot, scope)
    names = ['signal1', 'signal2']
    if scope == '2': channels = (1,2)
    if scope == '1': channels = (3,4)
    #channels = (3, 4)
    mult = (1, 1)
    data.setScopeProperties(names, channels, mult)

    cos, sin, t = data.signal1[data.time>50],data.signal2[data.time>50],data.time[data.time>50]
    # removes first bit of the data due to occasional noise at t = 0

    cos = mytools.lowpass_gfilter(t,cos,cutoff,return_mean=True)
    sin = mytools.lowpass_gfilter(t,sin,cutoff,return_mean=True)
    #apply lowpass filter to make the fits better
    print("stage 1 completed...")
    env = [max(cos)-min(cos), max(sin)-min(sin)]
    offset = [min(cos),min(sin)]
    #find ranges for the data

    cos = (cos-offset[0])/env[0]
    sin = (sin-offset[1])/env[1]
    #calibrate the data

    #find the phase error
    x_ym_list = [] #array of  xvals around the ymax pt
    weights_list = []
    for i in range(len(cos)):
        if sin[i] > 0.98:
            x_ym_list.append(cos[i])
            weights_list.append(np.abs((0.98-sin[i])/0.02))

    #find weighted average max pt.
    #x_ym = np.mean(x_ym_list)
    x_ym = np.dot(x_ym_list,weights_list)/sum(weights_list)
    #use cosine when sin = 1 to find phase error
    phasediff = np.arccos((x_ym-0.5)*2)-np.pi/2
    print("stage 2 completed...")
    return env, offset, phasediff

def interferometer(run, env, offset, phasediff, linfit = True, cutoff = 1, scope = '1', diam = 0.155575,  showPlot = False,
    writeFiles = False):
    """Interferometer analysis routine.

    Reads data from the scope data files and runs the analysis.  Inputting a
    proper diameter is necessary.

    Diameter is path length of interferometer in cm.  In oblate flux
    conservers it should be 50 cm.  For neck part of oblate conservers, it
    should be 15.5575 cm.

    run - str, mmddyy+'r'+shot number
    env - envelope output from densitycalibration, should be shot from same day
    offset - offset output from densitycalibration, should be shot from same day
    phasediff - phasediff output from densitycalibration , should be shot from same day
    linfit - whether or not to subtract a linear fit from the data. should set to false for
            calibration shots
    cutoff - float, lowpass filter cutoff frequency in MHz
    calib_shot - controls whether or not to subtract a linear fit to the density data
                 before t = 3 ms to remove phase drift. I.e. usually false for looking
                 at calibration data

    Written by Emma Suen-Lewis 2017, adapted from Tim Gray
    """
    print( 'Shot '+ run+': analyzing density')
    data = sdr.scope_data(run, scope)
    names = ['signal1', 'signal2']
    channels = (3, 4)
    mult = (1, 1)

    data.setScopeProperties(names, channels, mult)
    data.runname = run #save run name for folder naming in other functions
    data.cutoff = cutoff# save lowpass cutoff frequency

    t = data.time

    cos = mytools.lowpass_gfilter(t,data.signal1,cutoff,return_mean=True)
    sin = mytools.lowpass_gfilter(t,data.signal2,cutoff,return_mean=True)
    #load cos/sin data and apply lowpass filter

    cos = (cos-offset[0])/env[0]
    sin = (sin-offset[1])/env[1]
    #calibrate data with envelope and offsets

    #change ranges from [0,1] to [-1,1]
    for i in range(len(cos)):
        cos[i] = 2* cos[i]-1
        sin[i] = 2* sin[i]-1

    data.sig1cal = cos.copy()
    data.sig2cal = sin.copy()

    for i in range(len(cos)):
        cos[i] = (cos[i]/np.cos(phasediff)+(sin[i])*np.tan(phasediff))
    #remove effects of phase error

    #will output cosine and sine data (note that this doesn't include corrections
    #for low-frequency phase drift)
    data.cos = cos.copy()
    data.sin = sin.copy()

    phi =  np.arctan(sin/cos)
    #calculate phase in radians

    #get rid of jumps due to the range of arctan
    for i in range(len(phi)):
        if i < len(phi) - 1:
            if phi[i+1] - phi[i] > 3:
                phi[i+1:] = phi[i+1:] - np.pi
            if phi[i+1] - phi[i] < -3 :
                phi[i+1:] = phi[i+1:] + np.pi

    dphi = phi - np.mean(phi[0:10])
    #subtract out initial phase

    data.p = []
    #linear fit to data before firing
    if linfit == True:
        time_mask = t < 35
        p = np.polyfit(t[time_mask],ism.iter_smooth(dphi,loops=10,window_len=15)[time_mask],1)
        dphi = dphi - p[1] - p[0]*t
        data.p = p
    data.dphi = dphi

    # .155575 m is the ID of the small part of the flux conserver
    # mks
    #density = (dphi * 4 * pi * (3e8)**2 * 8.854e-                                                                                                                                        12 * 9.109e-31) /
    #((1.602e-19)**2 * 632.8e-9 * .155575)

    #cgs
    #density = (dphi * (3e10)**2 * 9.109e-28) / ((4.8032e-10)**2 * 632.8e-7 *
    #    diam)
    density=5.62e14*dphi/diam

    data.pathlength = diam
    data.density = density

    if showPlot or writeFiles:
        fig = figure(13, **ssxdef.f)
        fig.clear()
        a = axes()
        a.plot(data.time, data.density, 'k-')
        xlabel('Time (us)')
        ylabel('Density (#/cm$^3$)')
        title(data.shotname)
        if writeFiles:
            # make output directory
            fName = ssxutil.ssxPath('interferometer.png', 'output',
                data.runYear + '/' + data.runDate + '/' + run)
            dir, trash = os.path.split(fName)
            os.spawnlp(os.P_WAIT, 'mkdir', 'mkdir', '-p', dir)
            fig.savefig(fName)
        else:
            show()

    return data

def box():
    # makes a box and circle for the phase plots
    plt.plot([-1,1],[-1,-1],'r',[-1,1],[1,1],'r',[-1,-1],[-1,1],'r',[1,1],[-1,1],'r')
    t = np.arange(0,6.3,0.01)
    plt.plot(np.cos(t),np.sin(t),'r',linestyle = 'dashed')

def titlestr(run,cutoff,title):
    plt.suptitle(str(run)+' '+title+" (Lowpass Cutoff = "+str(cutoff*1000)+"kHz)", weight='bold')

def interferometerfigs(data, savefig = True, show = False, foldertitle = 'figures_interferometer',verbose = False):
    '''data - output of interferometer_new function
    verbose = True   will tell you what figure is being drawn

    Written by Emma Suen-Lewis 2017'''

    run = data.runname
    if verbose == True:
        print('Shot '+run+': starting figures')

    day = run[0:6]
    folder = foldertitle+'_'+day
    if savefig == True:
        if os.access(folder,os.F_OK) == False:
            os.mkdir(folder)

    cos = data.cos
    sin = data.sin
    dens = data.density
    t = data.time
    cutoff = data.cutoff
    p = data.p #coefficients for linear fit
    diam = data.pathlength

    #----------------------cos sin plot--------
    if verbose == True:
        print( '   starting plot: cosine/sine')

    plt.figure(figsize = (10,5))
    plt.axhline(1,color="grey")
    plt.axhline(0, color= "grey")
    plt.axhline(-1, color= "grey")
    plt.plot(t,cos,t,sin)
    titlestr(run,cutoff,'Cosine, Sine Data')
    plt.xlabel("Time $(\mu s)$",fontsize=12, weight='bold')
    plt.ylabel("$\\cos\\phi$ and $\\sin\\phi$")

    if savefig ==True:
        plt.savefig(folder+'\\'+run+'_cosine_sine.png')

    if show == True:
        plt.show()

    plt.close()

    #----------------phase space plot
    if verbose == True:
        print( '   starting plot: phase space')

    lim = 1.1

    plt.figure(figsize = (13,6))
    titlestr(run,cutoff,'Phase Space Data')
    plt.subplot(1,2,1)
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    plt.axhline(0,color="grey")
    plt.axvline(0,color="grey")
    box()
    plt.scatter(data.sig1cal,data.sig2cal,s=1,alpha=0.4,c = t, cmap = plt.cm.gist_rainbow)
    plt.axis("equal")
    plt.title("Before Phase Correction")

    plt.subplot(1,2,2)
    plt.axhline(0,color="grey")
    plt.axvline(0,color="grey")
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    box()
    plt.scatter(cos,sin,s=1,alpha=0.4,c = t, cmap = plt.cm.gist_rainbow)
    plt.axis("equal")
    plt.title("After Phase Correction")

    if savefig == True:
        plt.savefig(folder+'\\'+run+"_phase.png")
    if show == True:
        plt.show()

    plt.close()

    # ---------------density plot
    if verbose == True:
        print('   starting plot: density')

    plt.figure(figsize = (15,8))
    titlestr(run,cutoff,'Density Data')
    if p != []:
        plt.subplot(3,1,1)
        plt.axhline(0,color="grey")
        plt.plot(t,dens+(p[1]+p[0]*t)*5.62e14/diam)
        plt.plot(t,(p[1]+p[0]*t)*5.62e14/diam)
        plt.axvline(35,color="lightgrey",linestyle="dashed")
        plt.subplot(3,1,2)
    plt.axhline(0,color="grey")
    plt.plot(t,dens)
    plt.ylabel('Density (#/cm$^3$)',fontsize=12, weight='bold')

    if p!= []:
        #do analysis things
        plt.subplot(3,1,3)
        sm_dens = ism.iter_smooth(dens, loops = 50, window_len = 15)

        plt.axhline(0,color='grey')

        #find peak density
        peak = max(sm_dens)
        peakmask = np.abs(t - t[sm_dens == peak]) < 10
        peakavg = np.mean(sm_dens[peakmask])
        plt.axhline(peakavg,color="lightgrey",linestyle="dashed")

        for i in range(len(t)):
            if sm_dens[i] > peakavg/3:
                if t[i]>35:
                    start = t[i]
                    break
        plt.axvline(start,color="lightgrey",linestyle="dashed")

        indices = range(len(t))
        indices.reverse()
        for i in indices:
            if sm_dens[i] > peakavg/5*3:
                end = t[i]
                break
        plt.axvline(end,color="lightgrey",linestyle="dashed")

        plt.plot(t,sm_dens)
        plt.plot(t[peakmask],sm_dens[peakmask],color='red')

    if savefig == True:
        plt.savefig(folder+'\\'+run+"_density.png")
    if show == True:
        plt.show()

    plt.close()
