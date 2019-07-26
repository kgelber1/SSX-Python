from __future__ import division, print_function, absolute_import

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz
import itertools as it
from scipy.interpolate import interp1d


import interferometer_new_1 as ds
import density_Calshots as dcs
import ssx_basic as ssxd
import iter_smooth as ism
import basic_hiresmag as hdr
import basic_ids_1 as idsd
import mytools as my
import mjtools as mj
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from scipy import interpolate
from timeit import default_timer

"""
For reference:
matplotlib.figure.Figure.align_ylabels
"""


def get_nT(shot, t, ax1, ax2):
    """little plotter code to get the n and T

      -- KG 07/12/19
      """

    # minorLocator = AutoMinorLocator(10)       # leads to a single minor tick

    # Looks like the scope that is used for inferometer?
    scope_used='1'


    setting1 = '_wind-tunnel'
    setting2 = '_beta_Alfvenspeed'#'_WLH_GasDelay_550mus'
    setting3 = '_eos_windtunnel'
    title1 = r': WLH, 1 mW, 600 $\mu s$, Merging Configuration'

    env, offset, phasediff=ds.dens_calib(dcs.calshot(shot[:6]), scope= scope_used)
    a = env[0]/2
    b = env[1]/2

    def f(time, A, B): # this is your 'straight line' y=f(x)
        return A*time+B


    # plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)

    dens = ssxd.interferometer(shot, [a, b], scope = scope_used, showPlot=False)
    density= dens.density
    sm_density=ism.iter_smooth(density,loops=30, window_len=29)
    n = sm_density/(1e15)
    #popt, pcov = curve_fit(f, dens.time[0:2000], n[0:2000])
    #n = n + f(dens.time, *popt*1.3)
    timeN = dens.time
    ax1.plot(timeN, n, color='k',lw= 2)
    # plt.title(shot, fontsize=20, weight='bold')
    ax1.set_ylabel(r'n $(10^{15}\ cm^{-3})$',fontsize=20, weight='bold')
    # ax1.get_yaxis().set_label_coords(-0.11,0.6) # for aligning the y-labels in one line
    # plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
    # ax1.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 20)
    # ax1.tick_params(axis='x', which='minor', direction='in', length=5, width =1)

    ax1.set_xlim(t[0],t[-1])

    #########################################
    # ax2=plt.subplot(3,1,2)
    # ax2.text(0.07,0.92,'(b)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes)
    d=idsd.ids(shot)
    d.processIDS(times=[-2,125])
    timeT=d.time
    # This is where the errors happen?
    indices = np.where(d.kTFit.mask == False)[0] #Get indices of unmasked values
    Temp = d.kTFit.compressed() #Get unmasked values
    timeT = timeT[indices] #Adjust length of time array
    Terr = d.kTErr[indices]
    ax2.errorbar(timeT, Temp, Terr, fmt='None', ecolor='k',elinewidth=2,markeredgewidth=2,capsize=4)
    # ax2.plot(timeT, Temp, 'kx', color='k',ms = 8, mew=2)
    ax2.plot(timeT, Temp, color='k', linewidth=1)
    ax2.set_ylabel(r'T$_i\ (eV)$',fontsize=20, weight='bold')
    #ax2.set_xticklabels([])
    ax2.get_yaxis().set_label_coords(-0.11,0.6)
    plt.setp(ax2.spines.values(), linewidth=2)
    ax2.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 20)
    ax2.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
    #ax2.tick_params(axis='y', direction='in', length=7, width =2)
    plt.xlim(t[0],t[-1])
    # print(t[0], t[-1])



    #find the closest time t (regardless of repeats!)
    #and chuck the value at that time in
    # n_shifted = [n[my.tindex_center(timeN, time)] for time in t]
    # Temp_shifted = [interp_Temp[my.tindex_center(timeT, time)] for time in t ]
    interp_den = interp1d(timeN, n, kind='linear')
    interp_Temp = interp1d(timeT, Temp, kind='linear')
    plt.ylim(0,1.5*np.amax(interp_Temp(t)))

    return ax1, ax2, interp_den(t), interp_Temp(t)



def run_nBT(shots, day, t_min = 15, t_max = 100,  show = True, save = False, ylim = 35):
    """All of this based on Manjit's code, with a minimal amount of modifications.

    This is the main code to prodcue graphs of the temperature, density and
    magnetic feidl strength

      -- KG 06/28/19
      """

    minorLocator = AutoMinorLocator(10)       # leads to a single minor tick
    gs = gridspec.GridSpec(4,1)
    plt.rcParams['text.latex.preamble']=[r'\boldmath']

    # Looks like the scope that is used for inferometer?
    scope_used='1'

    path = 'data\\2019\\'+day+'\\Analyzed\\'

    setting1 = '_wind-tunnel'
    setting2 = '_beta_Alfvenspeed'#'_WLH_GasDelay_550mus'
    setting3 = '_eos_windtunnel'
    title1 = r': WLH, 1 mW, 600 $\mu s$, Merging Configuration'
    #title1 = ': WLH, 1 mW, 600 $\mu s$, coil scan at 25 kV'
    title2 = ': WLH, 1 mW, 600 $\mu s$, Merging Configuration'
    title3 = ': WLH, 1 mW, 600 $\mu s$, Merging Configuration'

    env, offset, phasediff=ds.dens_calib(dcs.calshot(day), scope= scope_used)
    a = env[0]/2
    b = env[1]/2
    # a = 1.312/2  for day  = '013017'
    # b = 1.234/2
    # a = 0.928/2
    # b = 0.978/2
    def f(time, A, B): # this is your 'straight line' y=f(x)
        return A*time+B

    for shot in shots:
        print( 'On Shot',shot)

        plt.close('all')
        # Adjust the spacing:
        fig=plt.figure(num=1,figsize=(8.5,10),facecolor='w',edgecolor='k')#, dpi=600)
        fig.subplots_adjust(top=0.95, bottom=0.11, left = 0.14, right=0.96, hspace=0.2)
        ax1=plt.subplot(3,1,1)

        plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)

        dens = ssxd.interferometer(day+'r'+str(shot), [a, b], scope = scope_used, showPlot=False)
        density= dens.density
        sm_density=ism.iter_smooth(density,loops=30, window_len=29)
        n = sm_density/(1e15)
        #popt, pcov = curve_fit(f, dens.time[0:2000], n[0:2000])
        #n = n + f(dens.time, *popt*1.3)
        timeN = dens.time
        plt.plot(timeN, n, color='k',lw= 2)
        plt.ylabel(r'n $(10^{15}\ cm^{-3})$',fontsize=20, weight='bold')
        plt.title(day+'r'+str(shot)+title1, fontsize=20, weight='bold')
        ax1.get_yaxis().set_label_coords(-0.11,0.6) # for aligning the y-labels in one line
        plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
        ax1.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 20)
        ax1.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
        ax1.xaxis.set_minor_locator(minorLocator)
        plt.xlim(t_min,t_max)

        #########################################
        ax2=plt.subplot(3,1,2)
        plt.text(0.07,0.92,'(b)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes)
        d=idsd.ids(day+'r'+str(shot))
        d.processIDS(times=[-2,125])
        timeT=d.time
        # This is where the errors happen?
        indices = np.where(d.kTFit.mask == False)[0] #Get indices of unmasked values
        Temp = d.kTFit.compressed() #Get unmasked values
        timeT = timeT[indices] #Adjust length of time array
        Terr = d.kTErr[indices]
        plt.errorbar(timeT, Temp, Terr, fmt='None', ecolor='k',elinewidth=2,markeredgewidth=2,capsize=4)
        plt.plot(timeT, Temp, 'kx', color='k',ms = 8, mew=2)
        plt.plot(timeT, Temp, color='k', linewidth=1)
        plt.ylabel(r'T$_i\ (eV)$',fontsize=20, weight='bold')
        #ax2.set_xticklabels([])
        ax2.get_yaxis().set_label_coords(-0.11,0.6)
        plt.setp(ax2.spines.values(), linewidth=2)
        ax2.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 20)
        ax2.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
        ax2.xaxis.set_minor_locator(minorLocator)
        #ax2.tick_params(axis='y', direction='in', length=7, width =2)
        plt.xlim(t_min,t_max)
        plt.ylim(0,ylim)

        #########################################
        ax3=plt.subplot(3,1,3)
        plt.text(0.07,0.92,'(c)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax3.transAxes)
        data=hdr.getquikData(day+'r'+str(shot))#x, y and z components of the 5th probe
        #calibration factors from lookup_4
        calib = [9.834933502238857272e+02, -1.263620982013806497e+03, -1.679900552773548725e+03]
        Bx=cumtrapz(data.unCalibData[0,4,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[0,4,:], 15), data.time)* calib[0]
        By=cumtrapz(data.unCalibData[1,4,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[1,4,:], 15), data.time)* calib[1]
        Bz=cumtrapz(data.unCalibData[2,4,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[2,4,:], 15), data.time)* calib[2]
        timeB=data.time-data.time[0]-2
        # timeB= data.time
        # make sure that the dimesntions match
        timeB = timeB[:len(Bx)]

        modB=np.sqrt((Bx**2)+(By**2)+(Bz**2))

        plt.plot(timeB, Bx,color='r', lw =1, label = 'B$_{x,6}$')
        plt.plot(timeB, By,color='g', lw =1, label = 'B$_{y,6}$')
        plt.plot(timeB, Bz,color='b', lw =1, label = 'B$_{z,6}$')
        plt.plot(timeB, modB,color='k', lw =2, label = '$|B|$')
        plt.legend().draggable()
        plt.ylabel(r'$|B|\ (G)$',fontsize=20, weight='bold')

        plt.setp(ax3.spines.values(), linewidth=2)
        ax3.get_yaxis().set_label_coords(-0.11,0.6)
        ax3.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 20)
        ax3.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
        ax3.xaxis.set_minor_locator(minorLocator)
        plt.xlim(t_min,t_max)
        plt.xlabel(r'$Time\ (\mu s)$',fontsize=20, weight='bold')

        ########## Saving Figure 1 ##################
        if show:
            plt.show()
        fName = path+day+'r'+str(shot)+setting1+'_plot.png'
        if save:
            fig.savefig(fName,dpi=600,facecolor='w',edgecolor='k')


def plot_nT(title,t_dens, den, t_Ti, Ti, ylim, d_err = None, t_err= None):
    """ As the name may sugest, plots density and temperatreu subplots with error bars

    ylim is the lim of the temp graph

    -- KG 07/19/19"""

    ax1= plt.subplot(2, 1, 1)
    plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)
    if(d_err != None):
        plt.errorbar(t_dens[::10], den[::10], d_err[::10], fmt='None', ecolor='k',elinewidth=.5,markeredgewidth=.5,capsize=.5, alpha  =0.2)
    plt.plot(t_dens, den, color='k',lw= 2)
    plt.ylabel(r'n $(10^{15}\ cm^{-3})$',fontsize=20, weight='bold')
    plt.title(title, fontsize=20, weight='bold')
    plt.xlim(15,100)

    ax2=plt.subplot(2,1,2)
    plt.text(0.07,0.92,'(b)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes)
    if(t_err != None):
        plt.errorbar(t_Ti, Ti, t_err, fmt='None', ecolor='k',elinewidth=1 ,markeredgewidth=1,capsize=2, alpha  =0.5)
    plt.plot(t_Ti, Ti, color='k', linewidth=1)
    plt.ylabel(r'T$_i\ (eV)$',fontsize=20, weight='bold')
    plt.xlim(15,100)
    plt.ylim(0,ylim)
    plt.show()



def plot_n(title,t_dens, den):
    # ax1= plt.subplot(2, 1, 1)
    # plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center')
    plt.plot(t_dens, den, color='k',lw= 2)
    plt.ylabel(r'n $(10^{15}\ cm^{-3})$',fontsize=20, weight='bold')
    plt.title(title, fontsize=20, weight='bold')
    plt.xlim(15,100)
    plt.show()

def get_stats(shots,day, show = True, save = True, ylim = 35):
    """
    A function to get the statistics for a series of shots
    Pass in the shots and the day, and it will plot the
    averaged density and temperature trace, as well as printing the
    density and temperature stats.

    There are some numbers hard-coded in, guesses of where the code should
    look for density and temperature peaks. You can also run run_nBT then
    look at all the outputs and updated these params if needed (block comments
    below)

      -- KG 06/28/19
    """

    num = len(shots)
    #stats I want:
    ne_t = np.zeros(num)# Time of ne peak
    ne_peak = np.zeros(num)# value of ne peak
    ne_pre = np.zeros(num)# value of ne before peak
    ne_post = np.zeros(num)# value of ne after peak (20 us average)
    t_dens = []
    den = []

    Ti_t = np.zeros(num)# Time of Ti p
    Ti_peak = np.zeros(num)# value of Ti p
    Ti_pre = np.zeros(num)# value of Ti before p
    Ti_post = np.zeros(num)# value of Ti after peak (20 us average)
    t_Ti = np.arange(-2,125)
    Ti = np.zeros(len(t_Ti))
    ave_over = np.zeros(len(t_Ti))


    scope_used='1'
    env, offset, phasediff=ds.dens_calib(dcs.calshot(day), scope= scope_used)
    a = env[0]/2
    b = env[1]/2

    for i,shot in enumerate(shots):
        #get density data:
        dens = ssxd.interferometer(day+'r'+str(shot), [a, b], scope = scope_used, showPlot=False)
        density= dens.density
        sm_density=ism.iter_smooth(density,loops=30, window_len=29)
        n = sm_density/(1e15)
        timeN = dens.time

        if i == 0:
            t_dens = timeN
            den = n
        else:
            den = [d + n[i] for i,d in enumerate(den,0)]


        # now get the peak between 20 and 30 us
        """ Edit density peak here """
        peak = np.array((20, 30))
        t_index, peakt, peakne= my.fix_array(peak, timeN, n)
        max_index = np.argmax(peakne)
        ne_t[i] =peakt[max_index]
        ne_peak[i] =peakne[max_index]

        #min in the 5mu before peak
        t_index, minT, minNe= my.fix_array(np.array((peak[0]-5, peak[0])), timeN, n)
        ne_pre[i] = np.min(minNe)

        #and the average value of the 20 mu after
        t_index, peakt, peakne= my.fix_array(np.array((peak[1], peak[1]+20)), timeN, n)
        ne_post[i] = np.average(peakne)

        # print(ne_t, ne_peak ,ne_pre , ne_post)
        ##########################################################

        #get temperature data
        d=idsd.ids(day+'r'+str(shot))
        d.processIDS(times=[-2,125])
        timeT=d.time
        indices = np.where(d.kTFit.mask == False)[0] #Get indices of unmasked values
        Temp = d.kTFit.compressed() #Get unmasked values
        timeT = timeT[indices] #Adjust length of time array
        Terr = d.kTErr[indices]

        # if i == 0:
        #     t_Ti = timeT
        #     Ti = Temp
        # print(timeT, t_Ti)
        j = 0 # index for the Ti data of the shot
        for k,t in enumerate(t_Ti):
            # jumping timesteps with missing values
            if(j>= len(timeT)):
                break
            if( np.absolute(timeT[j] - t) < .01):
                Ti[k] += Temp[j]
                ave_over[k] +=1
                # print(t, timeT[j])
                j+=1
            # Ti = [ti + Temp[i] for i,ti in enumerate(Ti) if i < len(Temp)]

        # now get the peak:
        """ Edit temperature peak here """
        t_index, peakt, peakTi= my.fix_array(np.array((35, 50)), timeT, Temp)
        max_index = np.argmax(peakTi)
        Ti_t[i] =peakt[max_index]
        Ti_peak[i] =peakTi[max_index]

        #the min in the 5mu before the peak
        minTi = my.local_min_before(Ti_t[i]-5, timeT, Temp)
        Ti_pre[i] = np.min(minTi)

        #and the average value after the peak
        t_index, peakt, peakti= my.fix_array(np.array((Ti_t[i]+5, Ti_t[i]+25)), timeT, Temp)
        Ti_post[i] = np.average(peakti)
        print("Shot", shot)

    #average
    den  = [d/num for d in den]
    for i in range(len(Ti)):
        if ave_over[i] > 0:
            Ti[i] = Ti[i]/ave_over[i]
            print(ave_over[i])
        else:
            Ti[i] = 0
    t_dens = t_dens[:len(den)]
    t_Ti= t_Ti[:len(Ti)]

    lens = np.sqrt(num)
    def stats(arr):
        return (np.mean(arr),np.std(arr, dtype=np.float64)/lens)

    if show:
        title = day + ' - averaged'
        plot_nT(title, t_dens, den, t_Ti, Ti, ylim)
        print("Density Stats:")
        print("\tAverage time of peak:\n\t %.1f +/- %2.1f us" %(stats(ne_t)))
        print("\tAverage Value of peak:\n\t %.1f +/- %2.1f e15" %(stats(ne_peak)))
        print("\tAverage value before peak:\n\t %.2f +/- %2.2f e15" %(stats(ne_pre)))
        print("\tAverage value after peak:\n\t %.1f +/- %2.1f e15" %(stats(ne_post)))

        print("Temp Stats:")
        print("\tAverage time of peak:\n\t %.1f +/- %2.1f us" %(stats(Ti_t)))
        print("\tAverage value of peak:\n\t %.1f +/- %2.1f eV" %(stats(Ti_peak)))
        # print(Ti_pre)
        print("\tAverage value before peak:\n\t %.1f +/- %2.1f eV" %(stats(Ti_pre)))
        print("\tAverage value after peak:\n\t %.1f +/- %2.1f eV" %(stats(Ti_post)))

    if save:
        #haven't wrote yet but you could add a function to save the data here
        pass


def get_stats_err(shots,day, show = True, save = True, ylim = 35):
    """
    A function to get the statistics for a series of shots -
    INCLUDES ERROR BARS.
    Pass in the shots and the day, and it will plot the
    averaged density and temperature trace, as well as printing the
    density and temperature stats.

    Code gets messy becuase Ti data is not masked, bad elements are removed
    so I had to re-fill in None values


    There are some numbers hard-coded in, guesses of where the code should
    look for density and temperature peaks. You can also run run_nBT then
    look at all the outputs and updated these params if needed (block comments
    below)

      -- KG 06/28/19
    """

    num = len(shots)
    #stats I want:
    all_ne = [i for i in range(num)]
    all_Ti = [i for i in range(num)]
    ne_t = np.zeros(num)# Time of ne peak
    ne_peak = np.zeros(num)# value of ne peak
    ne_pre = np.zeros(num)# value of ne before peak
    ne_post = np.zeros(num)# value of ne after peak (20 us average)
    t_dens = []
    den = []

    # Ti_stats = [np.zeros(num) for i in range(4)]
    Ti_t = np.zeros(num)# Time of Ti p
    Ti_peak = np.zeros(num)# value of Ti p
    Ti_pre = np.zeros(num)# value of Ti before p
    Ti_post = np.zeros(num)# value of Ti after peak (20 us average)
    t_Ti = np.arange(-2,125)
    Ti = np.zeros(len(t_Ti))
    ave_over = np.zeros(len(t_Ti))


    scope_used='1'
    env, offset, phasediff=ds.dens_calib(dcs.calshot(day), scope= scope_used)
    a = env[0]/2
    b = env[1]/2

    for i,shot in enumerate(shots):
        #get density data:
        dens = ssxd.interferometer(day+'r'+str(shot), [a, b], scope = scope_used, showPlot=False)
        density= dens.density
        sm_density=ism.iter_smooth(density,loops=30, window_len=29)
        n = sm_density/(1e15)
        timeN = dens.time

        if i == 0:
            t_dens = timeN
            den = n
        else:
            den = [d + n[i] for i,d in enumerate(den,0)]
        all_ne[i] = [x for x in n]

        # now get the peak between 20 and 30 us
        """ Edit density peak here """
        peak = np.array((20, 30))
        t_index, peakt, peakne= my.fix_array(peak, timeN, n)
        max_index = np.argmax(peakne)
        ne_t[i] =peakt[max_index]
        ne_peak[i] =peakne[max_index]

        #min in the 5mu before peak
        t_index, minT, minNe= my.fix_array(np.array((peak[0]-5, peak[0])), timeN, n)
        ne_pre[i] = np.min(minNe)

        #and the average value of the 20 mu after
        t_index, peakt, peakne= my.fix_array(np.array((peak[1], peak[1]+20)), timeN, n)
        ne_post[i] = np.average(peakne)

        # print(ne_t, ne_peak ,ne_pre , ne_post)
        ##########################################################

        #get temperature data
        d=idsd.ids(day+'r'+str(shot))
        d.processIDS(times=[-2,125])
        timeT=d.time
        indices = np.where(d.kTFit.mask == False)[0] #Get indices of unmasked values
        Temp = d.kTFit.compressed() #Get unmasked values
        timeT = timeT[indices] #Adjust length of time array
        Terr = d.kTErr[indices]

        j = 0 # index for the Ti data of the shot
        all_Ti[i] = [0 for x in t_Ti]
        for k,t in enumerate(t_Ti):
            # jumping timesteps with missing values
            if(j>= len(timeT)):
                break
            if( np.absolute(timeT[j] - t) < .01):
                Ti[k] += Temp[j]
                ave_over[k] +=1
                all_Ti[i][k] = Temp[j]
                j+=1

        # now get the peak:
        """ Edit temperature peak here """
        t_index, peakt, peakTi= my.fix_array(np.array((35, 50)), timeT, Temp)
        max_index = np.argmax(peakTi)
        Ti_t[i] =peakt[max_index]
        Ti_peak[i] =peakTi[max_index]

        #the min in the 5mu before the peak
        minTi = my.local_min_before(Ti_t[i]-5, timeT, Temp)
        Ti_pre[i] = np.min(minTi)

        #and the average value after the peak
        t_index, peakt, peakti= my.fix_array(np.array((Ti_t[i]+5, Ti_t[i]+25)), timeT, Temp)
        Ti_post[i] = np.average(peakti)
        print("Shot", shot)

    #average
    den  = [d/num for d in den]
    for i in range(len(Ti)):
        if ave_over[i] > 0:
            Ti[i] = Ti[i]/ave_over[i]
        else:
            Ti[i] = 0
            ave_over[i] = 1
    t_dens = t_dens[:len(den)]
    t_Ti= t_Ti[:len(Ti)]

    lens = np.sqrt(num)
    d_err = []
    t_err = []


    all_ne = [list(i) for i in zip(*all_ne)]
    all_Ti = list(it.zip_longest(*all_Ti, fillvalue = None))
    for i,t in enumerate(t_dens):
        x = np.std(all_ne[i][:])/lens
        d_err = np.append(d_err, x)

    for i in range(len(t_Ti)):
        x = np.nanstd(all_Ti[i][:])/lens
        t_err = np.append(t_err, x)

    def stats(arr):
        return (np.mean(arr),np.std(arr, dtype=np.float64)/lens)

    if show:
        title = day + ' - averaged'
        plot_nT(title, t_dens, den, t_Ti, Ti, ylim, d_err, t_err)
        print("Density Stats:")
        print("\tAverage time of peak:\n\t %.1f +/- %2.1f us" %(stats(ne_t)))
        print("\tAverage Value of peak:\n\t %.1f +/- %2.1f e15" %(stats(ne_peak)))
        print("\tAverage value before peak:\n\t %.2f +/- %2.2f e15" %(stats(ne_pre)))
        print("\tAverage value after peak:\n\t %.1f +/- %2.1f e15" %(stats(ne_post)))

        print("Temp Stats:")
        print("\tAverage time of peak:\n\t %.1f +/- %2.1f us" %(stats(Ti_t)))
        print("\tAverage value of peak:\n\t %.1f +/- %2.1f eV" %(stats(Ti_peak)))
        print("\tAverage value before peak:\n\t %.1f +/- %2.1f eV" %(stats(Ti_pre)))
        print("\tAverage value after peak:\n\t %.1f +/- %2.1f eV" %(stats(Ti_post)))

    if save:
        #haven't wrote yet but you could add a function to save the data here
        pass

def main():
    #change your params!
    day = '072419'
    first_shot = 11
    # last_shot = 44
    last_shot = 43
    bad_shots = [14]
    # bad_shots = [20,23,39]


    all_shots = np.arange(first_shot,last_shot+1)
    shots = [shot for shot in all_shots if shot not in bad_shots]
    # shots = [13]
    run_nBT(shots, day, t_min = 15, t_max = 100, show = False, save = True, ylim = 85)
    # get_stats_err(shots, day, ylim = 40)

if __name__ == '__main__':
    main()
