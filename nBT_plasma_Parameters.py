# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 11:08:14 2016
@author: Manjit Kaur
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz
# reload(ds)
# reload(dcs)
import interferometer_new_1 as ds
import density_Calshots as dcs
import ssx_basic as ssxd
import iter_smooth as ism
import basic_hiresmag as hdr
import basic_ids_1 as idsd
import mytools as my
import mjtools as mj
#from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy import interpolate
from timeit import default_timer
#import csv

def main():

    minorLocator = AutoMinorLocator(10)       # leads to a single minor tick
    #ml = MultipleLocator(5)
    #minorLocator   = MultipleLocator(10)

    gs = gridspec.GridSpec(4,1)
    # plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble']=[r'\boldmath']

    # where day
    day = '061219'
    shot_range = [6,15] ## 87 code broke
    # Looks like the scope that is used for inferometer?
    scope_used='1'

    t0 = 40
    tf = 70
    path = 'data\\2019\\'+day+'\\'+day+'-Analysed\\'
    BEGIN = 35
    END = 100
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

    for shot in range(shot_range[0],shot_range[1]+1):
        print( 'On Shot',shot)
    #    if shot==11 or 12: continue
        #########################################
        plt.close('all')
        fig=plt.figure(num=1,figsize=(8.5,10),facecolor='w',edgecolor='k')#, dpi=600)
        #plt.clf()
        fig.subplots_adjust(top=0.95, bottom=0.11, left = 0.14, right=0.96, hspace=0.2)

        ax1=plt.subplot(3,1,1)
        plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)

        density = np.zeros([20002])
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
        #ax1.set_xticklabels([])
        #ax1.set_yticks([0.5,1.0,1.5,2.0])
        plt.xlim(20,100)
        # plt.ylim(0, 2)

        #########################################
        ax2=plt.subplot(3,1,2)
        plt.text(0.07,0.92,'(b)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes)
        d=idsd.ids(day+'r'+str(shot))
        d.processIDS(times=[-2,125])
        timeT=d.time
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
        plt.xlim(20,100)
        plt.ylim(0,35)

        #########################################
        ax3=plt.subplot(3,1,3)
        plt.text(0.07,0.92,'(c)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax3.transAxes)
        data=hdr.getquikData(day+'r'+str(shot))#x, y and z components of the 5th probe
        ##### changing from 15x, 15y and 15z to 6x, 6y and 6z
        # Bx=cumtrapz(data.unCalibData[0,14,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[0,14,:], 15), data.time)/7.63e-4
        # By=cumtrapz(data.unCalibData[1,14,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[1,14,:], 15), data.time)/8.06e-4
        # Bz=cumtrapz(data.unCalibData[2,2,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[2,2,:], 15), data.time)/6.76e-4
        # I think the get moving average is accounting for some ofset...

        # Bx=cumtrapz(data.unCalibData[0,5,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[0,5,:], 15), data.time)/7.63e-4
        # By=cumtrapz(data.unCalibData[1,5,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[1,5,:], 15), data.time)/8.06e-4
        # Bz=cumtrapz(data.unCalibData[2,1,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[2,1,:], 15), data.time)/6.76e-4
        Bx=cumtrapz(data.unCalibData[0,5,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[0,5,:], 15), data.time)* -1.088718234303543795e+03
        By=cumtrapz(data.unCalibData[1,5,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[1,5,:], 15), data.time)* -1.119343007186590285e+03
        Bz=cumtrapz(data.unCalibData[2,1,:]-mj.get_gaussian_moving_avg(data.time, data.unCalibData[2,1,:], 15), data.time)* -1.422749719357707363e+03
        timeB=data.time-data.time[0]-2
        timeB= timeB[:-1]
    #    Bx=cumtrapz(data.unCalibData[0,14,1500:] - np.mean(data.unCalibData[0,14,4700:6500]), data.time[1500:])/7.63e-4
    #    By=cumtrapz(data.unCalibData[1,14,1500:] - np.mean(data.unCalibData[1,14,4700:6500]), data.time[1500:])/8.06e-4
    #    Bz=cumtrapz(data.unCalibData[2,2,1500:] - np.mean(data.unCalibData[2,2,4700:6500]), data.time[1500:])/6.76e-4
    #    timeB=data.time-data.time[0]-2
    #    timeB= timeB[1500:-1]
    ########## detrending the Bx, By, Bz signals ################
    #    poptx, pcovx = curve_fit(f, timeB, Bx)
    #    Bx = Bx - f(timeB, *poptx)
    #    popty, pcovy = curve_fit(f, timeB, By)
    #    By = By - f(timeB, *popty)
    #    poptz, pcovz = curve_fit(f, timeB, Bz)
    #    Bz = Bz - f(timeB, *poptz)
    #    plt.plot(timeB, Bx,color='r', lw =2)
    #    plt.plot(timeB, By,color='b', lw =2)
    #    plt.plot(timeB, Bz,color='g', lw =2)

        modB=np.sqrt((Bx**2)+(By**2)+(Bz**2))
        #poptz, pcovz = curve_fit(f, timeB, modB)
        #plt.plot(timeB, modB, linewidth=2)
    #    plt.plot(timeB, f(timeB, *poptz),color='g', lw =2)
    #    modB_diff = np.mean(data.unCalibData[0,14,4700:4700])
    #    modB = modB-modB_diff
        plt.plot(timeB, Bx,color='r', lw =1, label = 'B$_{x,6}$')
        plt.plot(timeB, By,color='g', lw =1, label = 'B$_{y,6}$')
        plt.plot(timeB, Bz,color='b', lw =1, label = 'B$_{z,6}$')
        plt.plot(timeB, modB,color='k', lw =2, label = '$|B|$')
        plt.legend().draggable()
        plt.ylabel(r'$|B|\ (T)$',fontsize=20, weight='bold')
        #ax.set_xticklabels([])
        #ax.xaxis.set_tick_params(labeltop='on',labelsize=20)
        plt.setp(ax3.spines.values(), linewidth=2)
        ax3.get_yaxis().set_label_coords(-0.11,0.6)
        ax3.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 20)
        ax3.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
        ax3.xaxis.set_minor_locator(minorLocator)
        plt.xlim(20,100)
        plt.xlabel(r'$Time\ (\mu s)$',fontsize=20, weight='bold')

        ########## Saving Figure 1 ##################
        plt.show()
        fName = path+day+'r'+str(shot)+setting1+'_plot.png'
        # fig.savefig(fName,dpi=600,facecolor='w',edgecolor='k')

    '''
    #########$$$$$$$$$$$############# Figure 2 stuff ##########$$$$$$$$$$$$$##############
        ################## Interpolating density data #######################
        start = default_timer()
        print("Interpolating density data...")
        dens_spline = interpolate.interp1d(timeN,n,kind = 'slinear')
    #    duration = default_timer() - start
    #    print "time: %s s" % str(duration)

        ################## Interpolating IDS data ###################
        print("Interpolating temperatures...")
        temp_spline = interpolate.interp1d(timeT,Temp,kind = 'cubic')
        #Interpolate temperature error
        print("Interpolating temperature errors...")
        temp_err = interpolate.interp1d(timeT,Terr)
        #Save the time for temperatures
        tindex, timeT = my.fix_time_bounds(np.array((BEGIN,END)),timeT)
        samplingFrequency = 65 #Use frequency at which B is sampled 65MHz
        numPoints = int((END - BEGIN)*samplingFrequency)
        TIME = np.linspace(timeT[0],timeT[-1],numPoints)

        ###################### Interpolating magnetics data ###################
        print("Interpolating B...")
        modB_spline = interpolate.interp1d(timeB, modB,kind = 'slinear')
        duration = (default_timer() - start)/60
        print "time: %s min" % str(duration)


       ######### fixing bounds for EOS ###########
        t_index, t = my.fix_bounds(np.array((t0, tf)), TIME)
        t0 = t_index[0]
        tf = t_index[-1]
        TIME = TIME[t0:tf]

        ############## beta and Alfven speed ##############
        density = dens_spline(TIME)
        Temperature = temp_spline(TIME)
        magB = modB_spline(TIME)#*10**-4

        #beta = 3*(10**4)*density*Temperature/(magB)**2
        beta = 1.26*(10**2)*(density*10**21)*(Temperature*1.6*10**-19)/(magB**2)
        v_alf = 0.069*magB/np.sqrt(density)

    ############ plotting beta parameter and Alfven speed ################
        fig2=plt.figure(num=2,figsize=(8.5,10),facecolor='w',edgecolor='k')#, dpi=600)
        fig2.subplots_adjust(top=0.92, bottom=0.14, left = 0.14, right=0.96, hspace=0)
        color = 'blue'
        ax1=plt.subplot(2,1,1)
        plt.plot(TIME, beta, 'kx', color='firebrick',ms = 3, mew=2)
        plt.ylabel(r'$beta$',fontsize=20, weight='bold')
        plt.title(day+'r'+str(shot)+title2, fontsize=20, weight='bold')
        ax1.set_xticklabels([])
        ax1.xaxis.set_minor_locator(minorLocator)# for creating minor ticks
        ax1.get_yaxis().set_label_coords(-0.1,0.6)# for arranging the labels in one line
        plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
        ax1.tick_params(axis='both', direction='in', length=10, width =3, labelsize = 16)
    #    plt.xlim(50.2,53.9)

        ax2=plt.subplot(2,1,2)
        plt.plot(TIME, v_alf, 'kx', color='firebrick',ms = 3, mew=2)
        plt.ylabel(r'$v_{Alfven}$',fontsize=20, weight='bold')
        ax2.xaxis.set_minor_locator(minorLocator)# for creating minor ticks
        ax2.get_yaxis().set_label_coords(-0.1,0.6)# for arranging the labels in one line
        plt.setp(ax2.spines.values(), linewidth=2)#changing the axis linewidth
        ax2.tick_params(axis='both', direction='in', length=10, width =3, labelsize = 16)
        ax2.tick_params(axis='y', direction='in', length=10, width =3)
    #    plt.xlim(50.2,53.9)
    #    plt.ylim(0, 0.5)
        plt.xlabel(r'$Time\ \ (\mu s)$',fontsize=20, weight='bold')
        #plt.show()
        #############Saving figure ################
        fName2 = path+day+'r'+str(shot)+setting2+'.png'
        fig2.savefig(fName2,dpi=600,facecolor='w',edgecolor='k')


    #########$$$$$$$$$$$############# Figure 2 stuff ##########$$$$$$$$$$$$$##############
    ##################### plotting the Equations of state ################
        quantity1 = Temperature*density**(-2/3)
        quantity3 = Temperature*(magB**2)/(density**2)
        quantity2=  Temperature/magB
        quantity1 = quantity1/max(quantity1)
        quantity2 = quantity2/max(quantity2)
        quantity3 = quantity3/max(quantity3)
        qErr2= 0.15*quantity1
        qErr3= 0.3*quantity2
        qErr1= 0.12*quantity3

        ############ plotting the MHD and double adiabats EOS ################
        fig3=plt.figure(num=3)
        fig3.subplots_adjust(top=0.92, bottom=0.14, left = 0.14, right=0.96, hspace=0)
        color = 'blue'
        ax1=plt.subplot(3,1,1)
        ax1.fill_between(TIME,quantity1*(1+0.12), quantity1*(1-0.12),alpha = 0.8, color = color)
        plt.plot(TIME, quantity1,color='k',zorder=10, lw=2)#, label=r'3D compression ($\gamma = 5/3$)')
        plt.ylabel(r'$P/n^{5/3}$',fontsize=20, weight='bold')
        plt.title(day+'r'+str(shot)+title3, fontsize=20, weight='bold')
        #plt.legend()
        ax1.set_xticklabels([])
        ax1.xaxis.set_minor_locator(minorLocator)# for creating minor ticks
        ax1.get_yaxis().set_label_coords(-0.1,0.6)# for arranging the labels in one line
        plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
        ax1.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 16)
        ax1.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
        #ax.tick_params(axis='y', direction='in', length=7, width =2)
    #    plt.xlim(t0, tf)
        #plt.ylim(y.min()-0.1,y.max()+0.1)

        ax2=plt.subplot(3,1,2)
        ax2.fill_between(TIME,quantity2*(1+0.14), quantity2*(1-0.14),alpha = 0.8, color = color)
        plt.plot(TIME, quantity2,color='k',zorder=10,lw=2)#, 'kx', color='firebrick',ms = 3, mew=2)
        plt.ylabel(r'${P}\slash{nB}$',fontsize=20, weight='bold')
        ax2.set_xticklabels([])
        ax2.xaxis.set_minor_locator(minorLocator)# for creating minor ticks
        ax2.get_yaxis().set_label_coords(-0.1,0.6)# for arranging the labels in one line
        plt.setp(ax2.spines.values(), linewidth=2)#changing the axis linewidth
        ax2.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 16)
        ax2.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
        #ax.tick_params(axis='y', direction='in', length=7, width =2)
    #    plt.xlim(t0, tf)
        #plt.ylim(0.25,0.8)
    ##################
        ax3=plt.subplot(3,1,3)
        ax3.fill_between(TIME,quantity3*(1+0.3), quantity3*(1-0.3),alpha = 0.8, color = color)
        plt.plot(TIME, quantity3,color='k',zorder=10,lw=2)#, 'kx', color='firebrick',ms = 3, mew=2)
        plt.ylabel(r'$PB^2 \slash{n^3}$',fontsize=20, weight='bold')
        ax3.get_yaxis().set_label_coords(-0.1,0.6)# for arranging the labels in one line
        plt.setp(ax3.spines.values(), linewidth=2)#changing the axis linewidth
        ax3.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 16)
        ax3.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
        #ax.xaxis.set_minor_locator(minorLocator)# for creating minor ticks
    #    plt.xlim(t0, tf)
        #plt.ylim(0, 0.5)
        plt.xlabel(r'$Time\ \ (\mu s)$',fontsize=20, weight='bold')
        plt.show()
        #############Saving figure ################
        fName3 = path+day+'r'+str(shot)+'_EOS'+setting3+'.png'
        fig3.savefig(fName3,dpi=600,facecolor='w',edgecolor='k')
'''


def main2():
    x = np.arange(0,10)
    y = np.arange(0,20)
    y = y[::2]
    x = x[:len(y)]
    # plt.plot(x,y)
    # plt.show()
    # minorLocator = AutoMinorLocator(10)       # leads to a single minor tick
    #ml = MultipleLocator(5)

    # gs = gridspec.GridSpec(4,1)
    # plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble']=[r'\boldmath']

    # where day
    day = '061219'
    shot_range = [6,15] ## 87 code broke
    # Looks like the scope that is used for inferometer?
    scope_used='1'

    # t0 = 40
    # tf = 70
    # path = 'data\\2019\\'+day+'\\'
    # BEGIN = 35
    # END = 100
    title1 = r': WLH, 1 mW, 600 $\mu s$, Merging Configuration'


    # env, offset, phasediff=ds.dens_calib(dcs.calshot(day), scope= scope_used)
    # a = env[0]/2
    # b = env[1]/2
    # a = 1.312/2  for day  = '013017'
    # b = 1.234/2
    # a = 0.928/2
    # b = 0.978/2
    # def f(time, A, B): # this is your 'straight line' y=f(x)
    #     return A*time+B

    #on;y plot one shot
    for shot in range(shot_range[0],shot_range[0]+1):
        print( 'On Shot',shot)
    #    if shot==11 or 12: continue
        #########################################
        plt.close('all')
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

        # fig=plt.figure(num=1,figsize=(8.5,10),facecolor='w',edgecolor='k')#, dpi=600)
        #plt.clf()
        # fig.subplots_adjust(top=0.95, bottom=0.11, left = 0.14, right=0.96, hspace=0.2)

        # ax1=plt.subplot(3,1,1)
        # plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)

        # density = np.zeros([20002])
        # dens = ssxd.interferometer(day+'r'+str(shot), [a, b], scope = scope_used, showPlot=False)
        # density= dens.density
        # sm_density=ism.iter_smooth(density,loops=30, window_len=29)
        # n = sm_density/(1e15)
        #popt, pcov = curve_fit(f, dens.time[0:2000], n[0:2000])
        #n = n + f(dens.time, *popt*1.3)
        # timeN = dens.time
        # ax1.plot(timeN, n, color='k',lw= 2)
        ax1.plot(x,y)
        # plt.ylabel(r'n $(10^{15}\ cm^{-3})$',fontsize=20, weight='bold')
        # plt.title(day+'r'+str(shot)+title1, fontsize=20, weight='bold')
        # ax1.get_yaxis().set_label_coords(-0.11,0.6) # for aligning the y-labels in one line
        # plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
        # ax1.tick_params(axis='both', direction='in', length=7, width =2, labelsize = 20)
        # ax1.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
        #ax1.set_xticklabels([])
        #ax1.set_yticks([0.5,1.0,1.5,2.0])
        # plt.show()
        # plt.ylim(0, 2)


if __name__ == '__main__':
    main()
    # main2()
