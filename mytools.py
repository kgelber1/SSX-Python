# import prime
import time
import numpy
import scipy
import pylab
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors
import colorsys
from scipy import optimize
from scipy.ndimage import gaussian_filter
import itertools
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pylab as plt
from scipy import signal
import math
from scipy.signal import blackman,bartlett,hanning,hamming


def normalized(array):# new array will vary in the range 0 to 1
    normalized_value= (array-min(array))/(max(array)-min(array))
    return normalized_value


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    y = signal.lfilter(b, a, data)
    return y

def windowed_fft(array,time,window='None'):#time should be in seconds
    #Size of array
    Nw = array.shape[0]
    #Calculate time step (assumed to be in seconds)
    dt = time[1]-time[0]
    #prefactor
    #print 'dt = ',dt
    prefactor = dt
    #Calculate array of frequencies, shift
    w = np.fft.fftfreq(Nw,dt)
    w0 = np.fft.fftshift(w)
    #make window
    #blackman window
    if window == 'blackman':
        bwin = blackman(Nw) #pretty good
    if window == 'hanning':
        bwin = hanning(Nw) #pretty good
    if window == 'hamming':
        bwin = hamming(Nw) #not as good
    if window == 'bartlett':
        bwin = bartlett(Nw) #pretty good
    if window == 'None':
        bwin = 1.0
    #Calculate FFT
    aw = prefactor*np.fft.fft(array*bwin)
    aw0 = np.fft.fftshift(aw)
    #Calcuate Phase
    phase = np.angle(aw)
    phase0 = np.fft.fftshift(phase)
    #Adjust arrays if not div by 2
    if not np.mod(Nw,2):
        w0 = np.append(w0,-w0[0])
        aw0 = np.append(aw0,-aw0[0])
        phase0 = np.append(phase0,-phase0[0])
    #Cut FFTs in half
    Nwi = Nw/2
    w2 = w0[Nwi:]
    aw2 = aw0[Nwi:]
    phase2 = phase0[Nwi:]
    comp = aw
    pwr = (np.abs(aw2))**2
    pwr2 = (np.abs(aw))**2
    mag = np.sqrt(pwr)
    cos_phase = np.cos(phase2)
    freq = w2
    freq2 = w
    return freq,freq2,comp,pwr,mag,phase2,cos_phase,dt # plt.loglog(freq, pwr) # freq in MHz

def semilogplot(fignum, x, y, linestyle, xlabel, ylabel, title, label, xlim, ylim, color,
                fontsize, linewidth, fName, marker, scatter = True, horzline = True, saveFig=False):
    '''
    fignum, x, y, xlabel, ylabel, title, xlim, ylim, color, fontsize, fName, horzline = True, saveFig=False
    '''
    fig=plt.figure(num = fignum,figsize=(8.5,6),facecolor='w',edgecolor='k')#, dpi=600)
    fig.subplots_adjust(top=0.94, bottom=0.11, left = 0.14, right=0.97, hspace=0)

    ax1=plt.subplot(1,1,1)
    # plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)
    if horzline ==True:
        ax1.axhline(linewidth=1, color='k', linestyle = 'dotted')
    if scatter ==True:
        plt.semilogx(x, y, marker = marker, color = color, linewidth = linewidth,label = label)
    else:
        plt.semilogx(x, y, linestyle = linestyle, color = color,linewidth = linewidth,label = label)

    # plt.semilogx(x, y, linestyle = linestyle, color = color,lw= linewidth)
    plt.ylabel(ylabel,fontsize=fontsize)#, weight='bold')
    plt.xlabel(xlabel,fontsize=fontsize)#, weight='bold')
    plt.title(title, fontsize=fontsize)#, weight='bold')
    # ax1.get_yaxis().set_label_coords(-0.11,0.5)
    plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
    ax1.tick_params(axis='x', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
    ax1.tick_params(axis='y', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='y', which='minor', direction='in', length=5, width =1)
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.show()
    if saveFig == True:
        # fName = day+'r'+str(shot)+'_parameters.png'
        fig.savefig(fName,dpi=600,facecolor='w',edgecolor='k')

def logplot(fignum, x, y, linestyle, xlabel, ylabel, title, label, xlim, ylim, color,
            fontsize, linewidth, fName, marker, scatter = True, horzline = True, saveFig=False):
    '''
    fignum, x, y, xlabel, ylabel, title, xlim, ylim, color, fontsize, fName, horzline = True, saveFig=False
    '''
    fig=plt.figure(num = fignum,figsize=(8.5,6),facecolor='w',edgecolor='k')#, dpi=600)
    fig.subplots_adjust(top=0.94, bottom=0.11, left = 0.14, right=0.97, hspace=0)

    ax1=plt.subplot(1,1,1)
    # plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)
    if horzline ==True:
        ax1.axhline(linewidth=1, color='k', linestyle = 'dotted')
    if scatter ==True:
        plt.loglog(x, y, marker = marker, color = color, linewidth = linewidth,label = label)
    else:
        plt.loglog(x, y, linestyle = linestyle, color = color,linewidth = linewidth,label = label)

    # plt.semilogx(x, y, linestyle = linestyle, color = color,lw= linewidth)
    plt.ylabel(ylabel,fontsize=fontsize)#, weight='bold')
    plt.xlabel(xlabel,fontsize=fontsize)#, weight='bold')
    plt.title(title, fontsize=fontsize)#, weight='bold')
    # ax1.get_yaxis().set_label_coords(-0.11,0.5)
    plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
    ax1.tick_params(axis='x', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
    ax1.tick_params(axis='y', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='y', which='minor', direction='in', length=5, width =1)
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.show()
    if saveFig == True:
        # fName = day+'r'+str(shot)+'_parameters.png'
        fig.savefig(fName,dpi=600,facecolor='w',edgecolor='k')

def semilogxplot(fignum, x, y, linestyle, xlabel, ylabel, title, label, xlim, ylim,
                color, fontsize, linewidth, fName, marker, scatter = True, horzline = True, saveFig=False):
    '''
    fignum, x, y, xlabel, ylabel, title, xlim, ylim, color, fontsize, fName, horzline = True, saveFig=False
    '''
    fig=plt.figure(num = fignum,figsize=(8.5,6),facecolor='w',edgecolor='k')#, dpi=600)
    fig.subplots_adjust(top=0.94, bottom=0.11, left = 0.14, right=0.97, hspace=0)

    ax1=plt.subplot(1,1,1)
    # plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)
    if horzline ==True:
        ax1.axhline(linewidth=1, color='k', linestyle = 'dotted')
    if scatter ==True:
        plt.semilogx(x, y, marker = marker, color = color, linewidth = linewidth,label = label)
    else:
        plt.semilogx(x, y, linestyle = linestyle, color = color,linewidth = linewidth,label = label)

    # plt.semilogx(x, y, linestyle = linestyle, color = color,lw= linewidth)
    plt.ylabel(ylabel,fontsize=fontsize)#, weight='bold')
    plt.xlabel(xlabel,fontsize=fontsize)#, weight='bold')
    plt.title(title, fontsize=fontsize)#, weight='bold')
    # ax1.get_yaxis().set_label_coords(-0.11,0.5)
    plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
    ax1.tick_params(axis='x', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
    ax1.tick_params(axis='y', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='y', which='minor', direction='in', length=5, width =1)
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.show()
    if saveFig == True:
        # fName = day+'r'+str(shot)+'_parameters.png'
        fig.savefig(fName,dpi=600,facecolor='w',edgecolor='k')

def semilogyplot(fignum, x, y, linestyle, xlabel, ylabel, title, label, xlim, ylim,
                color, fontsize, linewidth, fName, marker, scatter = True, horzline = True, saveFig=False):
    '''
    fignum, x, y, xlabel, ylabel, title, xlim, ylim, color, fontsize, fName, horzline = True, saveFig=False
    '''
    fig=plt.figure(num = fignum,figsize=(8.5,6),facecolor='w',edgecolor='k')#, dpi=600)
    fig.subplots_adjust(top=0.94, bottom=0.11, left = 0.14, right=0.97, hspace=0)

    ax1=plt.subplot(1,1,1)
    # plt.text(0.07,0.92,'(a)',fontsize=26, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)
    if horzline ==True:
        ax1.axhline(linewidth=1, color='k', linestyle = 'dotted')
    if scatter ==True:
        plt.semilogy(x, y, marker = marker, color = color, linewidth = linewidth,label = label)
    else:
        plt.semilogy(x, y, linestyle = linestyle, color = color,linewidth = linewidth,label = label)

    # plt.semilogx(x, y, linestyle = linestyle, color = color,lw= linewidth)
    plt.ylabel(ylabel,fontsize=fontsize)#, weight='bold')
    plt.xlabel(xlabel,fontsize=fontsize)#, weight='bold')
    plt.title(title, fontsize=fontsize)#, weight='bold')
    # ax1.get_yaxis().set_label_coords(-0.11,0.5)
    plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
    ax1.tick_params(axis='x', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
    ax1.tick_params(axis='y', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='y', which='minor', direction='in', length=5, width =1)
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.show()
    if saveFig == True:
        # fName = day+'r'+str(shot)+'_parameters.png'
        fig.savefig(fName,dpi=600,facecolor='w',edgecolor='k')

def doubleYplot(fignum, x, y1, y2, linestyle1, linestyle2, xlabel, ylabel1, ylabel2, title, xlim, ylim1, ylim2,
                color1, color2, fontsize, linewidth, fName, marker, scatter=True, horzline = False, saveFig=False):
    '''
    for example:
    tt.doubleYplot(1, a[:,0], a[:,1],a[:,2],'-','-','Stuffing flux (mWb)','C$(\delta n, \delta B_\parallel)$','Time ($\mu s$)',
    'Correlation coefficient and occurance time wrt stuffing flux', (0.7,1.3),(-0.7, -1),(43,55),'blue','red',18,2,'','h')

    '''
    fig=plt.figure(num = fignum,figsize=(8.5,6),facecolor='w',edgecolor='k')#, dpi=600)
    fig.subplots_adjust(top=0.94, bottom=0.11, left = 0.14, right=0.87, hspace=0)
    # print 'done'
    ax=plt.subplot(1,1,1)
    if horzline ==True:
        ax.axhline(linewidth=1, color='k', linestyle = 'dotted')
    # plt.plot(x, y1, color = color1,lw= 2,linestyle = linestyle2)
    if scatter ==True:
        plt.scatter(x, y1, marker = marker, s = 50*linewidth,color = color1,lw= linewidth)
    plt.plot(x, y1, linestyle = linestyle1, color = color1,lw= linewidth)
    # print 'done'
    plt.ylabel(ylabel1,fontsize=fontsize, color = color1)#, weight='bold')
    plt.xlabel(xlabel,fontsize=fontsize)#, weight='bold')
    plt.title(title, fontsize=fontsize)#, weight='bold')
    # ax1.get_yaxis().set_label_coords(-0.11,0.5)
    plt.setp(ax.spines.values(), linewidth=2)#changing the axis linewidth
    ax.tick_params(axis='x', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
    ax.tick_params(axis='y', direction='in', length=7, width =2, labelsize = fontsize-2, colors = color1)
    ax.tick_params(axis='y', which='minor', direction='in', length=5, width =1)#, colors = color1)
    ax.spines['left'].set_color(color1)
    plt.ylim(ylim1)
    plt.xlim(xlim)

    ax1=ax.twinx()#make second y axis (while duplicating/twin x)
    if scatter ==True:
        ax1.scatter(x, y2, marker = marker, s = 50*linewidth, color = color2,lw= linewidth)
    ax1.plot(x, y2, linestyle = linestyle2, color = color2,lw= linewidth)
    plt.ylim(ylim2)
    ax1.set_ylabel(ylabel2,fontsize=fontsize, color = color2)#, weight='bold')
    ax1.tick_params(axis='y', direction='in', length=7, width =2, labelsize = fontsize-2, colors = color2)
    ax1.tick_params(axis='y', which='minor', direction='in', length=5, width =1)#, colors = color2)
    ax1.spines['right'].set_color(color2)

    plt.show()
    if saveFig == True:
        # fName = day+'r'+str(shot)+'_parameters.png'
        fig.savefig(fName,dpi=600,facecolor='w',edgecolor='k')

def lineplot(fignum, x, y, linestyle, xlabel, ylabel, title, label, xlim, ylim, color,
            fontsize, linewidth, fName, marker, scatter = True, horzline = True, saveFig=False):
    '''
    fignum, x, y, xlabel, ylabel, title, xlim, ylim, color, fontsize, fName, horzline = True, saveFig=False
    '''
    fig=plt.figure(num = fignum,figsize=(8.5,6),facecolor='w',edgecolor='k')#, dpi=600)
    fig.subplots_adjust(top=0.94, bottom=0.11, left = 0.14, right=0.97, hspace=0)

    ax1=plt.subplot(1,1,1)
    if horzline ==True:
        ax1.axhline(linewidth=1, color='k', linestyle = 'dotted')
    if scatter ==True:
        plt.scatter(x, y, marker = marker, s = 50*linewidth, color = color,lw= linewidth, label = label)
    else:
        plt.plot(x, y, linestyle = linestyle, color = color,lw= linewidth, label = label)

    plt.ylabel(ylabel,fontsize=fontsize)#, weight='bold')
    plt.xlabel(xlabel,fontsize=fontsize)#, weight='bold')
    plt.title(title, fontsize=fontsize)#, weight='bold')
    # ax1.get_yaxis().set_label_coords(-0.11,0.5)
    plt.setp(ax1.spines.values(), linewidth=2)#changing the axis linewidth
    ax1.tick_params(axis='x', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='x', which='minor', direction='in', length=5, width =1)
    ax1.tick_params(axis='y', direction='in', length=7, width =2, labelsize = fontsize-2)
    ax1.tick_params(axis='y', which='minor', direction='in', length=5, width =1)
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.show()
    if saveFig == True:
        # fName = day+'r'+str(shot)+'_parameters.png'
        fig.savefig(fName,dpi=600,facecolor='w',edgecolor='k')


def straight_line(time, A, B): # this is your 'straight line' y=f(x)
    return A*time + B

def tanh_func(A, B):
    return A*np.tanh(B)

def three_moving_avg(X, Y, Z, N):
    '''
    X, Y, Z = lists to be averaged over a window of length N
    N = 125 means 1 mus
    Till now this is the best among other slding average tools
    len(result) = len(y) : which is the best feature
    '''
    if N==0:
        result1 = X
        result2 = Y
        result3 = Z
    else:
        sum1 = 0
        sum2 = 0
        sum3 = 0
        result1 = list( 0 for x in X)
        result2 = list( 0 for y in Y)
        result3 = list( 0 for z in Z)
        for i in range( 0, N ):
            sum1 = sum1 + X[i]
            result1[i] = sum1 / (i+1)
            sum2 = sum2 + Y[i]
            result2[i] = sum2 / (i+1)
            sum3 = sum3 + Z[i]
            result3[i] = sum3 / (i+1)
        for i in range( N, len(X) ):
            sum1 = sum1 - X[i-N] + X[i]
            result1[i] = sum1 / N
            sum2 = sum2 - Y[i-N] + Y[i]
            result2[i] = sum2 / N
            sum3 = sum3 - Z[i-N] + Z[i]
            result3[i] = sum3 / N
    return result1, result2, result3

def single_moving_avg(X, N):
    '''
    X, Y, Z = lists to be averaged over a window of length N
    N = 125 means 1 mus
    Till now this is the best among other slding average tools
    len(result) = len(y) : which is the best feature
    '''
    if N==0:
        result1 = X
    else:
        sum1 = 0
        result1 = list( 0 for x in X)
        for i in range( 0, N ):
            sum1 = sum1 + X[i]
            result1[i] = sum1 / (i+1)
        for i in range( N, len(X) ):
            sum1 = sum1 - X[i-N] + X[i]
            result1[i] = sum1 / N
    return result1


def center_point_avg(x, y, z, N):
    '''
    N should be an even integer
    '''
    if N==0:
        xm = x
        ym = y
        zm = z
    else:
        xm = np.zeros((len(x),))
        ym = np.zeros((len(x),))
        zm = np.zeros((len(x),))
        for i in range(len(x)):
            # xm[i] = np.sum(x[(i-N/2):(i+N/2)])/N
            # ym[i] = np.sum(y[(i-N/2):(i+N/2)])/N
            # zm[i] = np.sum(z[(i-N/2):(i+N/2)])/N
            xm[i] = np.mean(x[(i-N/2):(i+N/2)])
            ym[i] = np.mean(y[(i-N/2):(i+N/2)])
            zm[i] = np.mean(z[(i-N/2):(i+N/2)])
    return xm, ym, zm

#### two stage rotation #####
def rotation_moving_avg(x, y, z, xm, ym, zm, N):
    ######### Anticlockwise along z axis #############
    theta = np.arccos(xm/np.sqrt(xm**2+ym**2))
    theta[np.where(ym>0)]*=-1
    Bx_prime = np.multiply(x,np.cos(theta))-np.multiply(y,np.sin(theta))
    By_prime = np.multiply(x,np.sin(theta))+np.multiply(y,np.cos(theta))
    Bz_prime = z
    # # # ########### Clockwise along y axis #########
    # gamma = np.arctan2(Bz_prime, Bx_prime)-np.pi/2
    Bx_primem, By_primem, Bz_primem = three_moving_avg(Bx_prime, By_prime, Bz_prime, N)
    gamma = np.arctan2(zm, Bx_primem)-np.pi/2
    Bx_prime2 = np.multiply(Bz_prime,np.sin(gamma))+np.multiply(Bx_prime,np.cos(gamma))
    By_prime2 = By_prime
    Bz_prime2 = np.multiply(Bz_prime,np.cos(gamma))-np.multiply(Bx_prime,np.sin(gamma))

    return Bx_prime2, By_prime2, Bz_prime2

#### two stage rotation #####
def two_stage_rotation(x, y, z):
    ######### Anticlockwise along z axis #############
    theta = np.arccos(np.mean(x)/np.sqrt(np.mean(x)**2+np.mean(y)**2))
    if np.mean(y)>0:
        theta=-theta
    Bx_prime = np.multiply(x,np.cos(theta))-np.multiply(y,np.sin(theta))
    By_prime = np.multiply(x,np.sin(theta))+np.multiply(y,np.cos(theta))
    Bz_prime = z
    ########### Clockwise along y axis #########
    gamma = np.arctan2(np.mean(Bz_prime), np.mean(Bx_prime))-np.pi/2
    Bx_prime2 = np.multiply(Bz_prime,np.sin(gamma))+np.multiply(Bx_prime,np.cos(gamma))
    By_prime2 = By_prime
    Bz_prime2 = np.multiply(Bz_prime,np.cos(gamma))-np.multiply(Bx_prime,np.sin(gamma))

    return Bx_prime2, By_prime2, Bz_prime2

############### Rotation matrix ##############3
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    Example:
    vector = [3, 5, 0] # vector to be rotated
    axis = [4, 4, 1] # usually a unit vector
    theta = 1.2 ## in radians
    rot_matrix = rotation_matrix(axis,theta) # rotation matrix
    rotated_vector = np.dot(rot_matrix, vector)
    # [ 2.74911638  4.77180932  1.91629719]
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
def using_rot_mat(x, y, z):
    magB = np.sqrt(x**2+y**2+z**2)
    axis = [np.mean(x)/np.mean(magB), np.mean(y)/np.mean(magB), np.mean(z)/np.mean(magB)]
    theta = math.acos(np.mean(z)/np.mean(magB))
    rot_matrix = rotation_matrix(axis, theta/2)
    vector = np.zeros([len(x), 3])
    for i in range(len(x)):
        vector[i] = np.dot(rot_matrix, [x[i], y[i], z[i]])
    return vector[:,0], vector[:,1], vector[:,2]

def getcorr(a,b, mean = True):
    '''If already using a highpass filtered signal then don't use mean
    subtraction as highpass filter gives you high frequency fluctuations
    '''
    if mean:
        x = a - np.mean(a)
        y = b - np.mean(b)
    else:
        x = a
        y = b
    x_sd = (np.mean([n**2 for n in a]))**0.5
    y_sd = (np.mean([n**2 for n in b]))**0.5
    numerator = np.dot(x/x_sd,y/y_sd)/len(x)
    # x_sd = (np.mean([n**2 for n in a]))**0.5
    # y_sd = (np.mean([n**2 for n in b]))**0.5
    # return numerator/(x_sd*y_sd)
    return numerator

###### my filter #####
def butter_filter(array, cutoff, fs, order, btype):
    '''
    cutoff and fs (sampling frequency) are both in MHz
    for picoscope, fs = 125 (in MHz)
    '''
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype, analog=False)
    y_filtered = signal.lfilter(b, a, array)
    return y_filtered

def butterworth(interval, cutoff, fs = 65E6, order=5, btype = 'high'):
    #fs = 1/abs(interval[2]-interval[1])  #sampling_rate

    nyq = fs * 0.5

    stopfreq = float(cutoff)
    cornerfreq = 0.4 * stopfreq  # (?)

    ws = cornerfreq/nyq
    wp = stopfreq/nyq

    # for bandpass:
    # wp = [0.2, 0.5], ws = [0.1, 0.6]

    N, wn = scipy.signal.buttord(wp, ws, 3, 16)   # (?)

    # for hardcoded order:
    N = order

    b, a = scipy.signal.butter(N, wn, btype)   # should 'high' be here for bandpass?
    sf = scipy.signal.filtfilt(b, a, interval)
    return sf


###########################LOW-LEVEL FUNCTIONS#############################
def cust_butter_filter(cutoff, fs, btype='high',order=1):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs = 1e8, order=1, btype='high'):
    # fs= 65E6# 1/abs(data[2]-data[1])
    b, a = cust_butter_filter(cutoff, fs, order=order, btype='high')
    y_filtered = signal.lfilter(b, a, data)
    return y_filtered

def in_limits(arr,limits,return_bool=False,
              return_values=False):
    """
    Returns indices where ARR lies within range
    set by LIMITS.  Optionally return ARR[indices]
    or T/F array of length ARR.size

    ARR:
        1D array.
    LIMITS:
        2-element list of form [arrmin,arrmax]
        corresponding to desired range of ARR.
    RETURN_VALUES = False:
        Boolean controlling whether indices or
        ARR[indices] is returned.  Indices only
        returned by default.
    RETURN_BOOL = False:
        Boolean controlling whether condition
        returned instead of indices.  Indices only
        returned by default.

    RETURNS:
       INDICES (tuple of arrays - must be tuple for arr[indices]
       to function properly with ND arrays) or ARR[INDICES]
       (1D array with length dependent on LIMITS).
    """
    condition = np.logical_and(arr >= limits[0],arr <= limits[1])
    if return_bool:
        return condition
    indices = np.where(condition)
    if return_values:
        return arr[indices]
    else:
        return indices

#def get_corr(t,sig1,sig2,mode='same',normalized=1,optimize=1):
#    """
#    lag,corr = get_corr(t,sig1,sig2,mode='same',normalized=1,):
#    NORMALIZED uses (N-1)*sigma1*sigma2 to normalize
#        correlation function (standard deviation normalization)
#    OPTIMIZE calls OPTLENGTH on both signals so that the correlation
#        runs quickly.  NOTE: Do NOT run without this on raw probe
#        signals.  The number of points is absurdly un-optimized (odd,
#        maybe close to prime, etc) and the traces are huge (1-2 Msamples).
#
#    """
#
#    #n = len(sig1)
#    #correlate
#    corr = np.correlate(sig1,sig2,mode=mode)
#    #normalize
#    if normalized:
#        corr /= (len(t) - 1)*sig1.std()*sig2.std()
#    #    corr /= (1.0/(n-1))*sig1.std()*sig2.std()
#    #calculate lag array
#    dt = t[1]-t[0]
#    tau = dt*(np.arange(corr.size) - corr.size/2)
#    #integer division makes this agree with correlate
#    #correlate leaves off the most positive lag value if the number of points is even
#    return tau,corr


def get_gaussian_moving_avg(x,y,xsigma,**kwargs):
    """
    Wrapper for scipy.ndimage.gaussian_filter to calculate moving average
    with a Gaussian kernel/window.

    avg_signal = get_gaussian_moving_avg(signal,sigma)

    Parameters
    ----------
    x : array_like
        Evenly-spaced array corresponding to sample locations/times

    y : array_like
        Input array of samples to smooth/average.

    xsigma : int or list of ints
        Width parameter of Gaussian window, in units of `x`

    kwargs : various
        Parameters accepted by gaussian_filter1d.  Defaults
        seem appropriate for this application.

    Returns
    ----------
    avg_signal : mimics input `signal` type and shape

    """

    #figure out how many points correspond to xsigma:
    dx = x[1] - x[0]  # sample spacing
    sigma = int(xsigma/dx + 0.5) # round without calling np.round

    #print width of Gaussian in time/frequency
    width_factor = np.sqrt(2*np.log(2)) # convert sigma to HWHM
    #print "HWHM is %4.2e [x] or %d [points]."%(width_factor*xsigma,
    #                                           width_factor*sigma)
    #print "Corresponding HWHM in frequency space is %4.2e [1/x]"%(width_factor*1./(2*np.pi*xsigma))

    # call ndimage.gaussian_filter1d
    return gaussian_filter1d(y,sigma,**kwargs)


def find_closest(arr,val,value=True):
    """
    Returns value or index of ARR
    closest to VAL.
    """
    if value:
        return arr[absmin(arr-val,value=False)]
    else:
        return absmin(arr-val,value=False)



def tindex(timearr,timevalue,delta=2e-2):
    tind = np.where(np.absolute((timearr)-(timevalue))<delta)
    tind = tind[0][0]
    return tind

def tindex_center(timearr,timevalue):
    tind = np.argmin(abs((timearr)-(timevalue)))
    return tind

def tindex_min(timearr,timevalue):
    minval = min(abs((timearr)-timevalue))
    tind = np.where(abs((timearr)-(timevalue)) == minval)
    tind = tind[0][0]
    return tind

def tindex_near(timearr,timevalue,threshold):
    tinds = np.where(abs(timearr-timevalue) < threshold)
    return tinds

def firstzero(timearr,timevalue):
    for n in range(len(timearr)):
        tind = n
        if timearr[n] < 0: break
    return tind


def classic_lstsqr(x_list, y_list):
    """ Computes the least-squares solution to a linear matrix equation. """
    N = len(x_list)
    x_avg = sum(x_list)/N
    y_avg = sum(y_list)/N
    var_x, cov_xy = 0, 0
    for x,y in zip(x_list, y_list):
        temp = x - x_avg
        var_x += temp**2
        cov_xy += temp * (y - y_avg)
    slope = cov_xy / var_x
    y_interc = y_avg - slope*x_avg
    return (slope, y_interc)

def matrix_lstsqr(x, y):
    """ Computes the least-squares solution to a linear matrix equation. """
    X = np.vstack([x, np.ones(len(x))]).T
    return (np.linalg.inv(X.T.dot(X)).dot(X.T)).dot(y)

def Linear_fit_slope(x, m):
    """Computes the slope for a fixed intercept"""
    c=74
    return m*x+c


def fix_time_bounds(tshort,tlong):
    indices1 = [i for i,v in enumerate(tlong >= tshort[0]) if v]
    indices2 = [i for i,v in enumerate(tlong <= tshort[-1]) if v]
    indices = np.intersect1d(indices1,indices2)
    return indices, tlong[indices]

def fix_bounds(t_bounds, t_array):
    indices1 = [i for i, v in enumerate(t_array >= t_bounds[0]) if v]
    indices2 = [i for i, v in enumerate(t_array <= t_bounds[-1]) if v]
    indices = np.intersect1d(indices1, indices2)
    return indices, t_array[indices]

def fix_array(t_bounds, t_array, vals):
    indices1 = [i for i, v in enumerate(t_array >= t_bounds[0]) if v]
    indices2 = [i for i, v in enumerate(t_array <= t_bounds[-1]) if v]
    indices = np.intersect1d(indices1, indices2)
    if(len(vals) == len(t_array)):
        return indices, t_array[indices], vals[indices]
    else:
        print("error")
        return indices, t_array[indices], vals


def local_min_before(t_bound, t_array, vals):
    #1. find the index of the time bound
    #2. scan backwards though the array from the time array
    #   until you find somthing that is a local min.
    if(t_array[0]> t_bound):
        print("Error, bound is too low")
        return vals[0]
    indices1 = [i for i, v in enumerate(t_array >= t_array[0]) if v]
    indices2 = [i for i, v in enumerate(t_array <= np.array(t_bound)) if v]
    indices = np.intersect1d(indices1, indices2)
    reversed = indices[::-1]
    vals = vals[reversed]
    min = vals[0]
    for i in range(len(vals)-1):
        if vals[i] < min:
            min =vals[i]
        if vals[i] < vals[i+1]:
            return min
    return min




def ipybreak():
    """
    Retained only for backwards compatibility.
    Ctrl-C now interrupts even inside running subfunction.

    ipython 0.11 seems to have fixed the problem,
    as well as eliminating the wthread keyword that
    was the likely cause.
    """
    pass

def old_ipybreak():
    import IPython.Shell
    if IPython.Shell.KBINT:
        IPython.Shell.KBINT = False
        raise SystemExit

#http://code.activestate.com/recipes/134892/
class GetChUnix:
    def __init__(self):
        import tty, sys

    def __call__(self):
        import sys, tty, termios
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch

def getch():
    g = GetChUnix()
    return g()

###########################GRAPHICS#############################
def data2axes(xy,ax):
    """
    Converts a pair xy=(x,y) of data coordinates to
    axes coordinates.
    """
    disp2axes = ax.transAxes.inverted().transform
    data2disp = ax.transData.transform
    return disp2axes(data2disp(xy))

def axes2data(xy,ax):
    """
    Converts a pair xy=(x,y) of axes coordinates to
    data coordinates.
    """
    disp2data = ax.transData.inverted().transform
    axes2disp = ax.transAxes.transform
    return disp2data(axes2disp(xy))

def fig2data(xy,fig,ax):
    """
    Converts a pair xy=(x,y) of axes coordinates to
    data coordinates.
    """
    disp2data = ax.transData.inverted().transform
    fig2disp = fig.transFigure.transform
    return disp2data(fig2disp(xy))

def data2fig(xy,fig,ax):
    """
    Converts a pair xy=(x,y) of axes coordinates to
    data coordinates.
    """
    disp2fig = fig.transFigure.inverted().transform
    data2disp = ax.transData.transform
    return disp2fig(data2disp(xy))

def mygreen():
    """
    return (0.0,0.8,0.0)
    """
    return (0.0,0.8,0.0)

def mycolortriplet():
    """
    return ((0,0,0),(0,0,1),mygreen())
    """
    return ((0,0,0),(0,0,1),mygreen())

def mycmap():
    """
    Returns colormap that I'm using lately.
    (convenience mapping function)
    """
    #return yanceyct()
    return two_color_linearL_quadS()

def PGcmap():
    """
    Returns ListedColormap of green to black to pink

    Good for fluctuations (where
    you want to ignore the mean)
    """
    x = numpy.arange(256)/255.
    H = 0.2*numpy.tanh(4*numpy.pi*(x-0.5)) + 0.1
    H = H[::-1]
    #L = numpy.ones(256)*0.5
    #S = 2*(x-0.5)**2 +0.5
    #S = 1.5*numpy.abs(x-0.5) + 0.25
    S = numpy.ones(256)*0.5
    L = 1.5*numpy.abs(x-0.5) + 0.1
    rgb = numpy.zeros((256,3),float)
    for ii in numpy.arange(256):
        rgb[ii,:] = colorsys.hls_to_rgb(H[ii],L[ii],S[ii])

    return matplotlib.colors.ListedColormap(rgb,name='test')


def hvis(l=0.5,s=1.):
    H = numpy.linspace(0,1,num=100,endpoint=True)
    L = numpy.ones(100)*l
    S = numpy.ones(100)*s
    rgb = numpy.zeros((100,3),float)
    for ii in numpy.arange(100):
        rgb[ii,:] = colorsys.hls_to_rgb(H[ii],L[ii],S[ii])

    return matplotlib.colors.ListedColormap(rgb,name='test')


def two_color_flatS(H1=0.7,H2=0.9,S=0.75):
    """
    Returns a 256-color ListedColormap object.

    H1 is the hue applied to the lower half of the elements
    H2 is the hue applied to the upper half of the elements
    S is the saturation value for all elements

    The hue function uses a tanh to transition
    smoothly from H1 to H2 in the center of the map.
    Default values give a colormap of indigo
    to black to pink.

    Then lightness is linear [abs(x)], with a
    minimum of 0.1 in the middle of the map
    and a maximum of 0.9 at the extrema of the map.

    Good for fluctuation patterns (modes)
    (where you want to ignore the mean)
    """
    x = numpy.arange(256)/255.
    H =(H2-H1)*(0.5*numpy.tanh(4*numpy.pi*(x-0.5)) + 0.5) + H1
    S = numpy.ones(256)*S
    # following gives linear ramp from 0.1 to 0.9
    L = 1.6*numpy.abs(x-0.5) + 0.1
    rgb = numpy.zeros((256,3),float)
    for ii in numpy.arange(256):
        rgb[ii,:] = colorsys.hls_to_rgb(H[ii],L[ii],S[ii])

    return matplotlib.colors.ListedColormap(rgb,name='test')


def two_color_linearL_quadS(H1=0.6,H2=0.13,dS=0.75,dL=0.7):
    """
    Returns a 256-color ListedColormap object.

    H1 is the hue applied to the lower half of the elements
    H2 is the hue applied to the upper half of the elements
    dS is the range of saturation values, centered on 0.5
    dL is the range of lightness values, centered on 0.5

    The hue steps from H1 to H2 in the middle of the map.
    The saturation is quadratic [x**2], with a minimum in the
    middle of the map.
    Then lightness is linear [abs(x)], with a minimum in the middle
    of the map.


    Recall that H,L,S are all in the range [0,1].
    """
    x = numpy.arange(256)/255.
    H = numpy.zeros(256)
    H[:128] = H1
    H[128:] = H2
    S =  dS*4.0*(x-0.5)**2 + (0.5 - dS/2.)
    L = dL*2.0*numpy.abs(x-0.5) + (0.5 - dL/2.)
    rgb = numpy.zeros((256,3),float)
    for ii in numpy.arange(256):
        rgb[ii,:] = colorsys.hls_to_rgb(H[ii],L[ii],S[ii])

    return matplotlib.colors.ListedColormap(rgb,name='test')


def single_color(H=0.66,end_flatness=10,S=0.5):
    """
    Returns ListedColormap of a single hue
    (~linear lightness scale)

    H [0-1] specifies hue
    END_FLATNESS specifies the number of points at the ends
    used to transition to white and black.  If zero or > 127,
    a linear scale in full luminance is used.  Otherwise,
    a linear scale using the middle 80% of the luminance values
    is constructed and the first and last END_FLATNESS points are
    replaced by steeper lines leading from 0 and 1 to the original
    line at the appropriate points.

    In essence, this scale extends the range where the middle
    80% luminance values are used and reduces the range where
    the first and last 10% are used.
    """
    n = end_flatness
    x = numpy.arange(256)/255.
    H = numpy.ones(256)*H
    S = numpy.ones(256)*S
    if end_flatness > 127:
        L = x
    else:
        L = x*0.8 + 0.1
        if n:
            L[0:n] = numpy.linspace(0,L[n],num=n,endpoint=True)
            L[-n:] = numpy.linspace(L[-n],1,num=n,endpoint=True)
    rgb = numpy.zeros((256,3),float)
    for ii in numpy.arange(256):
        rgb[ii,:] = colorsys.hls_to_rgb(H[ii],L[ii],S[ii])

    return matplotlib.colors.ListedColormap(rgb,name='test')


def yanceyct():
    """Returns color map object corresponding to Yanceys 256-color
    black-green-yellow-red-white luminosity map."""
    ii_step = numpy.arange(256)/255.0
    H = numpy.mod(ii_step[::-1]*300,360)
    H[0:151] = (numpy.arange(151)*(60.0/150.0) + H[151])[::-1]
    H /= 360.  #python HLS has H in [0,1]
    L = (0.7*ii_step**(1.2) + 0.1) + 0.6*(numpy.exp(
            -numpy.log(2.0)*(ii_step/.95)**30) - 1.0)
    S = 0.5*ii_step + 0.5
    S[0] = 0
    L[0] = 0
    S[255] = 1
    L[255] = 1

    rgb = numpy.zeros((256,3),float)

    for ii in numpy.arange(256):
        rgb[ii,:] = colorsys.hls_to_rgb(H[ii],L[ii],S[ii])

    return matplotlib.colors.ListedColormap(rgb,name='yancey')


def myspectralr():

    H = numpy.array([ 0.69678715,  0.56354916,  0.44746377,  0.3128655 ,  0.19354839,
                      0.16666667,  0.12318841,  0.08226496,  0.03954802,  0.9812362 ,
                      0.93099788])
    L = numpy.array([ 0.47254904,  0.46862746,  0.58039217,  0.75490198,  0.77843139,
                      0.87450981,  0.77058825,  0.68627451,  0.60980393,  0.53921569,
                      0.31176472])
    S = numpy.array([ 0.34439834,  0.58158996,  0.42990656,  0.456     ,  0.82300885,
                      1.        ,  0.98290598,  0.975     ,  0.88944724,  0.64255321,
                      0.98742138])

    H = numpy.array([ 0.69678715,  0.56354916,  0.44746377,  0.3128655 ,  0.19354839,
                      0.18666667,  0.17318841,  0.12226496,  0.054802,  0.9812362 ,
                      0.93099788])
    L = numpy.array([ 0.47254904,  0.46862746,  0.58039217,  0.75490198,  0.843139,
                      0.97450981,  0.77058825,  0.68627451,  0.60980393,  0.53921569,
                      0.31176472])
    #S = numpy.array([ 0.34439834,  0.58158996,  0.42990656,  0.456     ,  0.82300885,
    #                  .9        ,  0.89290598,  0.875     ,  0.78944724,  0.64255321,
    #                  0.98742138])

    x = numpy.arange(11)/10.
    S *= numpy.abs(x-0.5)**0.5 + 0.2
    #S *= (x-0.5)**2+0.5
    #S *= (4*(x-0.5)**2)
    #S[4:7] = 1.2-S[4:7]
    #H[2:10] = 2*(numpy.linspace(H[2],H[10],num=8,endpoint=True))**2

    red = []
    green = []
    blue = []
    for ii in numpy.arange(11):
        r,g,b = colorsys.hls_to_rgb(H[ii],L[ii],S[ii])
        red.append((float(ii)/10.,r,r))
        green.append((float(ii)/10.,g,g))
        blue.append((float(ii)/10.,b,b))

    cdict = {"red":red,"green":green,"blue":blue}
    return matplotlib.colors.LinearSegmentedColormap("test",cdict)


def cvals2D(hl_arr,s_arr,cmap=mycmap()):
    """
    Returns color values corresponding to the
    combination of two arrays, one for color/lightness
    values and one for saturation.
    """
    colorarr = cmap(hl_arr)
    hcolorarr = matplotlib.colors.rgb_to_hsv(colorarr[...,:-1])
    #hcolorarr[...,1] = s_arr
    hcolorarr[...,2] = s_arr
    colorarr = matplotlib.colors.hsv_to_rgb(hcolorarr)
    return colorarr


def crosshair(ax,**kwargs):
    """
    Just draws a pair of lines through zero
    with the properties given by **kwargs.
    """
    if "linestyle" not in kwargs.keys():
        kwargs["linestyle"] = '--'
    if "color" not in kwargs.keys():
        kwargs["color"] = 'k'
    ax.plot([0,0],ax.get_ylim(),**kwargs)
    ax.plot(ax.get_xlim(),[0,0],**kwargs)
    return


def testmatrix(nx=100):
    """Returns diagonal matrix with elements equal to the row and column
    number for testing image orientation"""
    #return numpy.diag(numpy.arange(nx))
    a = numpy.arange(nx*nx)
    a = a.reshape((nx,nx))
    return a


def scatterimage(image,cmap=None,x=None,y=None,
                 aspect='equal'):
    """
    Use pylab.scatter to plot an array as an image.

    Allows setting the parameters of individual pixels,
    including color, alpha, etc.

    IMAGE can be a float or int array, as well as an
    array of RGB(A) values.
    """

    if image.ndim == 2:
        ny,nx = image.shape
        nc = 0
    elif image.ndim == 3:
        ny,nx,nc = image.shape
        if nc == 4:
            print ("Image read as (nx,nx,RGBA) array.")
            image = image.reshape(nx*ny,4)
        if nc == 3:
            print ("Image read as (nx,nx,RGB) array.")
            image = image.reshape(nx*ny,3)
    else:
        print ("Image dimensions must be (nx,ny), (nx,nx,RGB), or (nx,nx,RGBA)")
        return 0

    if x is None:
        x = numpy.arange(nx)
    if y is None:
        y = numpy.arange(ny)

    x,y = pylab.meshgrid(x,y)

    if cmap is None:
        cmap = cm.gist_heat

    ax = new_plot()
    ax.set_xlim(x.min(),x.max())
    ax.set_ylim(y.min(),y.max())
    ax.set_aspect(aspect)
    dw,dh = ax.transData.transform((x.max()-x.min(),y.max()-y.min()))
    dw /= float(nx-1)
    dh /= float(ny-1)
    m = ([(0,0),(dw,0),(dw,dh),(0,dh)],0) # see scatter docstring
    ss = ax.scatter(x,y,s=dw*dh*1.2,c=image,marker=m,cmap=cmap,edgecolors='none')
    ax.set_xlim(x.min(),x.max())
    ax.set_ylim(y.min(),y.max())
    pylab.show()
    return ss


def imview(image,x=None,xmin=None,xmax=None,y=None,ymin=None,ymax=None,
           ratio=0.8,cmap=None,oplot=0,ax=None,has_colorbar=1,
           aspect=None,interpolation='nearest',
           cbkwargs=dict(position="right",size="5%",pad=0.1),
           bare=False,symcb=None,**kwargs):
    """
    ax,cb,im = imview(image,x=None,xmin=None,xmax=None,y=None,ymin=None,ymax=None,
                      ratio=0.8,cmap=mycmap(),oplot=0,ax=None,has_colorbar=1,
                      cbkwargs=dict(position="right",size="5%",pad=0.1),
                      aspect='equal',interpolation='nearest',**kwargs)
    Returns AX,CB,IM; an axis object, a colorbar object, and an
    imshow object.
    Wrapper for imshow that implements x,y axes, proper scaling, and a colorbar.
    Uses origin='lower' to display image with (0,0) pixel in lower left corner.
    Gets tick labels from max/min values.  Sets aspect ratio based on
    xmax-xmin/ymax-ymin"""

    if x is None:
        # make index array spanning number of tiles,
        # center pixels on tile points
        x = numpy.arange(image.shape[1]+1) - 0.5
    if y is None:
        y = numpy.arange(image.shape[0]+1) - 0.5

    if xmin is None:
        xmin = x.min()
    if xmax is None:
        xmax = x.max()
    if ymin is None:
        ymin = y.min()
    if ymax is None:
        ymax = y.max()

    if aspect is None:
        if image.shape[0] == image.shape[1]:
            aspect = 'equal'
        else:
            aspect = 'auto'

    if cmap is None:
        if numpy.sign(image.min())*numpy.sign(image.max()) < 0:
            cmap = mycmap()
            symcb = True
        else:
            cmap = matplotlib.cm.gist_heat
            symcb = False

    if symcb:
        # if symmetric colorbar requested
        # set bounds to +- same value
        if "vmin" in kwargs.keys() and "vmax" in kwargs.keys():
            imax = kwargs["vmax"]
            imin = kwargs["vmin"]
            ext = max([abs(imin),abs(imax)])
            imin = -ext
            imax =  ext
        elif "vmin" in kwargs.keys():
            imin = kwargs["vmin"]
            imax = -imin
            print ("Color max set to -imin")
        elif "vmax" in kwargs.keys():
            imax = kwargs["vmax"]
            imin = -imax
            print( "Color min set to -imax")
        else:
            # otherwise symmetrize by making vmin/vmax
            # spaced at bigger of 2*image.min()
            # or 2*image.max(), centered on the mean
            imax = image.max()
            imin = image.min()
            if abs(imax) > abs(imin):
                imin = -imax
                print ("Color range set to +/- I_max")
            else:
                imax = -imin
                print ("Color range set to +/- I_min")
        print(imin,imax)
        kwargs["vmin"] = imin
        kwargs["vmax"] = imax

    #find max/min indices
    xindx = numpy.where(numpy.logical_and(x>=xmin,x<=xmax))[0]
    ixmin = xindx[0]
    ixmax = xindx[-1]+1
    yindx = numpy.where(numpy.logical_and(y>=ymin,y<=ymax))[0]
    iymin = yindx[0]
    iymax = yindx[-1]+1

    #plot
    if not ax:
        ax = new_plot(oplot=oplot)

    im = ax.imshow(image[iymin:iymax,ixmin:ixmax],origin='lower',
                   aspect=aspect,interpolation=interpolation,
                   cmap=cmap,**kwargs)
    dx = (xmax-xmin)/float(image.shape[1])
    dy = (ymax-ymin)/float(image.shape[0])
    im.set_extent((xmin-dx/2.,xmax+dx/2.,ymin-dy/2.,ymax+dy/2.))
    if has_colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes(**cbkwargs)
        cb = pylab.colorbar(im,cax=cax,cmap=cmap,)
    else:
        cb = None
        if bare:
            ax.set_axis_off()
            fig.subplots_adjust(top=1.0,bottom=0.0,left=0.0,right=1.0)

    pylab.show()
    return ax,cb,im


def imviewxy(x,y,image,**kwargs):
    """
    ax,cb,h = imview(image,x=x,y=y,**kwargs)
    Returns AX,CB,IM; an axis object, a colorbar object, and an imshow object.
    Wrapper for imshow that implements x,y axes, proper scaling, and a colorbar.
    Uses origin='lower' to display image with (0,0) pixel in lower left corner.
    Gets tick labels from max/min values.  Sets aspect ratio based on xmax-xmin/ymax-ymin"""

    return imview(image,x=x,y=y,**kwargs)


def arbitrary_line(slope, intercept):
    """Plot a line from slope and intercept"""
    ax = plt.gca()
    x_vals = np.array(ax.get_xlim())
    y_vals = intercept + slope * x_vals
    pylab.semilogy(x_vals, y_vals, '--')
    ax.set_yscale('log')
    ax.set_xscale('log')


def animate(images,cmap=None,slow=0,
            x=None,xmin=None,xmax=None,
            y=None,ymin=None,ymax=None,
            xlabel='',ylabel='',title='',
            show_fnum=0,show_time=False,
            t=None,t_format="%3.1e",
            t_units='s',
            has_colorbar=False,
            aspect=None,symcb=None,
            save_stills=False,
            play=False,
            save_prefix='/home/light/temp_movies/unk_',
            save_kwargs=dict(bbox_inches='tight'),
            xtpos=None,ytpos=None,
            xfpos=None,yfpos=None,
            textprops=dict(family='monospace',size=20),
            axisprops=dict(),**kwargs):
    """
    Uses the set_data() method of a pylab.imshow()
    object to update the image array quickly.

    PARAMETERS:
    -----------
    images: array, (ny,nx,nframes)
        3D array to animate
    cmap: colormap object
        Colormap to apply.  Default is mytools.mycmap()
    slow: float
        Slowing factor to reduce animation speed.
        Introduces a time.sleep(slow) line into
        the display loop.
    kwargs: dictionary
        Extra keyword arguments are passed directly to
        mytools.imview().
    """

    # Set up figure
    if pylab.get_fignums() == []:
        fig = pylab.figure(figsize=(13,8))
    else:
        fig = pylab.gcf()
        fig.clf()

    if not pylab.isinteractive():
        pylab.ion()

    if x is None:
        # make index array spanning number of tiles,
        # center pixels on tile points
        x = numpy.arange(images.shape[1]+1) - 0.5
    if y is None:
        y = numpy.arange(images.shape[0]+1) - 0.5

    if xmin is None:
        xmin = x.min()
    if xmax is None:
        xmax = x.max()
    if ymin is None:
        ymin = y.min()
    if ymax is None:
        ymax = y.max()

    if aspect is None:
        if images.shape[0] == images.shape[1]:
            aspect = 'equal'
        else:
            aspect = 'auto'

    if cmap is None:
        if numpy.sign(images.min())*numpy.sign(images.max()) < 0:
            cmap = mycmap()
            symcb = True
        else:
            cmap = matplotlib.cm.gist_heat
            symcb = False

    if symcb:
        # if symmetric colorbar requested
        # set bounds to +- same value
        if "vmin" in kwargs.keys() and "vmax" in kwargs.keys():
            imax = kwargs["vmax"]
            imin = kwargs["vmin"]
            ext = max([abs(imin),abs(imax)])
            imin = -ext
            imax =  ext
        elif "vmin" in kwargs.keys():
            imin = kwargs["vmin"]
            imax = -imin
            print ("Color max set to -imin")
        elif "vmax" in kwargs.keys():
            imax = kwargs["vmax"]
            imin = -imax
            print ("Color min set to -imax")
        else:
            # otherwise symmetrize by making vmin/vmax
            # spaced at bigger of 2*image.min()
            # or 2*image.max(), centered on the mean
            imax = images.max()
            imin = images.min()
            if abs(imax) > abs(imin):
                imin = -imax
                print ("Color range set to +/- I_max")
            else:
                imax = -imin
                print ("Color range set to +/- I_min")
        print (imin,imax)
        kwargs["vmin"] = imin
        kwargs["vmax"] = imax

    ax,cb,im = imview(images[...,0],cmap=cmap,aspect=aspect,
                      has_colorbar=has_colorbar,x=x,y=y,
                      **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title,y=1.02)

    if axisprops:
        ax.xaxis.set(**axisprops)
        ax.yaxis.set(**axisprops)

    # Iterate through frames, timing the process
    nframes = images.shape[2]
    if t == None:
        t = numpy.arange(nframes)
    frame = 0
    tstart = time.time()
    textsize = 20
    if show_fnum:
        if xfpos is None:
            xfpos = 1.05*ax.get_xlim()[1]
        if yfpos is None:
            yfpos = 1.05*ax.get_ylim()[1]
        fnum = ax.text(xfpos,yfpos,
                       "Frame %05d"%frame,
                       **textprops)
    if show_time:
        if xtpos is None:
            xtpos = 1.02*ax.get_xlim()[1]
        if ytpos is None:
            ytpos = 1.05*ax.get_ylim()[1]
        ftime = ax.text(xtpos,ytpos,"Time: "+
                        t_format %t[frame]+
                        " ("+t_units+")",
                        **textprops)

    try:
        while True:
            im.set_data(images[...,frame])
            # Add frame number to title if desired
            if show_fnum:
                fnum.set_text("Frame %05d"%frame)
            elif show_time:
                ftime.set_text("Time: "+
                               t_format %t[frame]+
                               " ("+t_units+")")
            pylab.draw()
            pylab.draw()

            # save each frame as png if desired
            if save_stills:
                pylab.savefig(save_prefix+"frame%05d" %frame+".png",
                              **save_kwargs)

            # evaluate next move
            if play:
                # Wait if desired
                time.sleep(slow)
                # iterate sequentially until nframes
                frame += 1
                if frame >= nframes:
                    break

            else:
                if frame == 0:
                    print ("To go forward, press 'f', to reverse, press 'd', "+
                           "to exit, press 'x'\n")

                # step through or exit depending on user input
                next = getch()

                # allow stepping forward or backwards
                if next == 'x':
                    break
                if next == 'f':
                    frame += 1
                if next == 'd':
                    if frame == 0:
                        frame = nframes - 1
                    else:
                        frame -= 1
                if frame >= nframes:
                    frame -= nframes

    except KeyboardInterrupt:
        return

    print ('FPS:' , nframes/(time.time()-tstart))
    if not pylab.isinteractive():
        pylab.ioff()

def animate_plots(x,lines,*morelines,**kwargs):
    """
    Uses the set_data() method of a Line2D
    object to update the array quickly.

    PARAMETERS:
    -----------
    x: array, (nx)
        1D array to plot against
    lines: array, (ny,nt)
        2D array to animate
    morelines: list of arrays, n x (ny,nt)
        other 2D arrays to animate simultaneously
    slow: float
        Slowing factor to reduce animation speed.
        Introduces a time.sleep(slow) line into
        the display loop.
    kwargs: dictionary
        Extra keyword arguments are passed directly to
        mytools.imview().
    """
    if not pylab.isinteractive():
        pylab.ion()

    local_kw = dict(slow=0,pad=0.,
                    xlabel='',ylabel='',title='',
                    show_fnum=0,show_time=False,
                    t=None,t_format="%3.1e",
                    t_units='s',
                    xlim=None,
                    ylim=None,
                    save_stills=False,
                    save_prefix='/home/light/temp_movies/unk_',
                    save_kwargs={},show_legend=False,
                    legkws=dict(),
                    labels=list(numpy.arange(len(morelines)+1).astype('|S1')) )

    # If these keywords are defined, assign them to
    # local variables and remove them from
    # the dictionary so that only plot kwargs remain
    for kw in local_kw.keys():
        if kw in kwargs.keys():
            # If it exists, define kwargs[key] = value as local variable
            if isinstance(kwargs[kw],numpy.ndarray):
                temp = kwargs[kw]
                exec (compile(kw+'= temp','<string>','single'))
            elif isinstance(kwargs[kw],basestring):
                exec (compile(kw+'='+'"'+str(kwargs[kw])+'"','<string>','exec'))
            else:
                exec (compile(kw+'='+str(kwargs[kw]),'<string>','single'))
            # Remove from kwargs dictionary
            kwargs.pop(kw)
        else:
            # Define key local_kw[key] = value as local variable
            if isinstance(local_kw[kw],basestring):
                exec (compile(kw+'='+'"'+str(local_kw[kw])+'"','<string>','exec'))
            else:
                exec (compile(kw+'='+str(local_kw[kw]),'<string>','single'))

    # If there are more lines to plot, count them
    n_extra = len(morelines)

    # Set up original plot
    # Make set ylim() so that all lines are in plot
    if ylim is None:
        ymin = lines.min()
        ymax = lines.max()
    else:
        ymin,ymax = ylim
    ax = new_plot()
    ax.plot(x,lines[:,0],label=labels[0],**kwargs)
    for lindex in range(n_extra):
        ax.plot(x,morelines[lindex][:,0],label=labels[lindex+1],**kwargs)
        if morelines[lindex].min() < ymin:
            ymin = morelines[lindex].min()
        if morelines[lindex].max() > ymax:
            ymax = morelines[lindex].max()
    ax.set_ylim((ymin-pad,ymax+pad))
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title,y=1.02)
    if show_legend:
        ax.legend(**legkws)

    # Iterate through frames, timing the process
    tstart = time.time()
    nframes = lines.shape[1]
    if t == None:
        t = numpy.arange(nframes)
        t_units = 'frames'
    try:
        for frame in range(nframes):
            l = ax.lines[0]
            l.set_data(x,lines[:,frame])
            for lindex in range(n_extra):
                l = ax.lines[lindex+1]
                l.set_data(x,morelines[lindex][:,frame])
            # Add frame number to title if desired
            if show_fnum:
                ax.set_title(title+"    Frame %05d"%frame)
            elif show_time:
                ax.set_title(title+"    Time: "+
                             t_format %t[frame]+
                             " ("+t_units+")")
            pylab.draw()
            # save each frame as png if desired
            if save_stills:
                pylab.savefig(save_prefix+"frame%05d" %frame+".png",
                              **save_kwargs)
            # Wait if desired
            time.sleep(slow)

    except KeyboardInterrupt:
        print("Interrupted by User.")
    print('FPS:', nframes/(time.time()-tstart))
    if not pylab.isinteractive():
        pylab.ioff()

def animate_histogram(data,slow=0,pad=0.,nbins=100,
                      xlabel='',ylabel='',title='',
                      show_fnum=0,show_time=False,
                      normed=True,
                      t=None,t_format="%3.1e",
                      t_units='s',
                      textprops=dict(family='monospace',size=20),
                      xlim=None,
                      ylim=None,
                      xfpos=None,
                      yfpos=None,
                      save_stills=False,
                      save_prefix='/home/light/temp_movies/unk_',
                      save_kwargs={},**kwargs):
    """
    Animates (...,nt) data into histogram,
    with each frame being the histogram of NCOUNT
    points at time t
    """

    nt = data.shape[-1]

    fig = myfig()
    ax = fig.add_subplot(111)
    n,bins,patches = ax.hist(data[...,0].flatten(),bins=nbins,
                                normed=normed,**kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if xlim is None:
        ax.set_xlim(0,data.max())
    else:
        ax.set_xlim(xlim)

    if show_fnum:
        if xfpos is None:
            xfpos = 0.9*ax.get_xlim()[1]
        if yfpos is None:
            yfpos = 1.1*ax.get_ylim()[1]
        fnum = ax.text(xfpos,yfpos,
                       "Frame %05d"%0,
                       **textprops)
    if show_time:
        if xtpos is None:
            xtpos = 0.9*ax.get_xlim()[1]
        if ytpos is None:
            ytpos = 1.1*ax.get_ylim()[1]
        ftime = ax.text(xtpos,ytpos,"Time: "+
                        t_format %t[0]+
                        " ("+t_units+")",
                        **textprops)

    try:
        for ii in range(nt):
            n,bins = numpy.histogram(data[...,ii].flatten(),bins,normed=normed)
            for rect,h in zip(patches,n):
                rect.set_height(h)
            if show_fnum:
                fnum.set_text("Frame %05d"%ii)
            if show_time:
                ftime.set_text("Time: "+
                               t_format %t[ii]+
                               " ("+t_units+")")
            fig.canvas.draw()
    except KeyboardInterrupt:
        print ("Interrupted by User.")


def myfigsize():
    """
    return (13,8)
    """
    return (13,8)

def myfig(newfig=False,oplot=False,mysize=False,**kwargs):
    """
    Generates new figure with myfig size
    or returns figure if it exists.
    """
    if newfig:
        if mysize:
            fig = pylab.figure(figsize=myfigsize(),**kwargs)
        else:
            fig = pylab.figure(**kwargs)
    else:
        if pylab.get_fignums() == []:
            if mysize:
                fig = pylab.figure(figsize=myfigsize(),**kwargs)
            else:
                fig = pylab.figure(**kwargs)
        else:
            fig = pylab.gcf()
            if not oplot:
                fig.clf()
    return fig

def mymoviesize():
    """
    return (12.80,9.60)

    for use as a 4:3 frame at 100 dpi
    """
    return (12.80,9.60)

def mymovie():
    """
    Generates new figure with mymovie size.
    """
    fig = pylab.figure(figsize=mymoviesize())
    fig.subplots_adjust(top=0.85)
    return fig

def new_plot(**kwargs):
    fig = myfig(**kwargs)
    ax = fig.add_subplot(111)
    return ax

def set_ultra_bold_rc():
    pylab.rcParamsDefault = pylab.rcParams
    pylab.rc("figure",figsize=[13,8])
    pylab.rc("axes", linewidth=3.5,labelsize=23,titlesize=26)
    pylab.rc("lines", markeredgewidth=2.0,linewidth=2.5)
    pylab.rc("legend",fontsize=20)
    pylab.rc("xtick",labelsize=26)
    pylab.rc("xtick.major",pad=20,size=20)
    pylab.rc("ytick",labelsize=26)
    pylab.rc("ytick.major",pad=20,size=20)
    pylab.rc("figure.subplot",bottom=0.3,left=0.2)
    pylab.rc("savefig",dpi=100)
    pylab.rc("axes",color_cycle=['b','g','r','c','m','y','k',(0.5,0,1),(1,0.5,0)])

def set_bold_rc():
    pylab.rcParamsDefault = pylab.rcParams
    pylab.rc("figure",figsize=[13,8])
    pylab.rc("axes", linewidth=2.5,labelsize=20,titlesize=23)
    pylab.rc("lines", markeredgewidth=1.0,linewidth=2.0)
    pylab.rc("legend",fontsize=18)
    pylab.rc("xtick",labelsize=23)
    pylab.rc("xtick.major",pad=14,size=18)
    pylab.rc("ytick",labelsize=23)
    pylab.rc("ytick.major",pad=14,size=18)
    pylab.rc("figure.subplot",bottom=0.15,left=0.15,
             top=0.9,right=0.9,hspace=0.25,wspace=0.25)
    pylab.rc("savefig",dpi=100)
    pylab.rc("axes",color_cycle=['b','g','r','c','m','y','k',(0.5,0,1),(1,0.5,0)])

def set_two_column_rc():
    pylab.rcParamsDefault = pylab.rcParams
    pylab.rc("figure",figsize=[3.375,2.31])
    pylab.rc("axes", linewidth=1.0,labelsize=10,titlesize=12)
    pylab.rc("lines", markeredgewidth=0.4,linewidth=1.2)
    pylab.rc("legend",fontsize=8,labelspacing=0.1,)
    pylab.rc("xtick",labelsize=10)
    pylab.rc("xtick.major",pad=6,size=8)
    pylab.rc("ytick",labelsize=10)
    pylab.rc("ytick.major",pad=6,size=8)
    pylab.rc("savefig",dpi=600)
    pylab.rc("axes",color_cycle=['b','g','r','c','m','y','k',(0.5,0,1),(1,0.5,0)])

def set_invert_rc():
    pylab.rcParamsDefault = pylab.rcParams
    rcParams = pylab.rcParams
    for key in pylab.rcParams.keys():
        if key.endswith('color'):
            if rcParams[key] == 'k':
                rcParams[key] = 'w'
                print ("Switched '"+key+"' from 'k' to 'w'")
            elif rcParams[key] == 'w':
                rcParams[key] = 'k'
                print ("Switched '"+key+"' from 'w' to 'k'")


def set_default_rc():
    pylab.rcdefaults()

def errorbartest(color,**kwargs):
    x = numpy.arange(50)
    y = numpy.sin(2*numpy.pi*x/25.)
    yerr = 0.5
    #fig = pylab.figure(1)
    #fig.clf()
    #ax = fig.add_subplot(111)
    ax = pylab.gca()
    ax.plot(x,y,**kwargs)

    myerrorbar(ax,x,y,yerr,color=color,**kwargs)

def myerrorbar_fuzzy(axis,x,y,yerr,color=(0,0,0),nlayers=20,label='',**kwargs):
    """Alternative to pylab.errorbar, with fading color intensity
    instead of discrete bars to denote width of uncertainty.  Uses
    multiple lines of same thickness offset by thickness to show
    fading. YERR gives the characteristic width of the distribution
    of data. Currently hard-coded to fade the input color to white
    using a Gaussian distribution of standard deviation YERR.
    """

    #make plots
    #if pylab.get_fignums() == []:
    #    fig = pylab.figure(figsize=(13,8))
    #else:
    #    fig = pylab.gcf()
    #    fig.clf()
    #axis = fig.add_subplot(111)


    #plot data so that plot data width is correct
    axis.plot(x,y,'w',**kwargs)
    y_data_width = axis.get_ylim()[1] - axis.get_ylim()[0]
    data_to_points = (axis.bbox.height/y_data_width)
    sigma_max = 3
    #determine linewidth
    width = 0.5*data_to_points*sigma_max*yerr/(2*nlayers + 1)  #total width of spread/number of lines
    #construct distribution
    colorindex = sigma_max*yerr*numpy.arange(nlayers)/(nlayers-1.0)  #want distribution to extend out to sigma_max
    distribution = numpy.exp(-(colorindex**2)/(2 * yerr**2)) #gaussian of height 1 (not normalized to area)

    #construct colors
    r0 = color[0]
    g0 = color[1]
    b0 = color[2]
    r = r0 + (1.0-r0)*distribution
    g = g0 + (1.0-g0)*distribution
    b = b0 + (1.0-b0)*distribution
    colors = []
    for layer in numpy.arange(nlayers):
        colors.append((r[layer],g[layer],b[layer]))

    colors.reverse()

    axis.plot(x,y,color=color,linewidth=width,label=label,**kwargs)
    for layer in numpy.arange(nlayers):
        shift = width*layer/data_to_points
        #plot layers of thick->thin lines, with increasing color intensity
        axis.plot(x,y+shift,color=colors[layer],linewidth=width,**kwargs)
        axis.plot(x,y-shift,color=colors[layer],linewidth=width,**kwargs)
        pylab.show()


#THICKNESS/COLOR DOESN'T GIVE GOOD IMAGE OF SPREAD
def myerrorbar_thickness(x,y,yerr,axis=None,color=(0,0,0),brightfactor=0.7,alpha=0.5,oplot=0,show_plot=0,**kwargs):
    """Alternative to pylab.errorbar, with fading color intensity
    instead of discrete bars to denote width of uncertainty.
    Uses successively thinner lines with darker colors overlaid
    to show fading.  YERR gives the characteristic width of the
    distribution of data.  Currently hard-coded to fade the input
    color to white using a Gaussian distribution of standard
    deviation YERR.
    """

    #build yerr array if necessary
    if numpy.asarray(yerr).size == 1:
        yerr = numpy.ones_like(y)*yerr
    else:
        assert yerr.size == y.size

    #make plots
    if axis == None:
        if pylab.get_fignums() == []:
            fig = pylab.figure(figsize=(13,8))
        else:
            fig = pylab.gcf()
            if not oplot:
                fig.clf()
        axis = fig.add_subplot(111)

    #plot data so that plot data width is correct
    axis.plot(x,y+yerr,alpha=0)
    axis.plot(x,y-yerr,alpha=0)
    y_data_width = axis.get_ylim()[1] - axis.get_ylim()[0]
    data_to_points = (axis.bbox.height/y_data_width)

    #determine linewidth
    sigma_max = 1
    envwidth = 0.5*data_to_points*sigma_max*yerr  #total width of spread

    #set color for envelope
    r0 = color[0]
    g0 = color[1]
    b0 = color[2]
    r = r0 + (1.0-r0)*brightfactor
    g = g0 + (1.0-g0)*brightfactor
    b = b0 + (1.0-b0)*brightfactor
    envcolor = (r,g,b)

    #plot thick envelope line
    axis.fill_between(x,y+yerr,y-yerr,facecolor=envcolor,edgecolor=envcolor,alpha=alpha)
    env = pylab.Rectangle((0,0,),1,1,fc=envcolor,alpha=alpha) #artist object representing envelope properties
    best = axis.plot(x,y,color=color,**kwargs) #line2d object representing best value properties

    if show_plot:
        pylab.show()

    return env,best










###########################ARRAY UTILITIES#############################
def get_log_series(xmin=1e-2,n_decades=3):
    """
    Returns array suitable for a logarithmically sampled axis,
    e.g. [0,1,2,5,10,20,50,100,200,500].
    """

    x = []

    for decade in range(n_decades):
        x.append(xmin*10**decade)
        x.append(2*xmin*10**decade)
        x.append(5*xmin*10**decade)

    x.append(xmin*10**n_decades)
    return numpy.asarray(x,float)


def in_limits(arr,limits,return_bool=False,
              return_values=False):
    """
    Returns indices where ARR lies within range
    set by LIMITS.  Optionally return ARR[indices]
    or T/F array of length ARR.size

    ARR:
        1D array.
    LIMITS:
        2-element list of form [arrmin,arrmax]
        corresponding to desired range of ARR.
    RETURN_VALUES = False:
        Boolean controlling whether indices or
        ARR[indices] is returned.  Indices only
        returned by default.
    RETURN_BOOL = False:
        Boolean controlling whether condition
        returned instead of indices.  Indices only
        returned by default.

    RETURNS:
       INDICES (tuple of arrays - must be tuple for arr[indices]
       to function properly with ND arrays) or ARR[INDICES]
       (1D array with length dependent on LIMITS).
    """
    condition = numpy.logical_and(arr >= limits[0],arr <= limits[1])
    if return_bool:
        return condition
    indices = numpy.where(condition)
    if return_values:
        return arr[indices]
    else:
        return indices


def absmin(arr,value=True):
    """
    Returns the of the ARR
    value closest to zero.
    return min(abs(ARR))
    """
    if value:
        return numpy.min(numpy.abs(arr))
    else:
        return numpy.argmin(numpy.abs(arr))


def absmax(arr,value=True):
    """
    Returns the of the ARR
    value from to zero.
    return max(abs(ARR))
    """
    if value:
        return numpy.max(numpy.abs(arr))
    else:
        return numpy.argmax(numpy.abs(arr))


def find_closest(arr,val,value=True):
    """
    Returns value or index of ARR
    closest to VAL.
    """
    if value:
        return arr[absmin(arr-val,value=False)]
    else:
        return absmin(arr-val,value=False)










###########################GAUSSIAN FITTING#############################
def gaussian1D(params):
    """
    func = gaussian((A,x0,width)) OR
    func = gaussian((A,x0,width,offset))

    OFFSET translates the zero point for the amplitude
    and scales the function so that the peak is still at A.
    """
    if len(params) == 4:
        A,x0,width,offset = params
        return lambda x: (A-offset)*numpy.exp(-((x-x0)/width)**2/2.) + offset
    elif  len(params) == 3:
        A,x0,width = params
        return lambda x: A*numpy.exp(-((x-x0)/width)**2/2.)
    else:
        print ("Parameter tuple to short or too long.")
        return None

def gauss_guess1D(x,y):
    """
    Uses moments of function to provide an intitial guess
    for fitting a Gaussian.
    """
    A = y.max()
    total = y.sum()
    mean = y.mean()
    x0 = (x*y).sum()/total
    width = numpy.sqrt(numpy.abs(
            numpy.sum((x-x0)**2 *y)/total))
    return A,x0,width

def fit_gaussian1D(x,y):
    """Returns gaussian function constructed with
    best fit parameters."""
    init_guess = gauss_guess1D(x,y)
    def errfunc(p,x,y):
        return y - gaussian1D(p)(x)

    fit,success = optimize.leastsq(errfunc,init_guess,args=(x,y))
    print ("Best fit parameters: ", fit)
    return gaussian1D(fit)

def gaussian2D(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*numpy.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def gauss_guess2D(data):
    """Returns (height, y0, x0, width_y, width_x)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    Y,X = numpy.indices(data.shape)
    x0 = (X*data).sum()/total
    y0 = (Y*data).sum()/total
    row = data[int(y0), :]
    col = data[:, int(x0)]
    width_x =  numpy.sqrt( numpy.abs(( numpy.arange(row.size)-x0)**2*row).sum()/row.sum())
    width_y =  numpy.sqrt( numpy.abs(( numpy.arange(col.size)-y0)**2*col).sum()/col.sum())
    height = data.max()
    return height, x0, y0, width_x, width_y

def fit_gaussian2D(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = gauss_guess2D(data)
    y,x = numpy.indices(data.shape)
    errorfunction = lambda p: numpy.ravel(gaussian2D(*p)(x,y) -
                                          data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def fit_gaussian(*args):
    print ("Deprecated in favor of fit_gaussian2D/1D.")
    #could use as chooser function instead
    pass









###########################CURL ROUTINES#############################
def get_curl(vx,vy,*dr):
    """
    Calculates the curl of a 2-D vector field.

    curl = get_curl(vx,vy,*dr)
    where *dr indicates a sample distance for
    each dimension of vx/vy.

    For 2D vx/vy, this is:
    curl = get_curl(vx,vy,[dx,dy])

    Can handle (ny,nx,nt,...) arrays, as long
    as the first two dimensions are y and x,
    respectively, and one sample distance is specified
    for each dimension of vx/vy.
    [This comes from python's row-major ordering]
    """

    dVxdy = numpy.gradient(vx,*dr)[0]
    dVydx = numpy.gradient(vy,*dr)[1]
    return dVydx - dVxdy

def get_unit_curl(vx,vy):
    """
    Calculates the curl of a 2D vector field,
    normalizing each vector to unit magnitude.
    Allows calculating the rotation tendency
    of a vector field with uneven magnitudes.

    curl = get_unit_curl(vx,vy)

    Can handle (ny,nx,nt,...) arrays, as long
    as the first two dimensions are y and x,
    respectively.
    [This comes from python's row-major ordering]
    """
    # normalize each vector to its magnitude
    mag = numpy.sqrt(vx**2 + vy**2)
    vx /= mag
    vy /= mag

    # reset to zero components that were divided
    # by zero magnitude
    if numpy.ma.isMaskedArray(vx):
        vx = numpy.ma.masked_invalid(vx)
        vy = numpy.ma.masked_invalid(vy)
    else:
        vx = numpy.ma.masked_invalid(vx).filled(fill_value=0.)
        vy = numpy.ma.masked_invalid(vy).filled(fill_value=0.)

    dVxdy = numpy.gradient(vx)[0]
    dVydx = numpy.gradient(vy)[1]
    return dVydx - dVxdy









###########################FLUCTUATION LEVEL ROUTINES#############################
def get_rms(t,signal_in,flims=[100.,1e6],sigma=[1],keep_mean=False):
    """Returns envelope amplitude of the input signal"""
    signal = signal_in.copy()
    # mask array if it has non-finite numbers
    if numpy.logical_not(numpy.isfinite(signal)).any():
        mask = True
        signal = numpy.ma.masked_invalid(signal)
    else:
        mask = False

    # convert integer to list
    try:
        temp = len(sigma)
    except TypeError:
        sigma = [sigma]

    if signal.ndim != len(sigma):
        if len(sigma) == 1:
            # assume the last dimension should be filtered
            kernel = numpy.zeros(signal.ndim)
            kernel[-1] = sigma[0]
    else:
        kernel = sigma
    # filter out low frequency components
    signal = ffilter(t,signal,flims[0],flims[1])
    # square and smooth
    signal = signal**2
    signal = gaussian_filter(signal,sigma=kernel)
    signal = numpy.sqrt(signal)

    if mask:
        return numpy.array(signal)
    else:
        return signal


def whiten_data(t,sig,noise_psd,psd_freq):
    """
    Divides out PSD of noise signal (passed by user)
    """
    from scipy.interpolate import interp1d

    # Calculate the freqs for the FT of the strain data:
    N = len(t)
    dt = t[1]-t[0]
    freqs = np.fft.rfftfreq(N, dt)

    # Calculate the interpolating function that we'll use to
    # resample the PSD onto the same frequency grid as the PSD:
    interp_psd = interp1d(psd_freq, noise_psd)

    # Whitening: transform to freq domain, divide by asd, then transform back,
    # taking care to get normalization right.

    sigf = np.fft.rfft(sig)
    wsigf = sigf / (np.sqrt(interp_psd(freqs) /dt/2.))
    wsig = np.fft.irfft(wsigf, n=N)
    return wsig

def compress_dynamic_range(t,signal,sigma,ratio=.8):
    """
    Runs an audio-style dynamic range compression on a
    signal, enhancing small amplitudes and reducing
    large ones.

    Uses get_rms to calculate an amplitude envelope.
    Constructs compressor functions for levels
    above/below envelope threshold (avg + one std).
    Applies compressor functions and returns
    compressed signal.
    """
    envelope = get_rms(t,signal,sigma=sigma)
    upper_threshold = envelope.mean() + envelope.std()/2.
    lower_threshold = envelope.mean() - envelope.std()/2.

    # normalize signal
    signal /= envelope.max()

    # convert amplitude to normalized dB
    #dBin = 20*log10(signal/signal.mean())
    # DOWNWARD COMPRESSION OF HIGH AMPLITUDES
    spread = envelope - upper_threshold
    newspread = spread*ratio
    newenvelope = upper_threshold + newspread
    gain = newenvelope

    # UPWARD COMPRESSION OF LOW AMPLITUDES


    #gain_range = pylab.linspace(-4*numpy.pi,4*numpy.pi,envelope.size)
    #gain = numpy.arctan(1. - gain_range)/numpy.pi + 0.5
    #gain[envelope > threshold] = ratio

    #??  need to pick this in a good way
    #gain = newenvelope/envelope

    # use arctangent transfer function
    # scale envelope to between -4pi and 4pi
    # want threshold to be at pi/2
    #envelope -= envelope.mean()
    #envelope = 0.5*numpy.pi*envelope/absmax(threshold)  # put threshold at pi/2
    #gain = numpy.arctan(envelope-threshold) # gives arctan output vs input curve

    return signal*gain


def root_drc(t,signal,sigma,power=5):
    """
    compressed_sig = mytools.root_drc(t,signal,sigma,power=5)

    Runs an audio-style dynamic range compression on a
    signal, enhancing small amplitudes by taking
    some fractional power of the signal amplitude
    to be the new amplitude.

    Returns zero-mean signal (fluctuations only).
    """
    envelope = get_rms(t,signal,sigma=sigma)
    upper_threshold = envelope.mean() + envelope.std()/2.
    lower_threshold = envelope.mean() - envelope.std()/2.

    # normalize signal
    sig_max = envelope.max()
    sig_mean = signal.mean()
    signal /= sig_max

    # subtract mean (redundant if done already)
    signal -= sig_mean

    return pow(numpy.abs(signal),1./float(power))*numpy.sign(signal)

def drc_max_gradient(images,):
    """
    Amplifies signal where image gradient is weak.
    (trying to avoid arbitrary velocimetry vectors
    where area is ~uniform in brightness)

    gain = 1./(1. + gmag), where gmag is the
    normalized magnitude of the image array gradient
    """
    grad = numpy.gradient(images)
    gmag = numpy.sqrt(grad[0]**2 + grad[1]**2)
    gmag /= absmax(gmag)
    gain = 1./(1. + gmag)

    return gain*images










###########################PLASMA-SPECIFIC ROUTINES#############################
def plasmaprops(n=1e19,Te_eV=3.2,Ti_eV=0.5,B_G=900.,mi_amu=40.,
                neutral_fraction=0.5):
    """Calculates common plasma parameters assuming quasineutrality
    and singly charged ions.

    Returns a dictionary of properties and values for the given parameters.
    N is the plasma density in units of m^{-3}
    Te and Ti are the electron and ion temperatures in eV
    B is the magnetic field strength in Gauss
    MI is the ion mass in units of the proton mass
    NEUTRAL_FRACTION is the ratio of neutral density to ion density

    Thermal speeds are calculated assuming particles have thermal kinetic
    energy = 3/2 kT.
    """

    props = {}
    q = 1.6e-19  #charge of electron
    me = 9.11e-31 #mass of electron
    twopi = 2*numpy.pi
    fourpi = 4*numpy.pi
    eps0 = 8.85e-12 #permittivity of free space
    u0 = fourpi*1e-7 #permeability of free space
    kB = 1.38e-23  #bolztmann constant

    #convert to SI units
    Te = q*Te_eV  # Joules
    Ti = q*Ti_eV
    B = B_G*1e-4  # Tesla
    mi = mi_amu*1.67e-27 # kg

    #BASIC
    props['T_e'] = Te_eV
    props['T_i'] = Ti_eV
    props['B'] = B_G
    props['n'] = n
    props['n_n'] = neutral_fraction*n
    props['mi'] = mi_amu

    #frequencies
    props['w_pe'] = numpy.sqrt(n * q**2/(me*eps0))
    props['w_pi'] = numpy.sqrt(n * q**2/(mi*eps0))
    props['f_pe'] = props['w_pe']/twopi
    props['f_pi'] = props['w_pi']/twopi
    props['w_ce'] = q*B/(me)
    props['w_ci'] = q*B/(mi)
    props['f_ce'] = props['w_ce']/twopi
    props['f_ci'] = props['w_ci']/twopi

    #lengths
    props['lambda_D'] = numpy.sqrt(eps0*Te/(n*q**2))

    #velocities
    props['v_the'] = numpy.sqrt(3*Te/me) #1/2 kT per degree of freedom
    props['v_thi'] = numpy.sqrt(3*Ti/mi)
    props['v_alfven'] = numpy.sqrt(B**2/(u0*n*mi))

    #dimensionless
    props['beta'] = 2*u0*n*(Te+Ti)/B**2
    props['N_D'] = (fourpi/3.)*n*props['lambda_D']**3
    props['logLambda'] = numpy.log(props['N_D']*3)


    #DERIVED

    #frequencies
    #collision frequencies estimated with v~v_th from Boyd & Sanderson, p 316
    sigma_prefactor = props['logLambda']*q**4/(twopi*eps0**2) #multiply by (1/mu*v**2)**2 to get sigma
    # then frequency is n*sigma*v
    props['nu_ee'] = n*props['v_the']*sigma_prefactor/(me*props['v_the']**2 )**2
    props['nu_ii'] = n*props['v_thi']*sigma_prefactor/(mi*props['v_thi']**2 )**2
    props['nu_ie'] = n*props['v_the']*sigma_prefactor*0.5/(mi*me*props['v_the']**4)
    props['nu_en'] = props['n_n']*1e-20*props['v_the']
    props['nu_in'] = props['n_n']*74e-20*props['v_thi']
    # Uses cross section for charge exchange:
    # pg. 712, J.B. Hasted, Charge Transfer and Collisional Detachment,
    # in Atomic and Molecular Processes, Ed. D.R. Bates, Pure and Applied Physics series, 1962

    # scaling only
    #props['nu_ee'] = n*(q**4)*props['logLambda']/(twopi*(eps0**2)*(me**2)*props['v_the']**3)
    #props['nu_ii'] = props['nu_ee']*numpy.sqrt(me*(Te**3)/(mi*(Ti**3)))
    #props['nu_ie'] = props['nu_ee']*(me*Te)/(mi*Ti)

    #velocities
    props['c_s'] = props['w_pi']*props['lambda_D']  #ion acoustic speed

    #lengths
    props['rho_le'] = numpy.sqrt(2./3.)*props['v_the']/props['w_ce'] #electron gyroradius
    props['rho_li'] = numpy.sqrt(2./3.)*props['v_thi']/props['w_ci'] #ion gyroradius
    props['rho_s'] = numpy.sqrt(2./3.)*props['c_s']/props['w_ci'] #ion sound gyroradius


    #dimensionless
    props['N_D'] = n*4.*numpy.pi*props['lambda_D']**3/3. #number of particles  in a Debye sphere.

    # should convert into class and then set attributes as well
    #for key,value in self.CineFileHeader.iteritems():
    #                            setattr(self,key,value)
    return props

def props2table(props,format_string="%3.1e"):
    """
    latexstring = props2table(props)

    Returns text to paste into a latex document
    containing the labels, numbers, and units
    for a plasmaprops dictionary.
    """
    # Assign custom (arbitrary) order to dictionary elements.

    key_list = ['B','T_e','T_i','beta','mi','n',
                'rho_le','rho_li','rho_s','lambda_D',
                'c_s','v_alfven','v_the','v_thi',
                'f_ce','f_ci','f_pe','f_pi',
                'w_ce','w_ci','w_pe','w_pi',
                'N_D','logLambda','n_n',
                'nu_ee', 'nu_en','nu_ie','nu_ii','nu_in']

    # List of latex names for items in key_list
    var_list = ['$B_0$','$T_e$','$T_i$','$\\beta$','$m_i$','$n$',
                '$\\rho_{e}$','$\\rho_i$','$\\rho_s$','$\\lambda_D$',
                '$c_s$','$v_A$','$v_{th,e}$','$v_{th,i}$',
                '$f_{ce}$','$f_{ci}$','$f_{pe}$','$f_{pi}$',
                '$w_{ce}$','$w_{ci}$','$w_{pe}$','$w_{pi}$',
                '$N_D$','ln$\Lambda$','$n_n$',
                '$\\nu_{ee}$', '$\\nu_{en}$','$\\nu_{ie}$','$\\nu_{ii}$','$\\nu_{in}$']

    # List of units for items in key_list
    units = ['G','eV','eV','','','m$^{-3}$',
             'm', 'm', 'm','m',
             'm/s','m/s','m/s','m/s',
             'Hz','Hz','Hz','Hz',
             'rad/s','rad/s','rad/s','rad/s',
             '','','m$^{-3}$',
             '1/s','1/s','1/s','1/s','1/s']

    # Make table
    table_str = "\hline \n Parameter & Value & Units \\ \hline \n"
    for ii,key in enumerate(key_list):
        #table_str += var_list[ii] + " & " + format_string%props[key] + " & " +units[ii] + " \\\\\n"
        table_str += var_list[ii] + " & " + format_string.format(props[key]) + " & " +units[ii] + " \\\\\n"

    table_str += "\hline"

    return table_str






#######################CORRELATION AND FFT ROUTINES#############################
def get_corr(t,sig1,sig2,detrend=pylab.detrend_linear,
             mode='same',normalized=1,optimize=1):
    """
    lag,corr = csdx.get_corr(t,sig1,sig2,detrend=pylab.detrend_linear,
                             mode='same',normalized=1,optimize=1):
    NORMALIZED uses (N-1)*sigma1*sigma2 to normalize
        correlation function (standard deviation normalization)
    OPTIMIZE calls OPTLENGTH on both signals so that the correlation
        runs quickly.  NOTE: Do NOT run without this on raw probe
        signals.  The number of points is absurdly un-optimized (odd,
        maybe close to prime, etc) and the traces are huge (1-2 Msamples).
    """
    #optimize lengths
    sig1 = optlength(sig1)
    sig2 = optlength(sig2)
    #detrend
    sig1 = detrend(sig1)
    sig2 = detrend(sig2)
    #correlate
    corr = numpy.correlate(sig1,sig2,mode=mode)
    #normalize
    if normalized:
        corr /= (t.size - 1)*sig1.std()*sig2.std()
    #calculate lag array
    dt = t[1]-t[0]
    tau = dt*(numpy.arange(corr.size) - corr.size/2)
    #integer division makes this agree with correlate
    #correlate leaves off the most positive lag value if the number of points is even

    return tau,corr

def corr_map(images,ref_pixel,verbose=False):
    """
    Makes map zero-lag correlation coefficient
    between a reference pixel and every other pixel.
    i.e. (nx,ny,nt) array becomes (nx,ny) array of
    correlation values between time series (nt)
    """
    ny,nx,nt = images.shape
    rx,ry = ref_pixel
    corr_map = numpy.zeros((ny,nx),float)
    for ii in range(nx):
        if verbose:
            if numpy.mod(ii,16)==0 : print ('Row: '+str(ii))
        for jj in range(ny):
            corr_map[jj,ii] = get_corr(numpy.arange(nt),
                                       images[ry,rx,:].squeeze(),
                                       images[jj,ii,:].squeeze(),
                                       mode='valid')[1]

    return corr_map

def corr_map_ext(images,ext_sig,zero_lag=True,
                 tlimits=[-100,100],verbose=False):
    """
    Makes map of zero-lag correlation coefficient
    between a set of images and an external time
    series of the same length.
    i.e. (nx,ny,nt) array becomes (nx,ny) array of
    correlation values between time series (nt)

    ASSUMES THAT SIGNALS ARE ON THE SAME TIME
    BASE / SAMPLED SIMULTANEOUSLY.
    (i.e. interpolate one to the times of the
    other before you call this routine.)

    If ZERO_LAG=False, then maps out maximum
    correlation within TLIMITS (in frames).
    """
    ny,nx,nt = images.shape
    corr_map = numpy.zeros((ny,nx),float)
    for ii in range(nx):
        if verbose:
            if numpy.mod(ii,16)==0 : print ('Row: '+str(ii))
        for jj in range(ny):
            if zero_lag:
                # get correlation coefficient for equal-length series
                corr_map[jj,ii] = get_corr(numpy.arange(nt),
                                           ext_sig,
                                           images[jj,ii,:].squeeze(),
                                           mode='valid')[1]
            else:
                # find absmax of correlation function within tlimits
                lag,corr = get_corr(numpy.arange(nt),
                                    ext_sig,
                                    images[jj,ii,:].squeeze(),
                                    mode='same')
                # limit to values with lag inside tlimits
                search_indices = in_limits(lag,tlimits)
                corr = corr[search_indices]
                # assign maximum +/- correlation value
                corr_map[jj,ii] = corr[absmax(corr,value=False)]

    return corr_map

def get_image_corr(t,images,images2=None):
    """
    corr_vs_lag = mytools.get_image_corr(t,images,images2=None)

    Returns correlation value (nt,nt) vs frame delay for an array
    of images (nx,ny,nt).

    Adding images2 allows cross correlation.

    Uses numpy.corrcoef to calculate one coefficient for
    each pair of frames.
    """
    try:
        nx = images.shape[1]
        ny = images.shape[0]
        nt = images.shape[2]
    except IndexError:
        print( "Image array must be of shape (nx,ny,nt).")
        return -1

    if isinstance(images2,numpy.ndarray):
        #double check shape
        if images2.shape != images.shape:
            print( "IMAGES2 must be the same shape as IMAGES.")
            return -1
        else:
            return numpy.corrcoef(images.flatten().reshape(nx*ny,nt),images2.flatten().reshape(nx*ny,nt),rowvar=0)
    else:
        return numpy.corrcoef(images.flatten().reshape(nx*ny,nt),rowvar=0)


def get_image_spectrum(t,images,row_avg=True,w=0,power=True,
                       plot=False,optimize=True):
    """
    f,spectrum = mytools.get_image_spectrum(t,images,**fftkwargs)

    Calculates a frequency power spectrum using the
    average frame correlation vs time given by get_image_corr.

    spectrum = FFT(corr)  #same shape as corr; FFT of each row
    """

    print ("Calculating frame correlation function...")
    corr = get_image_corr(t,images)
    # always a square matrix

    print( "Calculating FFT...")
    if optimize:
        nfft = optlength(corr[0,:].squeeze(),check=True)
    else:
        nfft = corr.shape[1]

    sigf = numpy.fft.fft(corr,n=nfft,axis=1)
    if power:
        sigf = numpy.abs(sigf)**2
    freq = numpy.fft.fftfreq(nfft,t[1]-t[0])
    freq = numpy.fft.fftshift(freq)
    sigf = numpy.fft.fftshift(sigf,axes=1)

    if w:
        print ("Smoothing with boxcar of width %d"%w)
        for ii,row in enumerate(sigf):
            sigf[ii,:] = smooth(row.squeeze(),window_len=w,
                                window='flat')

    if row_avg:
        print ("Averaging over lag combinations (rows).")
        sigf = sigf.mean(axis=0)

    return freq,sigf


# def optlength(signal,verbose=0,max_increment=20,check=False):
#     """
#     Uses J. Shipman's prime.py from NMT to check
#     for hard to factor numbers.
#
#     Returns more easily factored array length for faster FFTs.
#
#     if check:
#         return int(best_length)
#     else:
#         resizes 1D array to proper length
#     """
#     p = prime.Prime()  #Uses J. Shipman's prime.py from NMT
#
#     #calculate optimal number of time points (based on Tobin's optlength.pro)
#     length = signal.size
#     original_time = numpy.array(p.factorize(length)).sum()
#     best_length = length
#     best_time = original_time
#
#     for j in pylab.linspace(-max_increment,max_increment,num=2*max_increment + 1):
#         test_length = length+j
#         test_time = numpy.array(p.factorize(test_length)).sum() + numpy.abs(j)
#         #add j to weight towards smaller changes in size
#
#         if test_time <= best_time:
#             best_time = test_time
#             best_length = test_length
#         ipybreak()
#
#     if verbose:
#         if original_time != best_time:
#             print ('By changing '+str(length)+' to '+str(best_length)+',')
#             print (' FFT time will be reduced by a factor of '+str(original_time/best_time))
#         else:
#             print ('N = '+str(length)+' is already optimized')
#
#     # If we just want the new length, return it
#     if check:
#         return int(best_length)
#
#     #now adjust signal to be optimized length
#     if best_length == length:
#         return signal
#     elif best_length < length:
#         return signal[0:best_length]
#     else:
#         ndiff = best_length - length
#         return numpy.concatenate((signal,numpy.zeros(ndiff)))


def get_fft(t,signal,w=1,plot=False,log=True,optimize=True,window=True,**kwargs):
    """
    (freq, sigf) = get_fft(t,signal,log=0,w=1,plot=False,
    optimize=True,**kwargs)

    Parameters:
    ----------
    t:
        Time array (nt,)
    signal:
        Signal array (nt,)
    log:
        Plot on a log scale if True
    w:
        Width of boxcar smoothing filter
        (>1 for smoothing)
    plot:
        Plots result if True.
    optimize:
        Use optlength routine to
        truncate or pad arrays to speed
        up FFT calculation.
    **kwargs:
        Extra keywords passed directly
        to plot command.
    """
    import numpy
    if optimize:
        signal = optlength(signal)
    if window:
        hwindow = numpy.hanning(signal.size)
        signal *= hwindow
    sigf = numpy.fft.fft(signal)
    sigf = numpy.abs(sigf)**2
    indices = range(sigf.size)
    freq = numpy.fft.fftfreq(signal.size,t[1]-t[0])

    sigf = numpy.fft.fftshift(sigf)
    sigf = smooth(sigf,window_len=w)
    if window:
        sigf *= numpy.sqrt(8./3.)
    freq = numpy.fft.fftshift(freq)

    if plot:
        import pylab
        if log:
            pylab.semilogy(freq,sigf,**kwargs)
        else:
            pylab.plot(freq,sigf,**kwargs)

        pylab.xlabel('Frequency (Hz if t in seconds)')
        pylab.ylabel('Spectral Power')
        pylab.title('FFT Power Spectrum')

    return (freq,sigf)


def get_favg_spectrum(t,signal,fres=None,nd=None,demean=pylab.detrend_linear,
                      log=True,plot=False,optimize=True,
                      normalize=False,**kwargs):
    """
    Computes the frequency-averaged auto-spectral density as per the
    recipe in Bendat & Piersol, Section 11.5.5.

    f,Gxx = mytools.get_favg_spectrum(t,signal,fres=None,nd=None,
                                      demean=pylab.detrend_linear,
                                      log=True,plot=False,
                                      optimize=True,**kwargs)

    FRES is the desired frequency resolution of the averaged spectrum.

    ND is the number of frequency points averaged into each bin
        (equivalent to number of blocks for ensemble averaging)

    --> Must specify either FRES or ND.

    Mean is removed from signal using DEMEAN function.
        Default is pylab.detrend_linear().

    LOG determines whether to plot on a logarithmic scale.

    PLOT determines whether to plot the resultant spectrum.

    OPTIMIZE determines whether the length of the signal is optimized
        for fast FFT calculations.

    KWARGS are passed directly to the plot command.

    """
    if optimize:
        signal = optlength(signal)
    dt = t[1] - t[0]
    fsamp = 1./dt
    Tr = t.max() - t.min()
    df = 1./Tr
    Nr = signal.size

    if fres is None and nd is None:
        print ("Must specify either FRES or ND.")
        return 0,0

    if fres:
        nd = int(fres/df)
        print ("$n_d$ will be %d"%nd)
        fres = nd*df
        print( "Frequency resolution will be %2.1f Hz\n"%fres)
    else:
        fres = nd*df
        print ("Frequency resolution will be %2.1f Hz\n"%fres)

    signal = demean(signal)

    sigf = numpy.fft.fft(signal)
    # Gxx
    sigf = numpy.abs(sigf)**2
    # Gxx_hat if sampled at proper frequencies
    sigf = smooth(sigf,window_len=nd)

    # define center frequencies
    N = Nr/nd
    ii = (2*numpy.arange(0,N/2,dtype=float) + 1)*nd/2.
    ii0 = (1+nd)/2.
    freq = numpy.fft.fftfreq(Nr,dt)

    # extract center frequencies and corresponding average power
    #     Note that the 0th element excludes the zero frequency component,
    #     so the 0th frequency bin is averaged over (1,nd)
    favg = freq[ii.astype(int)]
    favg[0] = numpy.mean([freq[numpy.floor(ii0)],freq[numpy.ceil(ii0)]])
    Gxx = sigf[ii.astype(int)]
    Gxx[0] = sigf[1:nd+1].mean()

    if normalize:
        Gxx /= Gxx.sum()

    if plot:
        if log:
            pylab.semilogy(favg,Gxx,**kwargs)
        else:
            pylab.plot(favg,Gxx,**kwargs)

        pylab.xlabel('Frequency (Hz if t in seconds)')
        pylab.ylabel('Spectral Power')
        pylab.title('Auto-spectral Density, $n_d = %d$'%nd)

    return favg,Gxx


def lowpass_gfilter(t,signal,fc,return_mean=True):
    """
    Uses gaussian_filter to create a 'mean' signal,
    and returns either the mean or the fluctuations (remainder).

    filtsig = mytools.lowpass_gfilter(t,signal,fc,return_mean=True)
    """
    dt = t[1]-t[0]
    sigmat = 1./(2*numpy.pi*fc) # width of Gaussian smoothing kernel [s]
    sigmap = numpy.round(sigmat/dt)# width of gaussian smoothing kernel [points]
    moving_avg = gaussian_filter(signal,sigmap)  #sigma_f = 1kHz

    if return_mean:
        return moving_avg
    else:
        return signal - moving_avg


def ffilter(t,signal,fmin,fmax,**kwargs):
    """
    Choose subfunction depending on whether the signal is 1D
    or 3D (x,y,t).

    filtsig = ffilter(t,signal,fmin,fmax,notch=0,optimize=1,verbose=1,axis=-1,gsmooth=0)
    """

    # if masked, use raw data, then mask again
    if numpy.ma.isMaskedArray(signal):
        mask = signal.mask
        signal = signal.data
    else:
        mask = None

    if len(signal.shape) == 1:
        # Use old ffilter routine
        signal = _ffilter_1D(t,signal,fmin,fmax,**kwargs)
    else:
        signal = _ffilter_3D(t,signal,fmin,fmax,**kwargs)

    if mask is not None:
        signal = numpy.ma.masked_where(mask,signal)
    return signal

def _ffilter_1D(t,signal,fmin,fmax,notch=0,optimize=1,
                window=True,verbose=1,gsmooth=None):
    if optimize:
        signal = optlength(signal)
    nt = signal.size
    if window:
        hwindow = numpy.hanning(nt)
        signal *= hwindow
    sigf = numpy.fft.fftshift(numpy.fft.fft(signal))
    if window:
        sigf *= numpy.sqrt(8./3.)
    freq = numpy.fft.fftshift(numpy.fft.fftfreq(nt,t[1]-t[0]))

    if notch:  #currently not working
        if verbose:
            print ('Notch Filter Suppressing Frequency Band ['+str(fmin)+','+str(fmax)+'] Hz')
        plusfindex  = numpy.where(numpy.logical_and(freq>0,numpy.logical_and(freq <  fmin,freq >  fmax)))
        minusfindex = numpy.where(numpy.logical_and(freq<0,numpy.logical_and(freq > -fmin,freq < -fmax)))

    else: #bandpass is default
        if verbose:
            print ('Bandpass Filter Allowing Frequency Band ['+str(fmin)+','+str(fmax)+'] Hz')
        plusfindex  = numpy.where(numpy.logical_and(freq>0,numpy.logical_and(freq >  fmin,freq <  fmax)))
        minusfindex = numpy.where(numpy.logical_and(freq<0,numpy.logical_and(freq < -fmin,freq > -fmax)))

    #nfilter = plusfindex[0].size
    #hwindow = numpy.hanning(nfilter)
    #window = numpy.zeros(freq.size)
    #window[plusfindex] = hwindow
    #window[minusfindex]= hwindow

    window = numpy.zeros(freq.size)
    window[plusfindex] = 1.0
    window[minusfindex]= 1.0
    if gsmooth is not None:
        df = freq[1]-freq[0]
        nsmooth = numpy.floor(gsmooth/df/3.)
        window = gaussian_filter(window,nsmooth)
    signal = numpy.fft.ifft(numpy.fft.ifftshift(sigf*window))
    signal = numpy.real(signal) #abs(real())-abs() ~ 10^-24

    if nt > t.size:
        signal = signal[:t.size]
    else:
        signal = numpy.concatenate((signal,numpy.zeros(t.size-nt)))

    return signal

def _ffilter_3D(t,signal,fmin,fmax,notch=0,optimize=1,verbose=1):
    """Assumes signal.shape = (ny,nx,nt)."""
    nt = signal.shape[2]
    if optimize:
        best_length = optlength(numpy.arange(nt),check=True)

    signal = numpy.fft.fft(signal,n=best_length)
    freq = numpy.fft.fftfreq(best_length,t[1]-t[0])

    # Choose frequency ranges to suppress based on what to keep
    if notch:
        if verbose:
            print ('Notch Filter Suppressing Frequency Band ['+str(fmin)+','+str(fmax)+'] Hz')
        suppress = in_limits(numpy.abs(freq),[fmin,fmax],return_bool=True)
    else: #bandpass is default
        if verbose:
            print ('Bandpass Filter Allowing Frequency Band ['+str(fmin)+','+str(fmax)+'] Hz')
        suppress = numpy.logical_not(in_limits(numpy.abs(freq),[fmin,fmax],return_bool=True))

    # square window - could do better, but would have
    # to loop over.  extra window array takes up
    # too much memory
    signal[...,suppress] = 0.0
    signal[...,suppress] = 0.0
    signal = numpy.fft.ifft(signal)
    signal = numpy.real(signal) #abs(real())-abs() ~ 10^-24

    #if nt < t.size:
    #    t = t[:nt]
    if nt > t.size:
        signal = signal[:t.size]

    return signal


def _spectral_helper(x, y, NFFT=256, Fs=2, detrend=pylab.detrend_none,
        window=pylab.window_hanning, noverlap=0, pad_to=None, sides='default',
        scale_by_freq=None):
    #The checks for if y is x are so that we can use the same function to
    #implement the core of psd(), csd(), and spectrogram() without doing
    #extra calculations.  We return the unaveraged Pxy, freqs, and t.
    same_data = y is x

    #Make sure we're dealing with a numpy array. If y and x were the same
    #object to start with, keep them that way

    x = numpy.asarray(x)
    if not same_data:
        y = numpy.asarray(y)

    # zero pad x and y up to NFFT if they are shorter than NFFT
    if len(x)<NFFT:
        n = len(x)
        x = numpy.resize(x, (NFFT,))
        x[n:] = 0

    if not same_data and len(y)<NFFT:
        n = len(y)
        y = numpy.resize(y, (NFFT,))
        y[n:] = 0

    if pad_to is None:
        pad_to = NFFT

    if scale_by_freq is None:
        scale_by_freq = True

    # For real x, ignore the negative frequencies unless told otherwise
    if (sides == 'default' and numpy.iscomplexobj(x)) or sides == 'twosided':
        numFreqs = pad_to
        scaling_factor = 1.
    elif sides in ('default', 'onesided'):
        numFreqs = pad_to//2 + 1
        scaling_factor = 2.
    else:
        raise (ValueError("sides must be one of: 'default', 'onesided', or "
            "'twosided'"))

    # Matlab divides by the sampling frequency so that density function
    # has units of dB/Hz and can be integrated by the plotted frequency
    # values. Perform the same scaling here.
    if scale_by_freq:
        scaling_factor /= Fs

    if pylab.cbook.iterable(window):
        assert(len(window) == NFFT)
        windowVals = window
    else:
        windowVals = window(numpy.ones((NFFT,), x.dtype))

    step = NFFT - noverlap
    ind = numpy.arange(0, len(x) - NFFT + 1, step)
    n = len(ind)
    Pxy = numpy.zeros((numFreqs,n), numpy.complex_)

    # do the ffts of the slices
    for i in range(n):
        thisX = x[ind[i]:ind[i]+NFFT]
        thisX = windowVals * detrend(thisX)
        fx = numpy.fft.fft(thisX, n=pad_to)

        if same_data:
            fy = fx
        else:
            thisY = y[ind[i]:ind[i]+NFFT]
            thisY = windowVals * detrend(thisY)
            fy = numpy.fft.fft(thisY, n=pad_to)
        Pxy[:,i] = numpy.conjugate(fx[:numFreqs]) * fy[:numFreqs]

    # Scale the spectrum by the norm of the window to compensate for
    # windowing loss; see Bendat & Piersol Sec 11.5.2.
    Pxy *= 1 / (numpy.abs(windowVals)**2).sum()

    # Also include scaling factors for one-sided densities and dividing by the
    # sampling frequency, if desired. Scale everything, except the DC component
    # and the NFFT/2 component:
    Pxy[1:-1] *= scaling_factor

    #But do scale those components by Fs, if required
    if scale_by_freq:
        Pxy[[0,-1]] /= Fs

    t = 1./Fs * (ind + NFFT / 2.)
    freqs = float(Fs) / pad_to * numpy.arange(numFreqs)

    if (numpy.iscomplexobj(x) and sides == 'default') or sides == 'twosided':
        # center the frequency range at zero
        freqs = numpy.concatenate((freqs[numFreqs//2:] - Fs, freqs[:numFreqs//2]))
        Pxy = numpy.concatenate((Pxy[numFreqs//2:, :], Pxy[:numFreqs//2, :]), 0)

    return Pxy, freqs, t


def get_csd(x, y, NFFT=256, Fs=2, detrend_key='none',
        window=pylab.window_hanning, noverlap=0, pad_to=None, sides='default',
        scale_by_freq=None,return_nd=False,return_raw=False):
    """
    Returns (Pxy,f).
    Copy of pylab.csd without the attached plot.  Most of the
    function is the _spectral_helper() fucntion, except for the
    last if statement.  Pxy is an array of \~Sxy values for each
    block of NFFT points.

    Always takes blocks/FFTs along last axis of input arrays.
    """

    #The checks for if y is x are so that we can use the same function to
    #implement the core of psd(), csd(), and spectrogram() without doing
    #extra calculations.
    same_data = y is x

    #Make sure we're dealing with a numpy array. If y and x were the same
    #object to start with, keep them that way

    x = numpy.asarray(x)
    dims = list(x.shape)

    if not same_data:
        y = numpy.asarray(y)

    # zero pad x and y up to NFFT if they are shorter than NFFT
    if dims[-1]<NFFT:
        N = dims[-1]
        paddims = dims.copy()
        paddims[-1] = NFFT - N
        padding = numpy.zeros(paddims,dtype=float)
        x = numpy.concatenate((x,padding),axis=-1)

    if not same_data and len(y)<NFFT:
        N = dims[-1]
        paddims = dims.copy()
        paddims[-1] = NFFT - N
        padding = numpy.zeros(paddims,dtype=float)
        y = numpy.concatenate((y,padding),axis=-1)

    if pad_to is None:
        pad_to = NFFT

    if scale_by_freq is None:
        scale_by_freq = True

    # For real x, ignore the negative frequencies unless told otherwise
    if (sides == 'default' and numpy.iscomplexobj(x)) or sides == 'twosided':
        numFreqs = pad_to
        scaling_factor = 1.
    elif sides in ('default', 'onesided'):
        numFreqs = pad_to//2 + 1
        scaling_factor = 2.
    else:
        raise (ValueError("sides must be one of: 'default', 'onesided', or "
            "'twosided'"))

    # Matlab divides by the sampling frequency so that density function
    # has units of dB/Hz and can be integrated by the plotted frequency
    # values. Perform the same scaling here.
    if scale_by_freq:
        scaling_factor /= Fs

    # Generate container for output
    step = NFFT - noverlap
    ind = numpy.arange(0, x.shape[-1] - NFFT + 1, step)
    n_blocks = len(ind)
    Pxy_shape = list(x.shape) # Need to preserve shape except for frequency axis
    Pxy_shape[-1] = numFreqs  #    The number of frequencies is determined above
    block_shape = numpy.copy(Pxy_shape)
    Pxy_shape.append(n_blocks) # Need an index for each block
    Pxy = numpy.zeros(Pxy_shape, numpy.complex_)

    # Create window array
    if pylab.cbook.iterable(window):
        assert(window.shape == NFFT)
        windowVals = window
    else:
        windowVals = window(numpy.ones((NFFT,), x.dtype))
    # Calculate window power scaling factor
    window_sum = (numpy.abs(windowVals)**2).sum()

    # Tile to necessary size
    if len(block_shape) > 1:
        # Go through dimensions of block other than time,
        #    tiling the window array for each
        #Do this backwards so that taking the transpose at the end
        #    gives an array of the proper shape.
        for dim in block_shape[::-1][1:]:
            windowVals = numpy.tile(windowVals[...,numpy.newaxis],dim)
        windowVals = windowVals.transpose()

    # do the ffts of the slices
    for i in range(n_blocks):
        thisX = x[...,ind[i]:ind[i]+NFFT]
        thisX = windowVals * pylab.detrend(thisX,key=detrend_key,axis=-1)
        fx = numpy.fft.fft(thisX, n=pad_to, axis=-1)

        if same_data:
            fy = fx
        else:
            thisY = y[ind[i]:ind[i]+NFFT]
            thisY = windowVals * pylab.detrend(thisY,key=detrend_key,axis=-1)
            fy = numpy.fft.fft(thisY, n=pad_to, axis=-1)
        Pxy[...,i] = numpy.conjugate(fx[...,:numFreqs]) * fy[...,:numFreqs]

    # Scale the spectrum by the norm of the window to compensate for
    # windowing loss; see Bendat & Piersol Sec 11.5.2.
    Pxy *= 1. / window_sum

    # Also include scaling factors for one-sided densities and dividing by the
    # sampling frequency, if desired. Scale everything, except the DC component
    # and the NFFT/2 component:
    Pxy[1:-1] *= scaling_factor

    #But do scale those components by Fs, if required
    if scale_by_freq:
        Pxy[[0,-1]] /= Fs

    t = 1./Fs * (ind + NFFT / 2.)
    freqs = float(Fs) / pad_to * numpy.arange(numFreqs)

    if (numpy.iscomplexobj(x) and sides == 'default') or sides == 'twosided':
        # center the frequency range at zero
        freqs = numpy.concatenate((freqs[numFreqs//2:] - Fs, freqs[:numFreqs//2]))
        Pxy = numpy.concatenate((Pxy[...,numFreqs//2:, :], Pxy[...,:numFreqs//2, :]), -2)

    if return_raw:
        # return unaveraged (concatenated) block results
        if return_nd:
            return Pxy, freqs, n_blocks
        else:
            return Pxy, freqs

    #take this part from pylab.csd to average blocks
    if len(Pxy.shape) >= 2 and Pxy.shape[-1]>1: # Pxy is of shape (x.shape[:-1],NumFreqs,nd)
        Pxy = Pxy.mean(axis=-1)

    if return_nd:
        return Pxy, freqs, n_blocks
    else:
        return Pxy, freqs


def specgram(x, NFFT=256, Fs=2, detrend=pylab.detrend_none, window=pylab.window_hanning,
        noverlap=128, pad_to=None, sides='default', scale_by_freq=None):
    """
    Compute a spectrogram of data in *x*.  Data are split into *NFFT*
    length segements and the PSD of each section is computed.  The
    windowing function *window* is applied to each segment, and the
    amount of overlap of each segment is specified with *noverlap*.

    If *x* is real (i.e. non-complex) only the spectrum of the positive
    frequencie is returned.  If *x* is complex then the complete
    spectrum is returned.

    %(PSD)s

    Returns a tuple (*Pxx*, *freqs*, *t*):

         - *Pxx*: 2-D array, columns are the periodograms of
           successive segments

         - *freqs*: 1-D array of frequencies corresponding to the rows
           in Pxx

         - *t*: 1-D array of times corresponding to midpoints of
           segments.

    .. seealso::

        :func:`psd`
            :func:`psd` differs in the default overlap; in returning
            the mean of the segment periodograms; and in not returning
            times.
    """
    assert(NFFT > noverlap)

    Pxx, freqs, t = _spectral_helper(x, x, NFFT, Fs, detrend, window,
        noverlap, pad_to, sides, scale_by_freq)
    Pxx = Pxx.real #Needed since helper implements generically

    return Pxx, freqs, t

#---------------------------------------------------------------

#http://www.scipy.org/Cookbook/SignalSmooth

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """

    if x.ndim != 1:
        raise (ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise (ValueError, "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise (ValueError, "Window must be one of: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]



#----------------------------OTHER-----------------------------------

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''

    import numpy as n
    import scipy.interpolate
    import scipy.ndimage

    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print ("[congrid] dimensions error. " \
              "This routine currently only supports " \
              "rebinning to the same number of dimensions.")
        return None
    newdims = n.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )
            ipybreak()

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print ("Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported.")
        return None

#---------------------------------------------------------------
