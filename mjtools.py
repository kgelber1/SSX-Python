# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 16:40:10 2017

@author: Manjit
"""
from numpy import where, abs, min

import numpy as np
from scipy.ndimage import gaussian_filter1d



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


def get_corr(t,sig1,sig2,mode='same',normalized=1,optimize=1):
    """
    lag,corr = get_corr(t,sig1,sig2,mode='same',normalized=1,):
    NORMALIZED uses (N-1)*sigma1*sigma2 to normalize 
        correlation function (standard deviation normalization)
    OPTIMIZE calls OPTLENGTH on both signals so that the correlation
        runs quickly.  NOTE: Do NOT run without this on raw probe 
        signals.  The number of points is absurdly un-optimized (odd,
        maybe close to prime, etc) and the traces are huge (1-2 Msamples).
        
    """
    
    #n = len(sig1)
    #correlate
    corr = np.correlate(sig1,sig2,mode=mode)
    #normalize
    if normalized:
        corr /= (len(t) - 1)*sig1.std()*sig2.std()
    #    corr /= (1.0/(n-1))*sig1.std()*sig2.std()
    #calculate lag array
    dt = t[1]-t[0]
    tau = dt*(np.arange(corr.size) - corr.size/2) 
    #integer division makes this agree with correlate
    #correlate leaves off the most positive lag value if the number of points is even     
    return tau,corr


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
    tind = where(abs((timearr)-(timevalue))<delta)
    tind = tind[0][0]
    return tind
    
def tindex_min(timearr,timevalue):
    minval = min(abs((timearr)-timevalue))
    tind = where(abs((timearr)-(timevalue)) == minval)
    tind = tind[0][0]
    return tind
    
def tindex_near(timearr,timevalue,threshold):
    tinds = where(abs(timearr-timevalue) < threshold)
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
    c=83
    return m*x+c

    