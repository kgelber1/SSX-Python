import numpy as np
import smooth as sm
import gzip
import regex as re
from scipy import interpolate
import CheckTemp as ct
import csv
import matplotlib.pyplot as plt
from matplotlib import rc

__author__ = "Katie Gelber"
__version__ = "1.0"

"""Core VUV SSX functions.

Contains functions that read in, analyze and plot the data for the
Vacumn Ultra-Violet Spectrometer. THIS WHOLE FILE SHOULD NOT NEED TO BE EDITED
unless there is a bug, or you want to change how the analysis is done."""


def setXWin(min, max):
    '''This function sets the x window for any and all graphs plotted
    Also maxes and mins will be calulated within this window
    Parameters:
    min - int/float of the minimum x. Logically should be >=0
    max - int/float of maximum x, must be greater than min'''
    global xmin#Yes, I know global variables are not a good idea, but they are
    #very convient in this case.
    xmin = min
    global xmax
    xmax = max


def setYWin(minn, maxx):
    '''This function sets the y for the intensity curves
    Parameters:
    min - int/float of the minimum y. Logically should be >=0
    max - int/float of maximum y, must be greater than min

    for default handeling set both to 0'''
    global ymin
    ymin = minn
    global ymax
    ymax = maxx


def setChan(scopeStr, channel):
    '''This function sets the scope and channel from which the VUV data
    was collected, to be used
    Parameters:
    scopeStr - the Scope number  (ie 1). Int or String.
    channel -  the channel number (ie 2). Int or String. '''
    global chan
    chan = int(channel)
    global scope
    scope = str(scopeStr)

def setRes(Res, Out):
    ''''A quick funtion to set the DPI resolution, and if you want
    the file to be saved'''
    global res
    res = Res
    global out
    out = Out

def readin(fileLoc):
    '''Given a file location with then ending .txt.gz, this function will unzip
    and read the data into an array containing the time and signal data

    Parameters:
    fileLoc - a string pointing to the file location
    Returns:
    times - numpy array containing the time information, in microseconds
    vals - numpy array containing the voltage information, in volts'''
    #unzips
    with gzip.open(fileLoc, 'rb') as f:
        lines = f.read()
    print('Opening Zipped File')
    line = str(lines)
    text = line.split("\\n")
    i = 1

    data = re.split(r'\\t', text[i])
    time = float(data[0])
    val = float(data[chan])
    times = np.array(time)
    vals = np.array(val)

    while i < len(text)-2:
        i +=1
        data = re.split(r'\\t', text[i])
        time = float(data[0])
        val = float(data[chan])
        times = np.append(times, time)
        vals = np.append(vals, val)
    return (times, vals)



def readtxt(fileLoc):
    '''Given a file location with then ending .txt, this function will read
    the txt data into an array containing the time and signal data. Specific to the format
    outputted by the VUV system. 
    Parameters:
    fileLoc - a string pointing to the file location
    Returns:
    times - numpy array containing the time information, in microseconds
    vals - numpy array containing the voltage information, in volts'''
    f= open(fileLoc)
    lines = f.read()
    print("Reading Text File")
    f.close()

    text = lines.split("\n")
    i = 1
    data = re.split(r'\t', text[i])
    time = float(data[0])
    val = float(data[chan])
    times = np.array(time)
    vals = np.array(val)

    while i < len(text)-2:
        i +=1
        data = re.split(r'\t', text[i])
        time = float(data[0])
        val = float(data[chan])
        times = np.append(times, time)
        vals = np.append(vals, val)
    return (times, vals)


def findPath(run, end, date, base):
    '''This function will create the path to the file. This function
    is specialized to my usernae and the location of where I put the files.
    This assumes the folder of information has been copied to
    the desired directory (via ssh into ion.swarthmore.physics.edu).
    Parameters:
    run - an int/string of the run number
    end - a string of the ending, usually .txt.gz
    Returns:
    fileLoc - a string pointing to the file location'''
    run = str(run)
    base = base + "\\" + date + "\\"
    fileLoc = date + 'r' + run + "-scope" + scope
    fileLoc = base + fileLoc +end
    return fileLoc

def findMax(times, vals, maxt, mint, error):
    '''This function will find the maxmimum value within some time range
    Parameters:
    times - numpy array containing the time information, in microseconds
    vals - numpy array containing the y-values
    mint - the minimum value of the time window
    maxt - the maximum value of the time window

    Returns:
    time - the time at which the maximim occurs
    max - the maximum value found within that range.'''
    max = 0
    err = 0
    for i in range(0, len(vals)):
        if(times[i]>maxt):
            if(times[i]<mint):
                if(vals[i]>max):
                    max = vals[i]
                    err = error[i]
    return (abs(max - err), max)

def findMin(times, vals, maxt, mint):
    '''This function will find the maxmimum value within some time range
    Parameters:
    times - numpy array containing the time information, in microseconds
    vals - numpy array containing the y-values
    mint - the minimum value of the time window
    maxt - the maximum value of the time window
    Returns:
    min - the minimum value found within that range.'''
    min = 1000
    time = 0
    for i in range(0, len(vals)):
        if(times[i]>maxt):
            if(times[i]<mint):
                if(vals[i]<min):
                    min = vals[i]
                    time = times[i]
    return (min)

def openFile(run, end, date, base):
    '''A function that finds the file path,
       calls the appropriate readin file depending on the ending of
       the file and returns the data
    Parameters:
        run - an int/string of the run number
        end - a string of the ending, usually .txt.gz
    Returns:
        times - numpy array containing the time information, in microseconds
        vals - numpy array containing the voltage information, in volts'''

    run = str(run)#uses the global variables
    channel = str(chan)
    fileLoc = findPath(run,end, date, base)
    if('.gz' in end):
        times, vals = readin(fileLoc)
    else:
        times, vals =readtxt(fileLoc)
    return(times, vals)

def findNoisyAverage(arr):
    '''A function that converts that computes the average at each point in
       time.
    Parameters:
       array - an array containing the time and photocurrents from each run
    Returns:
        arr- numpy array containing the time information, followed by
       the photocurrent from several runs with the average appened on.
    '''
    ave = np.array([])
    for i in range(0, arr.shape[1]):
        sum = 0
        count = 1
        for j in range(1, arr.shape[0]):
        #IDK IF DATA IS NOISY:
            if(arr[j][i]>0):
                sum += arr[j][i]
                count+=1
        if(count > 1):
            count-=1
        ave  = np.append(ave, sum/count)
    arr = np.vstack((arr, ave))
    return arr


def readFile(runs, ending, gain, date, base):
    '''A wrapper funtion that reads in the file, applies smoothing
       converts voltage to photocurrent and appends it to an array
       for each run number given in the array of runs
    Parameters:
       run - an array of ints/strings of the run numbers
       ending - a string of the ending, usually .txt.gz
       gain - the gain set on the box, in microAmps/Volt
    Returns:
        arr- numpy array containing the time information, followed by
       the photocurrent from several runs with the average appened on.'''
    #Set up this way to ensure that the formating is desired.
    t, y =openFile(runs[0], ending, date, base)
    y  = sm.main(y, 23, 25)
    y = y*gain
    # diff =
    # t = t[0:len(y)]
    diff = (len(t) - len(y))/2
    diff = int(diff)

    t = t[diff:len(t)-diff]
    array = np.array([t, y])
    i = 1
    for i in range(1, len(runs)):
        t, vals =openFile(runs[i], ending, date, base)
        vals  = sm.main(vals, 23, 25)
        vals =vals*gain
        array = np.vstack((array, vals))
    return array

def addError(arr):
    '''A funtion that comuputs the standard error - the standard
       deviation and the standard deviation from the mean
    Parameters:
      arr- numpy array containing the time information, followed by
      the photocurrent from several runs with the average appened on.
    Returns:
       arr- returns the intial array with the standard deviation,
            the standard deviation from the mean, the upper error
            and the lower error for the set of runs'''

    std_dev = np.array([])#sigma
    std_dev_mean = np.array([])#sigma/srt(N)
    above = np.array([])
    below = np.array([])
    n = arr.shape[0]-1#N
    for i in range(0, arr.shape[1]):
        sum = 0
        for j in range(1, n-1):
            sum += (arr[j][i] - arr[n][i])**2
        std = (sum/(n-1))
        std = std**(1/2)
        std_mean = std/(n**(1/2))
        std_dev = np.append(std_dev, std)
        std_dev_mean = np.append(std_dev_mean, std_mean)
        above = np.append(above, arr[n][i] + std)
        below = np.append(below, arr[n][i] - std)

    arr = np.vstack((arr, std_dev))
    arr = np.vstack((arr, std_dev_mean))
    arr = np.vstack((arr, above))
    arr = np.vstack((arr, below))

    return arr

def plotEachPhoto(arr, title):
    '''A funtion that plots the photocurrent of every single run
    Parameters:
      arr- numpy array containing the time information, followed by
         the photocurrent from several runs
      title - a string contianing the title of the graph, ie 97 or 155
    Outputs:
       no returns, but it outputs a plot of the photocurrents
       The xmax and xmin are set by the global varibales'''
    t = arr[0,:]
    end = arr.shape[0] -1
    y = arr[end-4,:]
    up = arr[end-1,:]
    low = arr[end,:]
    err = y-low-30
    for i in range(1,end-4):
        y = arr[i,:]
        run = "Run" + str(i) + "VUV Data"
        plt.plot(t, y,  label=run)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14)

    plt.title(title)
    plt.xlabel('$\\bf{Time}$ ($\it{\mu s}$)', size = 14)
    plt.ylabel('$\\bf{Photo Current}$ ($\it{\mu A}$)', size = 14)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.legend()
    plt.show()
    # if(out):
    #     plt.savefig(dest, format='png', dpi=600)

def plotPhotoWithError(arr, title, dest):
    '''A funtion that plots the photocurrent
    Parameters:
      arr- numpy array containing the time information, followed by
         the photocurrent from several runs with the average appened
         and the error information
      title -  a string contianing the title of the graph, ie 97 or 155
    Outputs:
       no returns, but it outputs a plot of the photocurrent
       with error bars. The xmax and xmin are set by the global varibales'''
    t = arr[0,:]
    end = arr.shape[0] -1
    y = arr[end-4,:]
    up = arr[end-1,:]
    low = arr[end,:]
    up = up
    low = low
    err = y-low
    #if it's too choppy. I've smoothed the data
    #and added a +5 error to cover all the bases
    # up = sm.iter_smooth(up, 40,45) + 5
    # low = sm.iter_smooth(low, 40,45)-5
    plt.plot(t, y,  label="Averaged VUV Data", linewidth=2.0)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14)

    #make sure the array is centered
    diff = (len(t) - len(up))/2
    diff = int(diff)
    t = t[diff:len(t)-diff]

    #I could have used errorbars, but I did it by hand, plotting the upper, lower error
    #and filling between.
    plt.fill_between(t,up,low, interpolate=True, alpha = 0.1)
    plt.plot(t, up,  label="Upper Error", linestyle = "dashed")
    plt.plot(t, low,  label="Lower Error", linestyle = "dashed")
    plt.title(title)
    plt.xlabel('$\\bf{Time}$ ($\it{\mu s}$)', size = 14)
    plt.ylabel('$\\bf{Photo Current}$ ($\it{\mu A}$)', size = 14)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    # plt.legend()
    if(out):
        plt.savefig(dest, format='png', dpi=600)
    plt.show()

def mean(arr):
    '''A helper funtion that computes the average of a given array'''
    sum = 0
    for i in range(0, len(arr)):
        sum+=arr[i]
    ave = sum/len(arr)

    return ave

def cleanRatio(rat, t):
    '''Not currently in use!
    This code finds any rapid spikes in the ratio and smooths them over an adjustable window
    Param -
            rat - the ratio of 97.7/155 as an array of real numbers
            t - an array of times (-20 to 160 us)
    '''
    # dt = len(rat)/200
    #Takes the derivative
    dydx = np.diff(rat)/np.diff(t)
    if len(dydx) < len(t):
        dydx = np.append(dydx, dydx[1])#matches the length

    for j in range(0, 2):#Actuallly cleans twice over
        winSize = 500#adjustable window of ratio smoothing
        for i in range(0, len(rat)):
            if (dydx[i] > 0.2 or dydx[i] < -0.2):#if sudden change in slope
                if(i < winSize):#ensures that there won't be any out of bounds errors
                    start = 0
                else:
                    start = i -winSize
                if(i+winSize >= len(rat)) :
                    end = len(rat) -1
                else:
                    end = i+winSize
                # snippit = rat[start:end]
                if(i+10 >= len(rat)) :
                    Tend = len(rat) -1
                else:
                    Tend = i+100
                if(i > 100) :
                    sub = mean(rat[i-100:Tend])
                else:
                    sub = 0
                rat[i] = mean(rat[start:end]) - sub #The window over which the ratio is averaged
                #excludes the immediate area around the error point. For reference, there are
                #just over 20,000 points in the array
        plt.plot(t, dydx, label = str(j))
        dydx = np.diff(rat)/np.diff(t)
        if len(dydx) < len(t):
            dydx = np.append(dydx, dydx[1])
    #Plots the cleaned derivative
    plt.plot(t, dydx, label = "Final")
    plt.ylim(-1000, 1000)
    plt.legend()
    plt.show()
    return rat #returns the cleaned ratio and time

def findRatio(xArr, yArr, min):
    '''One of the bigger funtions that finds the ratio of the
       97 (denoted x) and the 155 (denoted y)
       Parameters:
         xArr- an array contianing all the information for the 97 photocurrent
         yArr- an array contianing all the information for the 97 photocurrent
         min- a boolean of if you want the minimum ratio to be found
       returns:
          arr - an array contiang the time,ratio, error, upper error limit'
              lower error limit, and the percent error'''
    t = yArr[0, :]
    xAve = xArr[xArr.shape[0]-5, :]
    yAve = yArr[yArr.shape[0]-5, 0:len(xAve)]

    #This discounts any zero or negative points (or points close enough to zero that
    #the code would take a very long time to analyze)
    #Setting to 0.1 vs 0.001 or 0 had no impact on the final result of Te, but slowed
    #down runtime considereably, hence 0.1
    xAve  = np.ma.masked_less_equal(xAve, 0.1)
    yAve  = np.ma.masked_less_equal(yAve, 0.1)
    ratio = xAve/yAve#97.7/155

    if(min == True):
        minR = findMin(t,ratio, xmin, xmax)
        print("Min Ratio of ", minR, "so max temp... ", ct.vuvLookupCurve(minR, 3))

    #This is all useless
    xErr =  xArr[xArr.shape[0]-3, :]
    xErr = abs(xErr-xAve)/xAve
    yErr = yArr[yArr.shape[0]-3, 0:len(xAve)]
    yErr =  abs(yErr-yAve)/yAve
    err = (xErr**2 + yErr**2)**(1/2)
    up = ratio + ratio *err
    low =ratio - ratio *err
    Pr  = err #percent error in ratio

    arr = np.vstack((t, ratio))
    arr = np.vstack((arr, err))
    arr = np.vstack((arr, up))
    arr = np.vstack((arr, low))
    arr = np.vstack((arr, Pr))#yes, this is reduntant. It was there
    #in case I needed to do something

    return arr

def plotRatioWithError(arr):
    '''A funtion that  plots the ratio with error
    Parameters:
      arr- numpy array the ratio information
    Outputs:
       no returns, but it outputs a plot of the ratio
       with error bars. The xmax and xmin are set by the global varibales'''
    t = arr[0,:]
    ratio = arr[1,:]
    up = arr[3,:]
    low = arr[4,:]

    plt.plot(t, ratio,  label="Ratio of 97.7/155")
    # plt.plot(t, up,  label="Upper Error", linestyle = "dashed")
    # plt.plot(t, low,  label="Lower Error", linestyle = "dashed")
    plt.fill_between(t,up,low, interpolate=True, alpha = 0.1)
    plt.title("Ratio of 97.7/155")
    plt.xlabel('$\\bf{Time}$ ($\it{\mu s}$)', size = 14)
    plt.ylabel('Ratio')

    plt.xlim(xmin, xmax)
    plt.ylim(0, 15)
    # plt.legend()
    plt.show()

def plotRatio(arr):
    '''A funtion that  plots the ratio without error. Probably
       not super useful
    Parameters:
      arr- numpy array the ratio information
    Outputs:
       no returns, but it outputs a plot of the ratio
       with error bars. The xmax and xmin are set by the global varibales'''
    t = arr[0,:]
    end = arr.shape[0] -1
    y = arr[end-3,:]

    plt.plot(t, y,  label="Ratio of 97/155", linewidth=2.0)
    plt.title("Ratio of 97/155")
    plt.xlabel('Time (us)')
    plt.ylabel('Ratio')

    plt.xlim(xmin, xmax)
    # plt.legend()
    plt.show()

def correct(up, low, t):
    '''Makes sure that the upper and lower limits are what they should be'''
    for i in range(0, len(up)):
        if(up[i]>t[i]):
            temp = up[i]
            up[i] = low[i]
            low[i] = temp
    return up,low

def plotTempwithError(arr, dest, title):
    '''A funtion that  plots the temperature with error propogated
       from the ratio
    Parameters:
      arr- numpy array the ratio information
    Outputs:
       no returns, but it outputs a plot of the ratio
       with error bars. The xmax and xmin are set by the global varibales'''
    t = arr[0,:]
    y = arr[1,:]
    end = arr.shape[0]-1

    up =arr[2, :]

    low =arr[3,:]
    up =2*y-low #becuase up hits zero, log zero breaks code, assign up with refernce to low
    up, low = correct(up, low, y)
    plt.plot(t, y,  label="Electron Temperature", linewidth=2.0)
    plt.plot(t, up,  label="Upper Error", linestyle = "dashed", alpha = 0.5)
    plt.plot(t, low,  label="Lower Error", linestyle = "dashed", alpha = 0.5)

    plt.fill_between(t,up,low, interpolate=True, alpha = 0.1)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.title(title)
    plt.xlabel('$\\bf{Time}$ ($\it{\mu s}$)', size = 14)
    plt.ylabel('$\\bf{Temperature}$ ($\it{e}$V)', size = 14)
    plt.xlim(xmin, xmax)

    plt.ylim(0, 12)
    if(out):
        plt.savefig(dest, format='png', dpi=600)
    plt.show()

def findTemp(arr, den, max):
    '''One of the bigger funtions that finds the temperature of the
       plasma based on the ratio
       Parameters:
         arr- an array contianing all the ratio information
         den- the density of the plasma (1e14, 5e14, 2e15, 5e15)
         max- a boolean of if you want the maximum ratio to be found
       returns:
          arr - an array contiang the time, teperature, upper error limit'
              lower error limit, and the percent error'''
    if('2e15' in den or '2 e15' in den or '2 e 15' in den):
        density = 2
    elif('5e14' in den or '5 e14' in den or '5 e 14' in den):
        density = 1
    elif ('1e14' in den or '1 e14' in den or '1 e 14' in den):
        density = 0
    elif ('5e15' in den or '5 e15' in den or '5 e 15' in den):
        density = 3
    else:
        print("INCORRECT DENSITY, only 5e15, 2e15, 5e14 or 1e14 are options")
        print("defaulting to density = 2e15")
        density = 2
    n = arr.shape[0]-1
    t = arr[0,:]
    y = arr[1,:]
    up = arr[2,:]
    low = arr[3,:]
    err = arr[4,:]

    temp  = ct.vuvLookupCurve(y, density)
    upperLimit = ct.vuvLookupCurve(up, density)
    lowerLimit = ct.vuvLookupCurve(low, density)

    tempArr = np.vstack((t,temp))
    tempArr = np.vstack((tempArr,upperLimit))
    tempArr = np.vstack((tempArr,lowerLimit))
    tempArr = np.vstack((tempArr,err))

    if(max == True):
        maxError, maxT = findMax(t,temp, xmin, xmax, lowerLimit)
        print("Maximum Temperature of ", maxT)
        print("with error of +- ", maxError)
    return tempArr

def outXLSX(x, y, title):
    '''Outputs an x and y array as a xlsx file'''
    title = title + ".xlsx"
    import pandas as pd
    import numpy as np
    xT = x.copy()
    yT = y.copy()
    array = np.vstack((xT, yT))
    array = np.transpose(array)
    out = pd.DataFrame(array)
    writer = pd.ExcelWriter(title, engine='xlsxwriter')
    out.to_excel(writer)
    writer.save()

def outCSV(x, y, err, title):
    '''This one outputs an x, y, and error column as a tab seperated
    value txt file. Changing the ending and the delimiter to a comma would outputs
    as csv'''
    xT = x.copy()
    yT = y.copy()
    eT = err.copy()
    array = np.vstack((xT, yT, eT))
    array = np.transpose(array)
    title = title + ".txt"
    np.savetxt(title, array, delimiter=" ")
