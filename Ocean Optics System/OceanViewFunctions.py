import numpy as np
import csv
'''Ocean View Analysis Helper Functions
   Katie Gelber - SSX
   Summer 2018

   This code will subtract a dark frame from a run of the ocean optics visible spectrum. The Dark frame and the Data frame
   be of the same exposure time. This code will still work if the dark frame and data have different exposure times, but
   the result will not be as high quality. All that is needed is for the file locations to be specified (as .txt files, the
   default of the ocean optics)'''






def importFile(fileLoc):
    '''This code will open up the txt files containing the data.
     This fuction will remove an header data and return just the raw data.

    Parameters:
            fileLoc - the file location as a string
    Returns:
            a two dimentional array of strings:
                the first column is the wavelength
                the second column is the intensities '''
    text_file = open(fileLoc, "r")
    lines = text_file.read().split('\n')
    text_file.close()
    i = 0
    for j in range(0,len(lines)):
        if "Begin" in lines[j]:
            i=j+1
            break
    line = str(lines[i])
    text = line.split()
    append = np.array([text[0], text[1]])
    store = np.array(append)
    i+=1
    line = str(lines[i])
    text = line.split()
    append = np.array([text[0], text[1]])
    store = np.array([[store[0], store[1]], [text[0], text[1]]])
    i+=1
    while(i +1<len(lines)):
        line = str(lines[i])
        text = np.array([line.split()])
        store = np.append(store, text, axis  = 0)
        i+=1

    return(store)

def subtract(data, dark, bias):

    '''This function actually subtracts the dark frame from the raw data. Dark and Data must be the same length, corresponding
    to the same numbers of pixels. There is no error checking for this, I just assumed that since the ocean optics device
    is being used for both runs that pixel numbers should be the same.

    Parameters:
            dark - a two dimensional array of strings corresponding to the dark frame
            data - a two dimensional array of strings corresponding to the raw data
    Returns:
            x_dat - an array of floats corresponding to the wavelegths
            y_dat - an array of floats corresponding to the dark subtracted intesities
    '''
    x_dat = np.array([])
    y_dat = np.array([])
    sub = True
    if(type(dark)==int):
        sub = False

    run = int(len(data))


    for i in range(0,run):
        x = float(data[i,0])
        a = float(data[i,1])
        if(sub):
            b = float(dark[i,1])
            b = b*bias
            s = float(a-b)
        else:
            s = a
        x_dat = np.append(x_dat, x)
        y_dat = np.append(y_dat, s)

    return (x_dat, y_dat)

def getPeaks(x, y, tolerance):
    '''This function find the peak (spectral lines) from within the data. This is more accurate
    and versitile than the Ocean Optics one, and it can be printed in a table.

    Parameters:
            x - an array of floats corresponding to the wavelegths in nm
            y - an array of floats corresponding to the subtracted intesities
            tolerance - a number, either an integer or a float that can be used to adjust the tolerance.
                    A larger value for tolerance will result in fewer peaks being selected
                    A lower value for tolerance will result in more peaks being selected
    Returns:
            peakWav - a numpy array containing the wavelegths corresponding to the selected peaks
            peakVals - a numpy array containing the relative intensities of the selected peaks
            line - an array of floats, the same legth as the number of wavelegths sampled.
                    Every value in the array is equal to the cutoff threshold.

    '''
    maxIn = max(y)
    tol = maxIn/10 * tolerance

    #Inefficently creates an array that is the same length as y, but just contains the cutoff
    #threshold so it can be plotted with the spectrum if desired
    line = np.array([])
    for i in x:
        line = np.append(line, tol)
    if(tol<0):
        tol = tol*-1
    peakWav = np.array([])
    peakVals = np.array([])
    for i in range(0, len(x)):
        if(x[i]>750):
            break
        if((x[i]< 350)):
            continue
        if((y[i]< 0)):
            continue

        j = len(peakWav)-1
        if(y[i]>tol):
            if(j >= 0):
                curIntens = y[i]
                if(peakWav[j] +3 > x[i]):
                #then this is probably the same peak
                #so if the exisiting value is less, replace it
                    if(peakVals[j]< y[i]):
                        peakWav[j]= x[i]
                        peakVals[j]= y[i]
                else:
                    #just append it
                    peakWav = np.append(peakWav, x[i])
                    peakVals= np.append(peakVals, y[i])
            else:
                #just append it
                peakWav = np.append(peakWav, x[i])
                peakVals= np.append(peakVals, y[i])

    return(peakWav, peakVals, line)



def plotTwo(yLimit, x_dat, y_dat, line, cutoff,x,y):
    '''This will plot two sets of data on the same plot (in different colors), useful for comparing which
    run or dark frame is better.

    Parameters:
            yLimit - the maximum height of the graph(s). Set to zero for default handeling
            x_dat - an array of floats corresponding to the wavelegths, in nm
            y_dat - an array of floats corresponding to the subtracted intesities
            Line- an array of numbers, all equal to the intensity of what counts as a peak
            Cutoff- a Boolean of if you want the cutoff line to be shown
            x - a second array of floats corresponding to the wavelegths, in nm
            y - a second array of floats corresponding to the intesities
    Returns:
            None, outputs a graph with the two sets of data plotted on the same graph

    '''
    if(yLimit == 0 ):
        if(type(y)!=int):
            if(max(y_dat)<max(y)):
                yLimit = max(y)*1.1
            else:
                yLimit = max(y_dat)*1.1
        else:
            yLimit = max(y_dat)*1.1
    import matplotlib.pyplot as plt
    # %matplotlib inline
    plt.xlim([380,750])
    plt.ylim([0,yLimit])
    plt.plot(x_dat,y_dat,  label="Spectrum Data")
    if(type(x)!=int):
        plt.plot(x,y,  label="Second Spectrum Data")
    if(cutoff):
        plt.plot(x_dat, line,  label='Cutoff')
    plt.title("Visible Spectrum Data")
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')
    plt.show()

def subplotTwo(yLimit, x_dat, y_dat, line, cutoff,x,y):
    '''This will plot two sets of data on a subplot, useful for comparing which
    run or dark frame is better, while being able to see the full set of data (no overlap).

    Parameters:
            yLimit - the maximum height of the graph(s)
            x_dat - an array of floats corresponding to the wavelegths, in nm
            y_dat - an array of floats corresponding to the subtracted intesities
            Line- an array of numbers, all equal to the intensity of what counts as a peak
            Cutoff- a Boolean of if you want the cutoff line to be shown
            x - a second array of floats corresponding to the wavelegths, in nm
            y - a second array of floats corresponding to the intesities
    Returns:
            None, outputs a graph with the two sets of data subplotted

    '''
    if(yLimit == 0 ):
        if(type(y)!=int):
            if(max(y_dat)<max(y)):
                yLimit = max(y)*1.1
            else:
                yLimit = max(y_dat)*1.1
        else:
            yLimit = max(y_dat)*1.1

    import matplotlib.pyplot as plt
    f, (ax1, ax2) = plt.subplots(2,1, sharex = True, sharey = False)
    ax1.plot(x_dat,y_dat)
    ax1.set_title("Run 1")
    ax2.plot(x,y)
    ax2.set_title("Run 2")
    ax1.set_xlim([380,750])
    ax1.set_ylim([0,yLimit])
    ax2.set_xlim([380,750])
    ax2.set_ylim([0,yLimit])
    plt.show()


def plotLog(x_dat, y_dat, line, cutoff, dest):
    '''This will plot two sets of data on a subplot, useful for comparing which
    run or dark frame is better, while being able to see the full set of data (no overlap).

    Parameters:
            yLimit - the maximum height of the graph(s)
            x_dat - an array of floats corresponding to the wavelegths, in nm
            y_dat - an array of floats corresponding to the subtracted intesities
            Line- an array of numbers, all equal to the intensity of what counts as a peak
            Cutoff- a Boolean of if you want the cutoff line to be shown
            dest - a string representing the location the file should be saved to
    Returns:
            None, outputs a graph with the intestiy vs wavelength on a log scale
    '''
    if(yLimit == 0 ):
        yLimit = max(y_dat)*1.1
    import matplotlib.pyplot as plt
    plt.xlim([380,750])
    plt.ylim([30,20000])
    plt.semilogy(x_dat,y_dat,  label="Spectrum Data", color = 'darkslategrey', linewidth=2.0)
    if(cutoff):
        plt.plot(x_dat, line,  label='Cutoff')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.title("Visible Spectrum Data")

    plt.xlabel('$\\bf{Wavelength}$ ($\it{nm}$)', size = 14)
    plt.ylabel('$\\bf{Relative}$ $\\bf{Intensity}$', size = 14)
    plt.savefig(dest, format='png', dpi=600)
    plt.show()

def plotOne(yLimit, x_dat, y_dat, line, cutoff, dest):
    '''This will plot two sets of data on a subplot, useful for comparing which
    run or dark frame is better, while being able to see the full set of data (no overlap).

    Parameters:
            yLimit - the maximum height of the graph(s)
            x_dat - an array of floats corresponding to the wavelegths, in nm
            y_dat - an array of floats corresponding to the subtracted intesities
            Line- an array of numbers, all equal to the intensity of what counts as a peak
            Cutoff- a Boolean of if you want the cutoff line to be shown
            dest - a string representing the location the file should be saved to
    Returns:
            None, outputs a graph with the intestiy vs wavelength on a stardard scale
    '''
    if(yLimit == 0 ):
        yLimit = max(y_dat)*1.1
    import matplotlib.pyplot as plt
    # %matplotlib inline
    plt.xlim([380,750])
    plt.ylim([0,yLimit])
    plt.semilogy(x_dat,y_dat,  label="Spectrum Data", color = 'darkslategrey', linewidth=2.0)
    if(cutoff):
        plt.plot(x_dat, line,  label='Cutoff')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.title("Visible Spectrum Data")

    plt.xlabel('$\\bf{Wavelength}$ ($\it{nm}$)', size = 14)
    plt.ylabel('$\\bf{Relative}$ $\\bf{Intensity}$', size = 14)
    plt.savefig(dest, format='png', dpi=600)
    plt.show()

def output(x,y, title):
    '''A function that writes output to an excel file
    parameters:
            x - a numpy array of integers or floats corresponing to the x values
            y - a numpy array of integers or floats corresponding to the y values
            title - a string containing the tile of the file
    Outputs:
            an excel file with the x and y data'''
    # print("yeet")
    title = title + "1.xlsx"
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
