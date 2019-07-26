import vuvFunctions as vuv
import numpy as np


__author__ = "Katie Gelber"
__version__ = "1.0"


"""Main VUV file. ONLY EDIT THIS FILE!!
This is a place to set the specific parameters of the run, such as date, gain, shot number, ect,
and control what analysis you want performed/outputted.

Throughout all these functions: x =97, y=155"""

##########EDIT ALL THIS FROM RUN TO RUN
date = "070918"#Run Date
base = r'C:\Users\Katie\Documents\SSX'#base directory

#Set 97 parameters
xGain = 50#uA/V
xStart = 7
xEnd = 20
xRuns = np.arange(xStart, xEnd +1)
#Or to ommit noisy data
# xRuns = np.array([7,8,9,10,11,12,13,14])

#Set 155 Parameters
yGain = 50#uA/V
yStart = 21
yEnd = 28
yRuns = np.arange(yStart, yEnd +1)
# yRuns = np.arange(yStart, yEnd +1)


Density = '5e15'#Density of the plasma
ending =  '.txt.gz'
factor = 3.04
xGain = xGain * factor #Note that this will give a photocurrent that's
#already been corrected

vuv.setXWin(10, 100)#Sets the time scale for ALL plots
vuv.setYWin(0, 120)#Sets the height of Photocurrent plots
vuv.setChan('1', 2)#Scope, Channel
vuv.setRes(600, False)#DPI resolution, if you want file saved

#######DON'T CHANGE THESE:
print("Opening 97")
xDat = vuv.readFile(xRuns, ending, xGain, date, base)
xDat = vuv.findNoisyAverage(xDat)
xDat = vuv.addError(xDat)

print("Opening 155")
yDat = vuv.readFile(yRuns, ending, yGain,  date, base)
yDat = vuv.findNoisyAverage(yDat)
yDat = vuv.addError(yDat)
#OK you can edit again

xTitle = date + " - $\\bf{97.7nm}$, shots: " + str(xStart) + "-" +  str(xEnd)
yTitle = date + " -$\\bf{155nm}$, shots: " + str(yStart) + "-" +  str(yEnd)


#Plots every single run, useful to ommit noisy runs
# vuv.plotEachPhoto(xDat,xTitle)
# vuv.plotEachPhoto(yDat, yTitle)


base = base + "\\Pictures\\" + date
dest  = base + "_97.png"
vuv.plotPhotoWithError(xDat, xTitle, dest)#PLOTS AVERAGE
dest = base + "_155.png"
vuv.plotPhotoWithError(yDat, yTitle, dest)#PLOTS AVERAGE

#DO NOT CHANGE:
ratio  = vuv.findRatio(xDat, yDat, False)

# vuv.plotRatioWithError(ratio)

#DO NOT CHANGE:
temp = vuv.findTemp(ratio, Density, True)#True will print the max temp within the x-window that the graph is displayed in


dest  =base + "_Te.png"
title = date + " $\\bf{Te}$, 97.7: " + str(xStart) + "-" +  str(xEnd) +" , 155: "+ str(yStart) + "-" + str(yEnd)
vuv.plotTempwithError(temp, dest, title)

#If you want the temp data as a CSV:
# t = temp[0,:]
# y = temp[1,:]
# low = temp[3,:]
# err = y -low
# vuv.outCSV(t,y,err,  "070918 Te")
