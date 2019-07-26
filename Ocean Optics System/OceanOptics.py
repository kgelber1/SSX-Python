import OceanViewFunctions as oo
from astropy.table import Table, Column

'''Ocean View Analysis
   Katie Gelber - SSX
   Summer 2018

   This code will subtract a dark frame from a run of the ocean optics visible spectrum. The Dark frame and the Data frame
   be of the same exposure time. This code will still work if the dark frame and data have different exposure times, but
   the result will not be as high quality. All that is needed is for the file locations to be specified (as .txt files, the
   default of the ocean optics)'''

#Specify the File location of the Data file:
datfile = r'D:\070518\070518_10s_HR4C31931_0012_dat.txt'
#Specify the File location of the dark file:
darkfile = r'D:\070518\070518_10s_HR4C31931_0009.txt'

#Don't change:
file = oo.importFile(datfile)
darkFile = oo.importFile(darkfile)
x, y = oo.subtract(file, darkFile, 1)#Last param is a bias for the dark frame

tolerance = .10#control what intensity counts as a peak. Lower tolerance = more peaks

#Don't change:
xpeak, ypeak, line = oo.getPeaks(x, y, tolerance)


# oo.output(xpeak,ypeak, "Peak Wavelengths") #Outputs to Excel
ylim =400 #Set to zero for defualt handeling. OTherwise set the maximum height
#of the graph

#Specify where the final file should be saved:
dest  = r'C:\Users\Katie\Documents\SSX\Pictures\070518_visible.png'
oo.plotOne(ylim, x, y,line, True, dest)#plots one photocurrent. The boolean controls
#if you want the line that determines what counts as a peak to run.
oo.plotLog(x, y,line, True, dest)#Same thing, but on a log axis


#Outputs table of peaks to Command Line
t = Table([xpeak, ypeak], names=('Wavelegth', 'Intensity'))
print(t)
print("total of ", len(t), "peaks counted")
##Outputs same peaks to an excel file for easy copy and pasting
# #
# #SECOND FRAME:###########
# file2 = oo.importFile(downfile)
# darkFile2 = oo.importFile(downdarkfile)
# xx, yy = oo.subtract(file2, darkFile2, 1)
# ylim =300 #Set to zero for defualt handeling
# # oo.plotTwo(ylim, x, y,line, True, xx, yy)
# oo.subplotTwo(ylim, x, y,line, True, xx, yy)#Useful for compaing dark frames!
