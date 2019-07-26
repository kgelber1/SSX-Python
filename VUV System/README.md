#VUV System - SSX Summer 2018

---------------------------------------------------------------------------------------------------------------------------------

To analyze the VUV data, you need to edit and run the Te.py file. 
Te.py does not take any parameters, it has to be edited to change things like the base, run date, shot numbers, density, 97 and 155 gain and other paramters. 
Then it has to be saved and run.
Te.py should be the only file that is edited from run to run


All the choices about how this code runs are done in Te.py. However, Te.py does not the analysis itself, it calls vuvFunctions.py which does the heavy lifting. For example of you want the average 97 photocurrent plotted, from Te.py you call vuv.plotPhotowithError, giving the function the appropriate paramters. Then vuvFunctions will take the data you sent and plot the data. VUVFUNCTIONS, CHECKTEMP and SMOOTH SHOULD NOT NEED TO BE CHANGED - unless there is a bug, or you want to change the way the code is analyzed. If you want the code to do or not do something, Te.py is the place to make those edits. 

vuvFunctions has two helper functions of it's own, CheckTemp which interpolates a curve from the PRISMSPECT data (which it's self was interpolated via RatioTempCurves.py - but the data has been saved so CheckTemp does not need to run RatioTempCurves). CheckTemp will convert the ratio to the temperature, handing it back to vuvFunctions. The other helper function is smooth.py, which I did not write, but found it very useful to smooth curves. vuvFunctions give smooth.py the data and smoothing parameters, smooth hands back smoothed data. 
