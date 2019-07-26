#Ocean Optics System - SSX Summer 2018

Similar to the Vuv data, this code has a main file, and a helper file with a series of functions that perform the actual analysis.
The main file is OceanOptics.py. This function will import data saved from the Ocean Optics and perform additional analysis in python, with more user control than the ocean optics system. OceanOptics.py shoudl be edited with things like the filename, darkfilename, optional second file name, optional second darkfile name, if you want a dark subtract, if you want a plot that compares two different runs side by side, if you want an excel sheet printed with the peaks, ect.

ONlY OCEANOPTICS.PY should be edited, except it is a lot more likely for there to be bugs in this code. OceanView functions contains all the functions that actually do the lifitng.


Future modifications could include making this a class (OOP) instead of just scattered functions,
but it's a small project
