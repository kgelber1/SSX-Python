import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import interpolate
from scipy.interpolate import interp1d

'''Use this to make a graph of the McPherson MagF grating efficiency'''

#data read in by hand
x = np.array([90,93,99, 110, 115, 125, 135,150,155, 175,215,255,310,345])
y = np.array([0.02,0.05,0.12, 0.23, 0.28, 0.320, 0.330, 0.325,0.320, 0.3,0.25,0.2,0.15,0.11])

#interpolate curve fit
f= sp.interpolate.interp1d(x,y, kind = 'cubic')
ynew = f(x)
print("97.7: ", f(97.7))
print("155: ", f(155))
print("so a good ratio is: ", f(155)/f(97.7) ) #prints the ideal ratio to use

#outputs plot
plt.plot(x,ynew)
plt.scatter(x,y)
plt.xlim(50,375)
plt.ylim(0, 0.4)


plt.xlabel('$\\bf{Wavelength}$ ($\it{nm}$)', size = 14)
plt.ylabel('$\\bf{Theoretical}$ $\\bf{Efficiency}$ ($\%$)', size = 14)

plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=14)
plt.title("McPherson Grating Efficiency")
plt.savefig('grating.png', format='png', dpi=600)
# plt.grid()
plt.show()
