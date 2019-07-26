import numpy as np
import BField_Plotter
# reload(BField_Plotter)
# import MagBField_Plotter_Luke
# reload(MagBField_Plotter_Luke)
# import process_mjmag_data as mj
import process_mag_data_4x4 as pmd
import iter_smooth as ism
import matplotlib.pylab as plt
import mytools as my
import ssx_data_read_basic as sdr

##############################################
#shot = '020618r2'


day = '061219'#''051917'#
shot = day+'r13'#'40'#
t0 = 35
tf = 150

# time,bdot,timeb,b,bmod, data = mj.process_mjmag_data(shot)
time,timeb,b,bmod, data = pmd.get_Bxy_vecRendering_format(shot)
timeb=timeb-timeb[0]-2

######### fixing bounds using a function called fix_bounds ###########
t_index, t = my.fix_bounds(np.array((t0, tf)), timeb)
t0 = t_index[0]
tf = t_index[-1]

sample_Freq = 10 # sampling frequency
N = 16
b = b[:,:,t0:tf]
s = np.arange(N) * 1.5# - 2.86
x = b[0,0:N,:].T[0:-1:sample_Freq]
y = b[1,0:N,:].T[0:-1:sample_Freq]
t = timeb[t0:tf][0:-1:sample_Freq]

path = 'Data\\2018\\'+day+'\\'+day+'-Analysed\\'


############### Plotting the Bx and By vector rendering  ##############################
B = BField_Plotter.BField_xy_Animator(shot,t,s,x,y) #Instantiate the BField_xy_Animator Object
# B.gen_Arrows(100) #Plot the vectors at some time, in this case t[10]
animation = B.make_animation() #Now create an animation
B.save_animation(path,animation) #Save the animation
plt.show()
