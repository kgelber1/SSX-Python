from __future__ import division, print_function, absolute_import

import numpy as np
import matplotlib.pylab as plt
import os


import MagBField_Plotter_Luke
import ssx_data_read_basic as sdr
import process_mag_data_4x4 as pmd
import vector_plotter as vp
import process_mjmag_data as mj


"""

    Class to interface with Luke's code to create vizualization
    of the magnetic field

    -- KG 06/18/19


 """

def main():
    # day = '083018'#''051917'#
    day ='061219'
    shot = day+'r13'#'40'#
    # t0 = 35
    # tf = 150

    # print("\n\n\n\t here?")
    time,bdot,timeb,b,bmod = mj.process_mjmag_data(shot)
    # time = pmd.get_time(shot)
    # data = pmd.get_data(shot)


    sample_Freq = 10 # sampling frequency


    b = b[:,:,t0:tf]
    s = np.arange(20) * 1.5# - 2.86
    x = b[0,0:20,:].T[0:-1:sample_Freq]
    y = b[1,0:20,:].T[0:-1:sample_Freq]
    t = timeb[t0:tf][0:-1:sample_Freq]

    path = os.getcwd()
    path = path + 'data/2019/'+day+'/'

    print(path)

    ############### Plotting the Bx and By vector rendering  ##############################
    # B = BField_Plotter.BField_xy_Animator(shot,t,s,x,y) #Instantiate the BField_xy_Animator Object
    B = BField_Plotter.BField_mag_Animator(shot,t,s,x,y) #Instantiate the BField_xyz_Animator Object
    # B.gen_Arrows(100) #Plot the vectors at some time, in this case t[10]
    animation = B.make_animation() #Now create an animation
    B.save_animation(path,animation) #Save the animation
    plt.show()

if __name__ == '__main__':
    main()
