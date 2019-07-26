# import urllib
# import numpy as np
# from scipy.ndimage.filters import convolve
# import moviepy.editor as mpy
# import vector_plotter as vp
#
#
#
# ##### MODEL
#
# def update(world):
#     """ spread the epidemic for one time step """
#     infect = infection(world['SIR'], infection_rate, incubation_rate)
#     disperse = dispersion(world['SIR'], dispersion_kernel, dispersion_rates)
#     world['SIR'] += dt*( infect + disperse)
#     world['t'] += dt
#
#
# # ANIMATION WITH MOVIEPY
#
#
# def world_to_npimage(world):
#     """ Converts the world's map into a RGB image for the final video."""
#     coefs = np.array([2,25,25]).reshape((3,1,1))
#     accentuated_world = 255*coefs*world['SIR']
#     image = accentuated_world[::-1].swapaxes(0,2).swapaxes(0,1)
#     return np.minimum(255, image)
#
# def make_frame(t):
#     """ Return the frame for time t """
#     while world['t'] < hours_per_second*t:
#         update(world)
#     return world_to_npimage(world)
#
#
# animation = mpy.VideoClip(make_frame, duration=25)
# # You can write the result as a gif (veeery slow) or a video:
# #animation.write_gif(make_frame, fps=15)
# animation.write_videofile('test.mp4', fps=20)


import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz
from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy
import vector_plotter as vp
import os

# DRAW A FIGURE WITH MATPLOTLIB

duration = 2
#

day = '061219r'
shot_num ='8'
# day = '050119r'
# shot_num = 3
shot = day+str(shot_num)

fig_mpl, ax = plt.subplots(1,figsize=(5,3), facecolor='white')
# data = [[0] *3 for i in range(num_probes)]

# ax.set_title("Testing")
# ax.set_ylim(-1.5,2.5)
# line, = ax.plot(times,Bx[0], lw=3)
ax, qq= vp.get_color_frame(shot, 0)
cbar = plt.colorbar(qq, cmap=plt.cm.jet)
cbar.set_label("Magnetic Field Strength (T)")
plt.xlabel("Y direction")
plt.ylabel("Z direction")
plt.xlim(-.04,.04)
plt.ylim(-.04,.04)
plt.title("Magnetic Data - Shot %s" %(shot))
# ANIMATE WITH MOVIEPY (UPDATE THE CURVE FOR EACH t). MAKE A GIF.

def make_frame_mpl(t):
    plt_new, qq = vp.get_color_frame(shot, t)# <= Update the curve
    ax = plt_new
    return mplfig_to_npimage(fig_mpl) # RGB image of the figure

def getFrame(t):
    day = '061219r'
    shot_num ='8'
    # day = '050119r'
    # shot_num = 3
    shot = day+str(shot_num)

    return vp.get_color_frame(shot, t)

animation =mpy.VideoClip(make_frame_mpl, duration=duration)
animation.write_gif("sinc_mpl.gif", fps=20)
