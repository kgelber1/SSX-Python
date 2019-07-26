from __future__ import division, print_function, absolute_import

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
from matplotlib.text import TextPath
from mpl_toolkits.mplot3d import proj3d
from matplotlib import animation
from timeit import default_timer
from scipy import interpolate
from matplotlib import gridspec
import os

import process_mjmag_data as mj
import iter_smooth as ism
import mytools as my
import ssx_data_read_basic as sdr
import process_mag_data_4x4 as pmd
import vector_plotter as vp
from  nBT import *

"""
This one creates an animation of B-field in one subplot, and then also plots
the density and temperature in two other subplots below

Based on:
 - Luke's code
 - https://stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d
 - https://matplotlib.org/gallery/mplot3d/pathpatch3d.html#sphx-glr-gallery-mplot3d-pathpatch3d-py

-xyz refers to the probe's coordinate system
-rzt refers to the coordinate system used by the SSX

j is usually the index corresponding to probe number (in this file at least)

-- KG 06/27/19
"""
##############################################
#shot = '020618r2'


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        """Starts as 2-d,g ets projected into 3-d"""
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


#
# class Circle3D(Circle):
#     def __init__(self, xs, ys, zs, *args, **kwargs):
#         Circle.__init__(self, *args, **kwargs)
#         # Circle.set_radius(self, rad)
#         self._verts3d = xs, ys, zs
#
#     def draw(self, renderer):
#         """Starts as 2-d,g ets projected into 3-d"""
#         xs3d, ys3d, zs3d = self._verts3d
#         xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
#         # self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
#         Circle.draw(self, renderer)

class BField_mag_Animator(object):

    def __init__(self, title,t,probe_locs,x,y,z, num_probes = 16):
        self.title = title
        self.Plot_title = title
        self.t = t
        self.num_probes = num_probes
        self.flags = []
        self.dead = []

        """Change between the calibration setup and the SSX setup...
             x -> z
             y -> theta
             z -> -r

             scale the r positions, and then scale the z,theta values

             r, z, theta will now need to be indexed using [i][j]

             -- KG 06/25/19"""
        max_pos = np.amax(probe_locs)

        x_min = np.amin(x)
        y_min = np.amin(y)
        z_min = np.amin(z)
        min = np.min([x_min,y_min,z_min])
        self.theta_min = min

        x_max = np.amax(x)
        y_max = np.amax(y)
        z_max = np.amax(z)
        max   = np.max([x_max,y_max,z_max])
        self.theta_max = max

        max = np.max(np.absolute([max,min]))
        # print(max)


        #scale everything between the ssx z-axis and the max B-vector
        scale = np.amax(probe_locs[:][0])/max
        self.scale = scale

        self.z_pos =[probe_locs[j][0] for j in range(num_probes)]
        self.theta_pos = [probe_locs[j][1]*(1/scale) for j in range(num_probes)]
        self.r_pos = [probe_locs[j][2]*-1 for j in range(num_probes)]

        self.r = [[zr[j]*scale*-1+self.r_pos[j] for j in range(len(zr))] for zr in z]
        self.z = [[xz[j]*scale+self.z_pos[j] for j in range(len(xz))] for xz in x]
        self.theta = [[yt[j]+self.theta_pos[j] for j in range(len(yt))] for yt in y]


    #Set up the figure
    def _set_Plot_title(self, title):
        self.Plot_title = self.title + title
    def _set_flags(self, flags):
        self.flags = flags
    def _set_dead(self, dead):
        self.dead = dead


    def __set_fig(self,i):
        plt.close('all')
        fig = plt.figure(figsize=(10,8),facecolor = 'white')
#        plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#        plt.rc('text', usetex=True)
#        plt.rc('lines', linewidth=2, color='k')
#        plt.rcParams['lines.linewidth'] = 1
        # gs = gridspec.GridSpec(1,1)
        gs = gridspec.GridSpec(5, 2)
        ax = fig.add_subplot(gs[0:4, :], projection='3d')
        ax.set_xlabel('\nProbe Axis - ($cm$)\n$B_r$',fontsize = 20)
        ax.set_ylabel('Flux conserver axis ($cm$)',fontsize = 20)
        ax.set_zlabel('$B_\\theta$ (Gauss)',fontsize = 20)
        plt.suptitle('Shot ' + self.title+'\n Time = %.2f $\mu s$'%(self.t[i]),fontsize = 20, weight = 'bold')
        ax.grid(True)

        r = self.r
        z = self.z
        theta = self.theta

        z_min = np.amin(self.z_pos)
        z_max = np.amax(self.z_pos)
        r_min = np.amin(self.r_pos)
        r_max = np.amax(self.r_pos)
        theta_min = self.theta_min
        theta_max = self.theta_max
        ax.set_xlim3d(left=1.1*r_min, right=1.1*r_max, emit=True, auto=False)
        ax.set_ylim3d(bottom=1.1*z_min, top=1.1*z_max, emit=True, auto=False)
        ax.set_zlim3d(bottom=1*theta_min, top=1*theta_max, emit=True, auto=False)

        # original
        # ax.view_init(elev=16.0, azim=-43.0)
        # ax.view_init(elev=30.0, azim=-20.0)
        ax.view_init(elev=26.0, azim=-25.0)

        ax1 = fig.add_subplot(gs[4, 0])
        plt.text(0.07,0.92,'(b)',fontsize = 18, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,)

        ax2 = fig.add_subplot(gs[4, 1])
        ax2.text(0.07,0.92,'(c)',fontsize = 18, weight='bold',horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes)

        ax1, ax2, self.ne, self.ti = get_nT(self.title, self.t, ax1, ax2)
        self.ax1 = ax1
        self.ax2 = ax2

        #Plot the probe axis
        #The coordinates from the SSX map to the python plots as xyz --> zxy
        for i in range(0, self.num_probes, 4):
            ax.plot([self.r_pos[i], self.r_pos[i+3]],[self.z_pos[i], self.z_pos[i]], [0,0] ,color='k', alpha=0.4, lw=3)

        # plot centerline
        ax.plot([0,0], [1*z_min,1.3*z_max],[0,0], color='black', alpha=0.7, lw=2, linestyle = "dashed")


        #TODO ADD OTHER OVERLAYS HERE (if you want)


        #add the semicircles
        phi = np.linspace(np.pi/6, np.pi*(5/6), 100)
        r = (16.05 + .63)/2
        x = r * np.sin(phi)
        y = r * np.cos(phi)/ self.scale
        # scale = 1/self.scale
        # radius = (16.05 + .63)/2
        cz = 17.78/2

        ax.plot(x,  [-cz for xx in x], y, label='East', color = "k", linestyle = "dashed", alpha = .3)
        ax.plot(x,  [cz for xx in x], y, label='West', color = "k", linestyle = "dashed", alpha = .3)

        ax.text(.9*r_min, z_min*1.1, .6*theta_max, "East", color='k')
        ax.text(.4*r_max, z_max*.9, .7*theta_max, "West", color='k')

        gs.tight_layout(fig)
        self.fig = fig
        self.ax = ax

    def gen_Arrows(self,i, num_probes = 16):

        self.__set_fig(i)
        ax = self.ax
        ax1 = self.ax1
        ax2 = self.ax2
        B_arrows = []

        r_pos = self.r_pos
        z_pos = self.z_pos
        theta_pos = self.theta_pos
        r = self.r
        z = self.z
        theta = self.theta
        flags = self.flags
        dead = self.dead
        z_min = np.amin(self.z_pos)
        z_max = np.amax(self.z_pos)


        for j in range(num_probes):
            if j in flags:
                B_arrows.append(Arrow3D([r_pos[j],r[i][j]], [z_pos[j],z[i][j]], [theta_pos[j], theta[i][j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="steelblue", alpha = .9))
            elif j in dead:
                B_arrows.append(Arrow3D([r_pos[j],r[i][j]], [z_pos[j],z[i][j]], [theta_pos[j], theta[i][j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="darkviolet", alpha = .9))
            else:
                B_arrows.append(Arrow3D([r_pos[j],r[i][j]], [z_pos[j],z[i][j]], [theta_pos[j], theta[i][j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="mediumblue", alpha = .9))
            ax.add_artist(B_arrows[j])
        """Leaving this in in case future people want to add interpolation back"""
        # magplot, = ax.plot([t[i]]*2, [ne[i], ti[i]], )
        # tck, u = interpolate.splprep([x[i,:],y[i,:],z[i:],], k=1)
        # tck, u = interpolate.splprep([r_pos[2:],r[i][2:],z[i][2:]], k=1)
        # u_fine = np.linspace(0,1,100)
        # r_fine, z_fine, t_fine = interpolate.splev(u_fine, tck)
#        magplot, = ax.plot(x_fine, y_fine, z_fine, '--', color="purple",lw = 2)
        # magplot, = ax.plot(r_fine, z_fine, t_fine, linestyle = '--', color="blue", alpha = 0,lw = 2)
        # ax1.add_artist[circle]


        magplot, = ax.plot([0,0], [1*z_min,1.3*z_max],[0,0], color='black', alpha=0.7, lw=2, linestyle = "dashed")
        self.B_arrows = B_arrows
        return magplot,


# animation function.  This is called sequentially
    def __animate(self,i,magplot):
        """ Update function """
        ax = self.ax
        B_arrows = self.B_arrows

        r_pos = self.r_pos
        z_pos = self.z_pos
        theta_pos = self.theta_pos
        r = self.r
        z = self.z
        theta = self.theta
        flags = self.flags
        dead = self.dead
        z_min = np.amin(self.z_pos)
        z_max = np.amax(self.z_pos)

        # dot = Circle3D
        for j in range(self.num_probes):
            #delete old frame
            B_arrows[j].remove()
            if j in flags:
                B_arrows[j] = Arrow3D([r_pos[j],r[i][j]],[z_pos[j],z[i][j]], [theta_pos[j],theta[i][j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="steelblue", alpha = .9)
            elif j in dead:
                B_arrows[j] = Arrow3D([r_pos[j],r[i][j]],[z_pos[j],z[i][j]], [theta_pos[j],theta[i][j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="darkviolet", alpha = .9)
            else:
                B_arrows[j] = Arrow3D([r_pos[j],r[i][j]],[z_pos[j],z[i][j]], [theta_pos[j],theta[i][j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="mediumblue", alpha = .9)
            ax.add_artist(B_arrows[j])
            plt.suptitle('Shot ' + self.Plot_title+'\n Time = %.2f $\mu s$'%(self.t[i]),fontsize = 30)

        ax1 = self.ax1
        ax2 = self.ax2
        # del ax1
        ax1.scatter(self.t[i], self.ne[i], color = 'red', alpha = .9)
        # del ax2
        ax2.scatter(self.t[i], self.ti[i], color = 'red', alpha = .9)
        # ax1.scatter(self.t[i], self.ti[i], color = 'red')
        """Leaving this in in case future people want to add interpolation back"""
        # tck, u = interpolate.splprep([r_pos[2:],r[i][2:],z[i][2:]], k=1)
        # u_fine = np.linspace(0,1,100)
        # r_fine, z_fine, theta_fine = interpolate.splev(u_fine, tck)
        # magplot.set_data(r_fine, z_fine)
        # magplot.set_3d_properties(theta_fine)
        # ax.plot()

        self.B_arrows = B_arrows
        return ax1,

    def make_animation(self):
        magplot = self.gen_Arrows(0)
        #self.gen_Arrows(0)
        fig = self.fig
        anim = animation.FuncAnimation(fig, self.__animate, frames=len(self.t), interval=20, blit=False,fargs=(magplot))
        # anim = animation.FuncAnimation(fig, self.__animate, frames=len(self.t), init_func = genA interval=20, blit=False)
        return anim




    def save_animation(self,path,anim):
        """Function to save the animation itsself.

           This function takes a bit longer, and needs a lot of cache (or else it
           is impossibly slow.) I recommend closing all other tasks on your computer
           if it is going very slowly, or turning up the sampeling frequency.

           --KG 07/15/19"""

        # resolution of saved graphic
        dpi = 250
        print ('dpi = %i'%(dpi))

        start = default_timer(); print("Creating B Field movie for "+self.title)
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        filename = path + 'B_VecRendering_'+self.title +'.mp4'
        anim.save(filename,'ffmpeg')
        duration = (default_timer() - start)/60; print("time: %s min" % str(duration))
        print("\n\nFile saved as: %s\n\n" %filename)
        # print
        # anim.save(filename, writer="ffmpeg")


def run(day, shot, t0 = 25, tf = 75, sample_Freq = 5, show = True):
    """ Main animating function.

        Update paramters like day and shot, as well as start and end times



        To make this code work with a different set up, you will need to:
         - update the calibration used in pmd (process_mag_data_4x4)
         - create new probe locations function (vp.get_probeLocs).
         That's it!


         --KG 07/15/19
        """
    # day = '061219'#'062618'#'101917'
    # shot = 13#85


    shot = day+'r'+ str(shot) #'40'#
    timeb,b = pmd.get_Bxy_vecRendering_format_lookup(shot)

    ######### fixing bounds using a function called fix_bounds ##########
    t_index, t = my.fix_bounds(np.array((t0, tf)), timeb)

    #Find what start and end index correspond to the right times
    t0 = t_index[0]
    tf = t_index[-1]
    b = b[:,:,t0:tf]


    probe_locs = vp.get_probeLocs_SSX_setup_cm(num_probes = 16)

    #convert things to time vs probes instead of probes vs time
    #and also thin the time
    x = b[0,:,:].T[0:-1:sample_Freq]
    y = b[1,:,:].T[0:-1:sample_Freq]
    z = b[2,:,:].T[0:-1:sample_Freq]

    t = timeb[t0:tf][0:-1:sample_Freq]
    # path = os.getcwd() + '\\' + 'Magnetic_animations\\'
    path = os.getcwd() + '\\data\\2019\\' +day+'\\Analyzed\\'

    ##############################################
    ######## Plotting the magnetic field vector rendering ############

    Bmag = BField_mag_Animator(shot,t,probe_locs,x,y,z)#Instantiate the BField_xy_Animator Object
    Bmag._set_flags([8,10,12]) #some probes have questionable directions
    Bmag._set_dead([0])#probe has dead z- component
    # Bmag._set_Plot_title("\nWest Gun lagging by 1 $\\mu$s")
    # Bmag._set_Plot_title("\nWest Gun leading by 5 $\\mu$s")
    animat = Bmag.make_animation() #Now create an animation
    # animation.save('test.mp4', writer="ffmpeg")
    Bmag.save_animation(path,animat) #Save the animation
    if show:
        plt.show()

def main():
    """Just a place to specifiy variables"""
    day ='072419'
    shot = 26

    sample_Freq = 5# sampling frequency - turn up for faster animations
    t0 = 25
    tf = 55

    run(day, shot, t0, tf, sample_Freq, show = True)



if __name__ == '__main__':
    main()
