from __future__ import division, print_function, absolute_import

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import proj3d
from matplotlib import animation
from timeit import default_timer
from scipy import interpolate
from matplotlib import gridspec
import os

# reload(mj)
import process_mjmag_data as mj
import iter_smooth as ism
import mytools as my
import ssx_data_read_basic as sdr
import process_mag_data_4x4 as pmd
import vector_plotter as vp

"""
Code to create a 3-d animation of the magetic field
Based on Luke's code and https://stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d

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

class BField_mag_Animator(object):

    def __init__(self,title,t,probe_locs,x,y,z, dir =2, num_probes = 16):
        self.title = title
        self.Plot_title = title
        self.t = t
        self.num_probes = num_probes
        self.dir = dir
        # for i in range(0,len(x)):
        #     fx = interpolate.interp1d(probe_locs[xprobes],x[i,xprobes],kind = 'cubic')
        #     fy = interpolate.interp1d(probe_locs[yprobes],y[i,yprobes],kind = 'cubic')
        #     x[i,:]= fx(s)
        #     y[i,:]= fy(s)

        """There is no change between the calibration setup...

             but x,y,z will now need to be indexed using [i][j]
             instead of [i,j]

             -- KG 07/09/19        """
        max_pos = np.amax(probe_locs)
        # print(max_pos)
        if dir ==  0:
            self.y_pos =[probe_locs[j][1] for j in range(num_probes)]
            self.ybound_up =  1.2*np.amax(self.y_pos)
            self.ybound_low = 1.2*np.amin(self.y_pos)
            yscale = self.ybound_up/np.amax(y)
            self.scale = yscale
            # print(self.scale)

            self.z_pos =[probe_locs[j][2] for j in range(num_probes)]
            self.zbound_up =  1.2*np.amax(self.z_pos)
            self.zbound_low = 1.2*np.amin(self.z_pos)
            zscale = self.zbound_up/np.amax(z)


            self.x_pos = [probe_locs[j][0]*(1/yscale) for j in range(num_probes)]


            self.x = [[xx[j]+self.x_pos[j] for j in range(len(xx))] for xx in x]
            self.y = [[yy[j]*yscale+self.y_pos[j] for j in range(len(yy))] for yy in y]
            self.z = [[zz[j]*zscale+self.z_pos[j] for j in range(len(zz))] for zz in z]


            # print("\n", x)
        elif dir ==  1:


            self.x_pos =[probe_locs[j][0] for j in range(num_probes)]
            self.xbound_up =  1.2*np.amax(self.x_pos)
            self.xbound_low = 1.2*np.amin(self.x_pos)
            xscale = self.xbound_up/np.amax(x)
            self.scale = xscale

            self.z_pos =[probe_locs[j][2] for j in range(num_probes)]
            self.zbound_up =  1.2*np.amax(self.z_pos)
            self.zbound_low = 1.2*np.amin(self.z_pos)
            zscale = self.zbound_up/np.amax(z)


            self.y_pos = [probe_locs[j][1]*(1/xscale) for j in range(num_probes)]

            # print(self.x_pos)
            # print(self.y_pos)
            self.x = [[xx[j]*xscale+self.x_pos[j] for j in range(len(xx))] for xx in x]
            self.y = [[yy[j]+self.y_pos[j] for j in range(len(yy))] for yy in y]
            self.z = [[zz[j]*zscale+self.z_pos[j] for j in range(len(zz))] for zz in z]

        else:

            self.x_pos =[probe_locs[j][0] for j in range(num_probes)]
            self.xbound_up =  1.2*np.amax(self.x_pos)
            self.xbound_low = 1.2*np.amin(self.x_pos)
            xscale = self.xbound_up/np.amax(x)

            self.y_pos =[probe_locs[j][1] for j in range(num_probes)]
            self.ybound_up =  1.2*np.amax(self.y_pos)
            self.ybound_low = 1.2*np.amin(self.y_pos)
            yscale = self.ybound_up/np.amax(y)
            self.scale = yscale

            self.z_pos = [probe_locs[j][2]*(1/yscale) for j in range(num_probes)]

            self.x = [[xx[j]*xscale+self.x_pos[j] for j in range(len(xx))] for xx in x]
            self.y = [[yy[j]*yscale+self.y_pos[j] for j in range(len(yy))] for yy in y]
            self.z = [[zz[j]+self.z_pos[j] for j in range(len(zz))] for zz in z]



    #Set up the figure
    def _set_Plot_title(self, title):
        self.Plot_title = self.title + title
    def __set_fig(self,i):
        plt.close('all')
        fig = plt.figure(figsize=(10,8),facecolor = 'white')
#        plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#        plt.rc('text', usetex=True)
#        plt.rc('lines', linewidth=2, color='k')
#        plt.rcParams['lines.linewidth'] = 1
        gs = gridspec.GridSpec(1,1)
        ax = fig.add_subplot(gs[0], projection='3d')

        plt.suptitle('Shot ' + self.title+'\n Time = %.2f $\mu s$'%(self.t[i]),fontsize = 20, weight = 'bold')
        ax.grid(True)

        x = self.x
        y = self.y
        z = self.z

        ax.view_init(elev=30.0, azim=-20.0)
        ax.view_init(elev=97.0, azim=0.0)
        # ax.view_init(elev=-85.0, azim=-90.0)




        # original
        # ax.view_init(elev=16.0, azim=-43.0)

        #plot the centerline
        if self.dir == 0:
            ax.plot([np.amin(x),np.amax(x)],[0,0], [0,0], color='black', alpha=0.7, lw=2, linestyle = "dashed")
            ax.set_xlim3d(left=np.amin(x), right=np.amax(x), emit=True, auto=False)
            ax.set_ylim3d(bottom=self.ybound_low, top=self.ybound_up, emit=True, auto=False)
            ax.set_zlim3d(bottom=self.zbound_low, top=self.zbound_up, emit=True, auto=False)

            ax.set_xlabel('Bx (Gauss)',fontsize = 20)
            ax.set_ylabel('Probe axis (cm)',fontsize = 20)
            ax.set_zlabel('Probe axis (cm)',fontsize = 20)

        elif self.dir == 1:
            ax.plot([0,0],[np.amin(y),np.amax(y)],  [0,0],color='black', alpha=0.7, lw=2, linestyle = "dashed")
            ax.set_xlim3d(left=self.xbound_low, right=self.xbound_up, emit=True, auto=False)
            ax.set_ylim3d(bottom=np.amin(y), top=np.amax(y), emit=True, auto=False)
            ax.set_zlim3d(bottom=self.zbound_low, top=self.zbound_up, emit=True, auto=False)

            ax.set_xlabel('Probe axis (cm)',fontsize = 20)
            ax.set_ylabel('By (Gauss)',fontsize = 20)
            ax.set_zlabel('Probe axis (cm)',fontsize = 20)

        else:
            ax.plot([0,0], [0,0],[np.amin(z),np.amax(z)], color='black', alpha=0.7, lw=2, linestyle = "dashed")
            ax.set_xlim3d(left=self.xbound_low, right=self.xbound_up, emit=True, auto=False)
            ax.set_ylim3d(bottom=self.ybound_low, top=self.ybound_up, emit=True, auto=False)
            ax.set_zlim3d(bottom=np.amin(z), top=np.amax(z), emit=True, auto=False)

            ax.set_xlabel('Probe axis (cm)',fontsize = 20)
            ax.set_ylabel('Probe axis (cm)',fontsize = 20)
            ax.set_zlabel('Bz (Gauss)',fontsize = 20)

        gs.tight_layout(fig)
        self.fig = fig
        self.ax = ax

    def gen_Arrows(self,i, num_probes = 16):
        # print(index)
        self.__set_fig(i)
        ax = self.ax
        B_arrows = []

        x_pos = self.x_pos
        y_pos = self.y_pos
        z_pos = self.z_pos
        x = self.x
        y = self.y
        z = self.z

        # i --> j
        # index --> i
        #vectors
        for j in range(num_probes):
            B_arrows.append(Arrow3D([x_pos[j],x[i][j]], [y_pos[j],y[i][j]], [z_pos[j], z[i][j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="blue", alpha = .9))
            ax.add_artist(B_arrows[j])

        if self.dir == 0:
            magplot, = ax.plot([np.amin(x),np.amax(x)],[0,0], [0,0], color='black', alpha=0.7, lw=2, linestyle = "dashed")
        elif self.dir == 1:
            magplot, = ax.plot([0,0],[np.amin(y),np.amax(y)],  [0,0],color='black', alpha=0.7, lw=2, linestyle = "dashed")
        else:
            magplot, = ax.plot([0,0], [0,0],[np.amin(z),np.amax(z)], color='black', alpha=0.7, lw=2, linestyle = "dashed")

        self.B_arrows = B_arrows
        #plt.show()
        return magplot,


# animation function.  This is called sequentially
    def __animate(self,i,magplot):
        """ Update function """
        ax = self.ax
        B_arrows = self.B_arrows

        x_pos = self.x_pos
        y_pos = self.y_pos
        z_pos = self.z_pos
        x = self.x
        y = self.y
        z = self.z
        # print(len(s))

        for j in range(self.num_probes):
            # ax.remove_artist(B_arrows[j])
            #delete old frame
            B_arrows[j].remove()
            B_arrows[j] = Arrow3D([x_pos[j],x[i][j]],[y_pos[j],y[i][j]], [z_pos[j],z[i][j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="blue", alpha = .9)
            ax.add_artist(B_arrows[j])
            plt.suptitle('Shot ' + self.Plot_title+'\n Time = %.2f $\mu s$'%(self.t[i]),fontsize = 30)

        # tck, u = interpolate.splprep([x_pos[2:],r[i][2:],z[i][2:]], k=1)
        # u_fine = np.linspace(0,1,100)
        # r_fine, z_fine, theta_fine = interpolate.splev(u_fine, tck)
        # magplot.set_data(r_fine, z_fine)
        # magplot.set_3d_properties(theta_fine)
        self.B_arrows = B_arrows
        if self.dir == 0:
            magplot, = ax.plot([np.amin(x),np.amax(x)],[0,0], [0,0], color='black', alpha=0.7, lw=2, linestyle = "dashed")
        elif self.dir == 1:
            magplot, = ax.plot([0,0],[np.amin(y),np.amax(y)],  [0,0],color='black', alpha=0.7, lw=2, linestyle = "dashed")
        else:
            magplot, = ax.plot([0,0], [0,0],[np.amin(z),np.amax(z)], color='black', alpha=0.7, lw=2, linestyle = "dashed")
        return magplot,

    def make_animation(self):
        magplot = self.gen_Arrows(0)
        #self.gen_Arrows(0)
        fig = self.fig
        anim = animation.FuncAnimation(fig, self.__animate, frames=len(self.t), interval=20, blit=False,fargs=(magplot))
        # anim = animation.FuncAnimation(fig, self.__animate, frames=len(self.t), init_func = genA interval=20, blit=False)
        #plt.show()
        return anim




    def save_animation(self,path,anim):
        dpi = 250
        print ('dpi = %i'%(dpi))

        start = default_timer(); print("Creating B Field movie for "+self.title)
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        filename = path + 'B_VecRendering_calib'+self.title +'.mp4'
        anim.save(filename,writer=writer,dpi=dpi)
        duration = (default_timer() - start)/60; print("time: %s min" % str(duration))
        print("\n\nFile saved as: %s\n\n" %filename)


def main():
    # day = '061219'#'062618'#'101917'
    # shot = 13#85

    day ='050119'
    shot = 23

    dir = 2

    if day == '050119':
        if (shot >= 16 and shot <=20):
            dir = 0
        elif (shot >= 21 and shot <=25):
            dir = 1
        elif(shot >= 11 and shot <=15):
            dir = 2
            # print(dir)
        else:
            print('Program will run as if it was in Z-dir. If this is a mistake,please check the shotnumber')


    shot = day+'r'+ str(shot) #'40'#

    timeb,b = pmd.get_Bxy_vecRendering_format(shot)
    # timeb=timeb-timeb[0]-2
    # timeb = time

    ######### fixing bounds using a function called fix_bounds ###########
    t0 = 5
    tf = 90
    t_index, t = my.fix_bounds(np.array((t0, tf)), timeb)

    #just trying to find what start and end index correspond to the right times
    t0 = t_index[0]
    tf = t_index[-1]
    sample_Freq = 20 # sampling frequency
    b = b[:,:,t0:tf]
    # b 

    # probe_locs = vp.get_probeLocs_SSX_setup(num_probes)
    probe_locs = vp.get_probeLocs_calib_setup(dir =2, num_probes = 16)

    #convert things to time x probes instead of probes x time
    #and also thin the time
    x = b[0,:,:].T[0:-1:sample_Freq]
    y = b[1,:,:].T[0:-1:sample_Freq]
    z = b[2,:,:].T[0:-1:sample_Freq]

    t = timeb[t0:tf][0:-1:sample_Freq]

    path = os.getcwd() + '\\' + 'Magnetic_animations\\'

    ##############################################
    ######## Plotting the magnetic field vector rendering ############

    Bmag = BField_mag_Animator(shot,t,probe_locs,x,y,z, dir)#Instantiate the BField_xy_Animator Object
    Bmag._set_Plot_title("\nCalibration in Y-dir")
    # Bmag._set_Plot_title("\nWest Gun leading by 5 $\\mu$s")
    animation = Bmag.make_animation() #Now create an animation
    # Bmag.save_animation(path,animation) #Save the animation

    plt.show()

if __name__ == '__main__':
    main()
