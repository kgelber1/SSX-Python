import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib import animation
from timeit import default_timer
from scipy import interpolate
from matplotlib import gridspec

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs
    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

class BField_mag_Animator(object):

    def __init__(self,title,t,s,x,y):
        self.title = title
        self.t = t
        self.s = s
        probes = np.arange(len(s))
        xprobes = np.delete(probes,[3])
        yprobes = np.delete(probes,[12])

        for i in range(0,len(x)):
            fx = interpolate.interp1d(s[xprobes],x[i,xprobes],kind = 'cubic')
            fy = interpolate.interp1d(s[yprobes],y[i,yprobes],kind = 'cubic')
            x[i,:]= fx(s)
            y[i,:]= fy(s)

        self.x = x
        self.y = y
        self.zeros = np.zeros(len(s))

    #Set up the figure
    def __set_fig(self,i):
        plt.close('all')
        fig = plt.figure(figsize=(10,8),facecolor = 'white')
#        plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#        plt.rc('text', usetex=True)
#        plt.rc('lines', linewidth=2, color='k')
#        plt.rcParams['lines.linewidth'] = 1
        gs = gridspec.GridSpec(1,1)
        ax = fig.add_subplot(gs[0], projection='3d')
        ax.set_xlabel('Flux conserver axis ($cm$)',fontsize = 20)
        ax.set_ylabel('$B_x$ (Gauss)',fontsize = 20)
        ax.set_zlabel('$B_y$ (Gauss)',fontsize = 20)
        plt.suptitle('Shot ' + self.title+'\n Time = %.2f $\mu s$'%(self.t[i]),fontsize = 20, weight = 'bold')
        ax.grid(True)

        x = self.x
        y = self.y
        tol = 2
        ax.set_xlim3d(left=min(self.s) - tol, right=max(self.s)+tol, emit=True, auto=False)
        ax.set_ylim3d(bottom=np.amin(x), top=np.amax(x), emit=True, auto=False)
        ax.set_zlim3d(bottom=np.amin(y), top=np.amax(y), emit=True, auto=False)
        ax.view_init(elev=16.0, azim=-43.0)

        #Plot the probe axis
        #The coordinates from the SSX map to the python plots as xyz --> zxy
        ax.plot([min(self.s)-tol,max(self.s) +tol], [0,0], [0,0], color='black', alpha=0.7, lw=3) #plot probe axis

        #Plot back of can
        ax.plot([0,0], [0,0],[np.amin(x)-tol,np.amax(x) + tol], color='black', alpha=0.7, lw=2)
        ax.plot([0,0], [np.amin(x)-tol,np.amax(x) + tol],[0,0], color='black', alpha=0.7, lw=2)

        gs.tight_layout(fig)
        self.fig = fig
        self.ax = ax

    def gen_Arrows(self,index):
        self.__set_fig(index)
        ax = self.ax
        B_arrows = []
        s = self.s
        x = self.x
        y = self.y
        #vectors
        for i in range(0,len(s)):
            B_arrows.append(Arrow3D([s[i],s[i]], [0,x[index,i]], [0,y[index,i]], mutation_scale=20, lw=2, arrowstyle="-|>", color="c"))
            ax.add_artist(B_arrows[i])
        #print np.shape(x[i,:])
        tck, u = interpolate.splprep([s[2:],x[index,2:],y[index,2:]], k=1)
        u_fine = np.linspace(0,1,100)
        x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)
#        magplot, = ax.plot(x_fine, y_fine, z_fine, '--', color="purple",lw = 2)
        magplot, = ax.plot(x_fine, y_fine, z_fine, linestyle = '--', color="blue", alpha = .3,lw = 2)
        self.B_arrows = B_arrows
        #plt.show()
        return magplot,

# animation function.  This is called sequentially
    def __animate(self,i,magplot):
        ax = self.ax
        B_arrows = self.B_arrows
        s = self.s
        x = self.x
        y = self.y
        for j in range(0,len(s)):
            B_arrows[j].remove()
            B_arrows[j] = Arrow3D([s[j],s[j]], [0,x[i,j]], [0,y[i,j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="red")
            ax.add_artist(B_arrows[j])
            plt.suptitle('Shot ' + self.title+'\n Time = %.2f $\mu s$'%(self.t[i]),fontsize = 30)

        tck, u = interpolate.splprep([s[2:],x[i,2:],y[i,2:]], k=1)
        u_fine = np.linspace(0,1,100)
        x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)
        magplot.set_data(x_fine, y_fine)
        magplot.set_3d_properties(z_fine)
        self.B_arrows = B_arrows
        return magplot,

    def make_animation(self):
        magplot = self.gen_Arrows(0)
        #self.gen_Arrows(0)
        fig = self.fig
        anim = animation.FuncAnimation(fig, self.__animate, frames=len(self.t), interval=20, blit=False,fargs=(magplot))
        #plt.show()
        return anim

    def save_animation(self,path,anim):
        dpi = 250
        #print ('dpi = %i'%(dpi))
        writer = animation.writers['ffmpeg'](fps=30)
        start = default_timer(); print("Creating B Field movie for "+self.title)
        anim.save('B_VecRendering_'+self.title +'.mp4',writer=writer,dpi=dpi)
        duration = (default_timer() - start)/60; print("time: %s min" % str(duration))
