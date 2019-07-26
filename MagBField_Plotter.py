import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib import animation
from timeit import default_timer
from scipy import interpolate

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
        self.title=title
        self.t = t
        self.s = s

        # print("\n\n\n\nall good in inti")
        for i in range(0,len(x)):
            fx = interpolate.interp1d(s[np.arange(len(s))!=3],x[i,np.arange(len(s))!=3],kind = 'slinear')
            fy = interpolate.interp1d(s[np.arange(len(s))!=12],y[i,np.arange(len(s))!=12],kind = 'slinear')
            x[i,:]= fx(s)
            y[i,:]= fy(s)

        self.x = x
        self.y = y
        self.zeros = np.zeros(len(s))

    #Set up the figure
    def __set_fig(self,i):
        plt.close('all')
        fig = plt.figure(figsize=(12,10))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('Flux Conserver Axis (cm)')
        ax.set_ylabel('$B_x$ (Gauss)')
        ax.set_zlabel('$B_y$ (Gauss)')
        ax.set_title('$|B|$ at t = %.2f $\mu s$'%(self.t[i]),loc = 'center')
        ax.grid(True)

        x = self.x
        y = self.y
        tol = 2
        ax.set_xlim3d(left=min(self.s) - tol, right=max(self.s)+tol, emit=True, auto=False)
        ax.set_ylim3d(bottom=np.amin(x), top=np.amax(x), emit=True, auto=False)
        ax.set_zlim3d(bottom=np.amin(y), top=np.amax(y), emit=True, auto=False)
#        ax.view_init(elev=7.0, azim=-8.0)

        #Plot the probe axis
        #The coordinates from the SSX map to the python plots as xyz --> zxy
        ax.plot([min(self.s)-tol,max(self.s) +tol], [0,0], [0,0], color='black', alpha=0.8, lw=3) #plot probe axis
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
            B_arrows.append(Arrow3D([s[i],s[i]], [0,x[index,i]], [0,y[index,i]], mutation_scale=20, lw=2, arrowstyle="-|>", color="b"))
            ax.add_artist(B_arrows[i])
#        tck, u = interpolate.splprep([s,x[index,:],y[index,:]],k=1)
#        u_fine = np.linspace(0,1,100)
#        x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)
#        magplot, = ax.plot(x_fine, y_fine, z_fine, 'g',lw = 2)
        self.B_arrows = B_arrows
        plt.show()
        return

# animation function.  This is called sequentially
    def __animate(self,i):
        ax = self.ax
        B_arrows = self.B_arrows
        s = self.s
        t = self.t
        x = self.x
        y = self.y
        for j in range(0,len(s)):
            B_arrows[j].remove()
            B_arrows[j] = Arrow3D([s[j],s[j]], [0,x[i,j]], [0,y[i,j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="b")
            ax.add_artist(B_arrows[j])
            ax.set_title('$|B|$ at t = %.2f $\mu s$'%(t[i]),loc = 'center')

        tck, u = interpolate.splprep([s,x[i,:],y[i,:]])
        u_fine = np.linspace(0,1,100)
        x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)
#        magplot.set_data(x_fine, y_fine)
#        magplot.set_3d_properties(z_fine)
        self.B_arrows = B_arrows
        return

    def make_animation(self):
        self.gen_Arrows(0)
        fig = self.fig
        anim = animation.FuncAnimation(fig, self.__animate, frames=len(self.t), interval=20, blit=False)
        plt.show()
        return anim

    def save_animation(self,path,anim):
        dpi = 175
        writer = animation.writers['ffmpeg'](fps=30)
        start = default_timer()
        print("Creating BField Magnitude movie...")
        anim.save(self.title+'_BField.mp4',writer=writer,dpi=dpi)
        duration = default_timer() - start
        print "time: %s s" % str(duration)
