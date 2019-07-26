import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib import animation
from timeit import default_timer
from scipy import interpolate
#import smooth as sm
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

class BField_xy_Animator(object):

    def __init__(self,title,t,s,x,y):
        self.title = title
        self.t = t
        self.s = s
        probes = np.arange(len(s))
        # xprobes = np.delete(probes, [6,14,15])
        # yprobes = np.delete(probes, [11,12])
        xprobes = probes
        yprobes = probes
        for i in range(0,len(x)):
            fx = interpolate.interp1d(s[xprobes], x[i, xprobes], kind = 'cubic')
            fy = interpolate.interp1d(s[yprobes], y[i, yprobes], kind = 'cubic')
            x[i,:] = fx(s)
            y[i,:] = fy(s)
        self.x = x
        self.y = y
        self.zeros = np.zeros(len(s))

    #Set up the figure
    def __set_fig(self,i):
        plt.close('all')
        fig = plt.figure(figsize=(10,8),facecolor = 'white')
#        ax = fig.add_subplot(111, projection='3d')
        gs = gridspec.GridSpec(1,1)
#        plt.rc('font', **{'family':'serif', 'serif':['Computer Modern']})
#        plt.rc('text', usetex = True)
#        plt.rc('lines', linewidth = 2, color = 'k')
#        plt.rcParams['lines.linewidth'] = 1

        ax = fig.add_subplot(gs[0], projection='3d')
        ax.set_xlabel('Flux conserver axis ($cm$)',fontsize = 20)
        ax.set_ylabel('$B_x$ (Gauss)',fontsize = 20)
        ax.set_zlabel('$B_y$ (Gauss)',fontsize = 20)
        plt.suptitle('Shot ' + self.title+'\n Time = %.2f $\mu s$'%(self.t[i]),fontsize = 20, weight = 'bold')
        ax.grid(True)

        x = self.x
        y = self.y
        tol = 2
        ax.set_xlim3d(left=min(self.s) - tol, right=max(self.s) + tol, emit=True, auto=False)
        ax.set_ylim3d(bottom=np.amin(x), top=np.amax(x), emit=True, auto=False)
        ax.set_zlim3d(bottom=np.amin(y), top=np.amax(y), emit=True, auto=False)

        #Plot the probe axis
        #The coordinates from the SSX map to the python plots as xyz --> zxy
        ax.plot([min(self.s)-tol,max(self.s) +tol], [0,0], [0,0], color='black', alpha=0.8, lw=3) #plot probe axis
        #### plotting the back end of SFC ########
        ax.plot([0,0], [0,0], [np.amin(x)-tol, np.amax(x)+tol], color  = 'k', alpha =0.8, lw =2)
        ax.plot([0,0], [np.amin(x)-tol, np.amax(x)+tol], [0,0], color  = 'k', alpha =0.8, lw =2)
        self.fig = fig
        self.ax = ax

    def gen_Arrows(self,index):
        self.__set_fig(index)
        ax = self.ax
        Bx_arrows = []
        By_arrows = []
        s = self.s
        x = self.x
        y = self.y
        #vectors
        for i in range(0,len(s)):
            Bx_arrows.append(Arrow3D([s[i],s[i]], [0,x[index,i]], [0,0],mutation_scale=20,lw=2,arrowstyle="-|>",color="r"))
            By_arrows.append(Arrow3D([s[i],s[i]], [0,0], [0,y[index,i]], mutation_scale=20, lw=2, arrowstyle="-|>", color="b"))
            ax.add_artist(Bx_arrows[i])
            ax.add_artist(By_arrows[i])
        self.Bx_arrows = Bx_arrows
        self.By_arrows = By_arrows
        #plt.show()

# animation function.  This is called sequentially
    def __animate(self,i):
        ax = self.ax
        Bx_arrows = self.Bx_arrows
        By_arrows = self.By_arrows
        s = self.s
#        t = self.t
        x = self.x
        y = self.y
        for j in range(0,len(s)):
            Bx_arrows[j].remove()
            By_arrows[j].remove()
            Bx_arrows[j] = Arrow3D([s[j],s[j]], [0,x[i,j]], [0,0], mutation_scale=20, lw=2, arrowstyle="-|>", color="r")
            By_arrows[j] = Arrow3D([s[j],s[j]], [0,0], [0,y[i,j]], mutation_scale=20, lw=2, arrowstyle="-|>", color="b")
            ax.add_artist(Bx_arrows[j])
            ax.add_artist(By_arrows[j])
            plt.suptitle('Shot '+self.title+'\n Time = %.2f $\mu s$'%(self.t[i]),fontsize = 30)
            #ax.set_title('$B_x$  and $B_y$ at t = %.2f $\mu s$'%(t[i]),loc = 'center')
        self.Bx_arrows = Bx_arrows
        self.By_arrows = By_arrows

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
        print("Creating Bxy Field movie for "+self.title)
#        anim.save(path+ '(dpi='+str(dpi)+').avi',writer=writer,dpi=dpi)
        anim.save('Bxy_VecRendering_'+self.title+'.mp4',writer=writer,dpi=dpi)
#        anim.save('_BXY_.mp4',writer=writer,dpi=dpi)
        duration = (default_timer() - start)/60
        print "time: %s min" % str(duration)
