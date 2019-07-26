import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots(1,1)
x=np.linspace(np.pi,4*np.pi,100)
N=len(x)
ax.set_xlim(len(x))
ax.set_ylim(-1.5,1.5)
line,  = ax.plot([],[],'o-')

def init():
    line.set_ydata(np.ma.array(x[:], mask=True))
    return line,

def animate(i, *args, **kwargs):
    y=np.sin(x*i)
    line.set_data(np.arange(N),y)            # update the data
    return line,

ani = animation.FuncAnimation(fig, animate, init_func=init,
     frames=100, interval=10, blit= False, repeat = False)
ani.save('2osc.mp4', writer="ffmpeg")
fig.show()
