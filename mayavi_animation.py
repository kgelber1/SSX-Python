# import numpy as np
# import mayavi.mlab as mlab
# # from enthought.mayavi import mlab
# # from enthought.mayavi.mlab import quiver3d,savefig,view,clf
#
# import  moviepy.editor as mpy
#
# duration= 2 # duration of the animation in seconds (it will loop)
#
# # MAKE A FIGURE WITH MAYAVI
#
# fig_myv = mlab.figure(size=(220,220), bgcolor=(1,1,1))
# X, Y = np.linspace(-2,2,200), np.linspace(-2,2,200)
# XX, YY = np.meshgrid(X,Y)
# ZZ = lambda d: np.sinc(XX**2+YY**2)+np.sin(XX+d)
#
# # ANIMATE THE FIGURE WITH MOVIEPY, WRITE AN ANIMATED GIF
#
# def make_frame(t):
#     mlab.clf() # clear the figure (to reset the colors)
#     mlab.mesh(YY,XX,ZZ(2*np.pi*t/duration), figure=fig_myv)
#     return mlab.screenshot(antialiased=True)
#
# animation = mpy.VideoClip(make_frame, duration=duration)
# animation.write_mp4("sinc.gif", fps=20)


from mayavi import mlab
import numpy as np

def V(x, y, z):
    """ A 3D sinusoidal lattice with a parabolic confinement. """
    return np.cos(10*x) + np.cos(10*y) + np.cos(10*z) + 2*(x**2 + y**2 + z**2)
X, Y, Z = np.mgrid[-2:2:100j, -2:2:100j, -2:2:100j]
mlab.contour3d(X, Y, Z, V)
