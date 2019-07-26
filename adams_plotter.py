def vector_movie_v0 (timestamp,tmin=50.0,tmax=250.0,nskip=100,integrate=False,showplot=False):



    from frc_data_class import FrcData,FrcDataOld

    from enthought.mayavi import mlab

    from enthought.mayavi.mlab import quiver3d,savefig,view,clf

    from scipy import arange,ndarray,array,zeros
    from scipy.integrate import cumtrapz

    import numpy



    filename = 'mag.x1.'+timestamp+'.csv'

    shot = FrcData(filename)

    if tmin > tmax:

        temp = tmin

        tmin = tmax

        tmax = temp

    imin = tmin*1e-6/shot.deltat

    imax = tmax*1e-6/shot.deltat

    n = imax-imin

    time = shot.timepoints[imin:imax]

    Bx = ndarray((16,n),float)

    By = ndarray((16,n),float)

    Bz = ndarray((16,n),float)



    Bdot = {}

    # for ch in arange(16):
    #
    #     for axis in ['x','y','z']:
    #
    #         filename = 'mag.'+axis+str(ch+1)+'.'+timestamp+'.csv'
    #
    #         shot = FrcData(filename)
    #
    #         dB = shot.voltage[imin:imax]*calfactor
    #
    #         Bdot[axis] = dB/numpy.max(abs(dB))
    #
    #
    #
    #     if integrate:
    #
    #         Bx[ch,1:] = cumtrapz(Bdot['x'],x=time)
    #
    #         By[ch,1:] = cumtrapz(Bdot['y'],x=time)
    #
    #         Bz[ch,1:] = cumtrapz(Bdot['z'],x=time)
    #
    #     else:
    #
    #         Bx[ch,:] = Bdot['x']
    #
    #         By[ch,:] = Bdot['y']
    #
    #         Bz[ch,:] = Bdot['z']



    modB = numpy.sqrt(Bx**2 + By**2 + Bz**2)

    Bmin = numpy.min(modB)

    Bmax = numpy.max(modB)

    modB = modB/Bmax

    x = zeros(16)

    y = zeros(16)

    z = arange(16)



    # if showplot:
    #     pass
    # else:
    #     #mlab.options.offscreen=True
    #     pass


    box = 0.5

    fig = quiver3d(x,y,z,Bx[:,0],By[:,0],Bz[:,0],mode='arrow',scale_mode='scalar',resolution=20,color=(0.1,1.0,0.1),vmin=Bmin,vmax=Bmax,scalars=modB[:,0],scale_factor=0.9,extent=(-box,box,-box,box,-9,9))

    view(azimuth=0,elevation=90,distance=30)

    source = fig.mlab_source

    frames = range(int(round(n/nskip)))

    for j in frames:

        i = j*nskip

        source.set(u=Bx[:,i],v=By[:,i],w=Bz[:,i],scalars=modB[:,i])

        savefile = 'anim'+str(j)+'.png'

        #print 'Creating image '+str(j)
        #savefig(savefile)
    clf()


if __name__ == '__main__':
    vector_movie_v0('dt2009.10.29_hr21.54',tmin=100,tmax=150,nskip=4)
