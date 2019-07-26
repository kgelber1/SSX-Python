import numpy as np
import scipy as sp
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
import os


#And from the SSX files:
import basic_hiresmag as hdr
import mjtools as mj
import ssx_data_read_basic as sdr

"""
#######################################
USE get_data() as the general function.
#######################################
(Lookup calibration)

also paired with that is get_Bxy_vecRendering_format_lookup() for the
animation functions



TODO: get the matrix thing working.
if you have the matrix thing working, use getMagData()


--KG 07/25/19
"""

def get_just_time(shot):
    magdata = hdr.getMagData(shot)
    return magdata.time

def get_data(shot, num_probes = 16):
    """Gets the data from the 4x4 probe setup at ALL TIMES

    Returns data, the calibrated, integrated magentic field data
    using lookup calibration.

    Data has the dimensions of [probes][direction]
        - probes is the number of probes, from 0 to 15 typically
        - direction is the direction of the magnetic field
            - 0 is x-dir
            - 1 is y-dir
            - 2 is z-dir

    ie Bx of probe 5 (index starting from 0) can be found:

        Bx = data[5][0]

    ASSUMES THE CALIBRATION FILE IS IN magcalib
    (feel free to change that tho)"""

    data = [[0] *3 for i in range(num_probes)]
    magdata = hdr.getMagData(shot)

    calibFile = 'calib-050119-4x4_lookup_final.txt'
    calibFile = os.path.join(os.getcwd() + '\\magcalib\\' , calibFile)#change where caliv file is here
    calibFactor = np.loadtxt(calibFile)

    for probe in range(num_probes):
        Bx_all=cumtrapz(magdata.unCalibData[0,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[0,probe,:], num_probes-1), magdata.time)
        Bx = Bx_all * calibFactor[probe][0]

        By_all=cumtrapz(magdata.unCalibData[1,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[1,probe,:], num_probes-1), magdata.time)
        By = By_all * calibFactor[probe][1]

        Bz_all=cumtrapz(magdata.unCalibData[2,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[2,probe,:], num_probes-1), magdata.time)
        Bz = Bz_all * calibFactor[probe][2]

        data[probe] = [Bx, By, Bz]


    return magdata.time, data


def get_Bxy_vecRendering_format_lookup(shot, fix = False, num_probes = 16):
    """Basically get_data() but returns data in a form that is compatible
       with the b-field animations found in anim_bfield_merging files


       Returns data, the calibrated, integrated magentic field data
       using lookup calibration.

    --KG 07/25/19
    """

    data = [[0] *3 for i in range(num_probes)]
    # B25 = [[0] *num_probes for i in range(3)]

    magdata = hdr.getquikData(shot)
    timeB = magdata.time

    calibFile = 'calib-050119-4x4_lookup_final.txt'
    calibFile = os.path.join(os.getcwd() + "\\magcalib\\", calibFile)#change where caliv file is here
    calibFactor = np.loadtxt(calibFile)


    B25 = np.zeros([3,num_probes, len(timeB)-1])


    for probe in range(num_probes):
        Bx_all=cumtrapz(magdata.unCalibData[0,probe,:]-mj.get_gaussian_moving_avg(timeB,
        magdata.unCalibData[0,probe,:], num_probes-1), timeB)
        Bx = Bx_all * calibFactor[probe][0]

        By_all=cumtrapz(magdata.unCalibData[1,probe,:]-mj.get_gaussian_moving_avg(timeB,
        magdata.unCalibData[1,probe,:], num_probes-1), timeB)
        By = By_all * calibFactor[probe][1]

        Bz_all=cumtrapz(magdata.unCalibData[2,probe,:]-mj.get_gaussian_moving_avg(timeB,
        magdata.unCalibData[2,probe,:], num_probes-1), timeB)
        Bz = Bz_all * calibFactor[probe][2]

        B25[0, probe, :] = Bx
        B25[1, probe, :] = By
        B25[2, probe, :] = Bz


    if fix:
        B25[2, 9, :] = B25[2, 9, :]* 100
        B25[2, 11, :] = B25[2, 11, :]* 100
        B25[1,12, :] = B25[1,12, :]* 100

    Bmod25 = np.zeros([num_probes, len(timeB)-1])
    for j in range(num_probes):
        Bmod25[j,:] = np.sqrt(B25[0,j,:]**2+B25[1,j,:]**2+B25[2,j,:]**2)

    return (timeB,B25)


def get_Bxy_vecRendering_format(shot, fix = False, num_probes = 16):
    """
    Returns data, the calibrated, integrated magentic field data
    using matrix calibration.

    UNFINISHED
    becuase of the matrix thing. But this function would return getMagData()
    data but in a version compatible with the animations
    in anim_bfield_merging.py

    --KG 07/25/19


    """

    data = [[0] *3 for i in range(num_probes)]
    # B25 = [[0] *num_probes for i in range(3)]

    magdata = hdr.getquikData(shot)
    timeB = magdata.time

    calibFile = 'calib-050119-4x4_matrix.txt'
    calibFile = os.path.join(os.getcwd() +  '\\magcalib\\', calibFile)#change where caliv file is here


    print("\n\nHERE")
    calibFactor = np.loadtxt(calibFile)

    # calibFactor[14] = [0,0,0]
    # calibFactor[6] = [0,0,0]
    # calibFactor[7] = [0,0,0]

    B25 = np.zeros([3,num_probes, len(timeB)-1])
    B_vec = np.zeros([3,num_probes, len(timeB)-1])
    win_size = 25

    for probe in range(num_probes):
        for i in range(3):
            B_vec[i, probe, :]=cumtrapz(magdata.unCalibData[i,probe,:]-mj.get_gaussian_moving_avg(timeB,
            magdata.unCalibData[i,probe,:], win_size), timeB)


    for t in range(len(timeB)-1):
        for probe in range(num_probes):
            C = np.reshape(calibFactor[probe], (3,3))
            #put into 3 x1 matric
            V = [[B_vec[0,probe,t]],[ B_vec[1,probe,t]],[ B_vec[2,probe,t]]]

            B = np.dot(C, V)
            Bx, By, Bz = B

            B25[0, probe, t] = Bx
            B25[1, probe, t] = By
            B25[2, probe, t] = Bz

    Bmod25 = np.zeros([num_probes, len(timeB)-1])

    for j in range(num_probes):
        Bmod25[j,:] = np.sqrt(B25[0,j,:]**2+B25[1,j,:]**2+B25[2,j,:]**2)

    return (timeB,B25)

def getMagData(shot, t0 = -1, tf = -1, num_probes = 16):
    """The main function that should be used

       Gets the data from the 4x4 probe setup at ALL TIMES

       Returns Bx,By and Bz, the calibrated, integrated magentic field data
       using a matrix calibration

       Data has the dimensions of [probes][time]
            - probes is the number of probes, from 0 to 15 typically

       assumes the calibration file is stored in magcalib folder


        --KG 07/10/19"""

    magdata = hdr.getquikData(shot)
    timeB = magdata.time

    if t0 != -1 and tf !=-1:
        #trim to desired time
        t_index, t = my.fix_bounds(np.array((t0, tf)), timeB)
        magdata = magdata[:,:,t0:tf]

    calibFile = 'calib-050119-4x4_matrix.txt'
    calibFile = os.path.join(os.getcwd() + '\\magcalib\\', calibFile)#change where caliv file is here
    calibFactor = np.loadtxt(calibFile)

    B = np.zeros([3,num_probes, len(timeB)-1])

    win_size = 25

    for probe in range(num_probes):
        for i in range(3):
            B_vec[i, probe, :]=cumtrapz(magdata.unCalibData[i,probe,:]-mj.get_gaussian_moving_avg(timeB,
            magdata.unCalibData[i,probe,:], win_size), timeB)


    for t in range(len(timeB)-1):
        for probe in range(num_probes):
            C = np.reshape(calibFactor[probe], (3,3))
            #put into 3 x1 matric
            V = [[B_vec[0,probe,t]],[ B_vec[1,probe,t]],[ B_vec[2,probe,t]]]

            B_data = np.dot(C, V)
            B[0, probe,t], B[1, probe,t], B[2, probe,t] = B_data

    return B


def get_data_slice(shot, t, num_probes = 16):
    """Gets the data from the 4x4 probe setup AT A SPECIFIC SLICE OF TIME


       Returns data, the calibrated, integrated magentic field data

       Data has the dimensions of [probes][direction][time]
            - probes is the number of probes, from 0 to 15 typically
            - direction is the direction of the magnetic field
                - 0 is x-dir
                - 1 is y-dir
                - 2 is z-dir

       ie Bx of probe 5 (index starting from 0) can be found:

            Bx = data[5][0]

        Very time inefficent, so don't plan to use this often

        --KG 07/25/19
    """

    data = [[0] *3 for i in range(num_probes)]
    magdata = hdr.getMagData(shot)

    calibFile = 'calib-050119-4x4_lookup_final.txt'
    calibFile = os.path.join(os.getcwd() + '\\magcalib\\', calibFile)#change where caliv file is here
    calibFactor = np.loadtxt(calibFile)

    for probe in range(num_probes):
        Bx_all=cumtrapz(magdata.unCalibData[0,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[0,probe,:], num_probes-1), magdata.time)
        Bx = Bx_all[int(t)] * calibFactor[probe][0]

        By_all=cumtrapz(magdata.unCalibData[1,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[1,probe,:], num_probes-1), magdata.time)
        By = By_all[int(t)] * calibFactor[probe][1]

        Bz_all=cumtrapz(magdata.unCalibData[2,probe,:]-mj.get_gaussian_moving_avg(magdata.time,
        magdata.unCalibData[2,probe,:], num_probes-1), magdata.time)
        Bz = Bz_all[int(t)] * calibFactor[probe][2]

        data[probe] = [Bx, By, Bz]


    return data
