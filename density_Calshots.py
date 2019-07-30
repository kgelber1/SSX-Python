# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 10:46:03 2018

@author: Manjit

This function contains the density calibration shots for each day
"""

def calshot(day):
    if day == '101718':
        return '101718r2'
    if day =='101618':
        return '101618r3'
    if day == '100218':
        return '100218r5'
    if day == '081518':
	    return '081518r2'
    if day == '083018':
	    return '083018r4'
    if day == '082118':
	    return '082118r5'
    if day == '081418':
        return '081418r4'
    if day == '081717':
        return '081717r2'
    if day == '101917':
        return '101917r5'
    if day == '070918':
        return '070918r5'
    if day == '071218':
        return '071218r1'
    if day == '101917':
        return '101917r5'
    if day == '030118':
        return '030118r1'
    if day == '101317':
        return '101317r5'
    if day == '101117':
        return '101117r4'
    if day == '091818':
        return '091818r4'
    if day == '061219':
        return '061219r19'
    if day == '062619':
        return '062619r3'
    if day == '071019':
        return '071019r2'
    if day == '071819':
        return '071819r2'
    if day == '072419':
        return '072419r3'
    if day == '073019':
        return '073019r5'
    else:
        print('Day/shot not found')
