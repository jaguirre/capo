# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 12:55:35 2016

@author: jaguirre
"""

import numpy as np
import pylab as plt
from astropy import constants as c
from astropy import units as u
import FarIRTools as fir
import aipy as a
import capo

lambda0 = c.c/(1420.*u.MHz)

def tau_to_k(tau,z):
    Y = fir.drdnu(z,lambda0)
    return (tau/Y).to(1./u.Mpc)
    
def cable_delay(cable_len,n=1.5):
    # Cable delay for echoes is 2 * cable_len / (c/n)
    return (2 * cable_len / (c.c/n)).to(u.ns)

bl_length = np.linspace(14,350,num=1000)*u.m
frequency = np.linspace(50,250,num=1000)*u.MHz

delays = (bl_length/c.c).to(u.ns)
redshifts = 1420.*u.MHz/frequency - 1.
