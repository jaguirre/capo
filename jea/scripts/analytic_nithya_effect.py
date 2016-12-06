# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:56:50 2016

@author: jaguirre
"""

import numpy as np
import healpy as hp
from astropy import units as u
from astropy import constants as c
import pylab as plt

import wignerpy._wignerpy as wp

""" The argument goes that we're interested in what couples to 
l_sky = m_sky = 0, the constant mode of the sky """

l_beam = np.arange(0,3)
l_fringe = np.arange(0,3)

for l_b in l_beam:
    for l_f in l_fringe:
        for m_b in np.arange(-l_b,l_b+1):
            for m_f in np.arange(-l_f,l_f+1):
                print l_b, 0, l_f
                print m_b, 0, m_f
                fac1 = wp.wigner3j(l_b, 0, l_f,0,0,0)
                fac2 = wp.wigner3j(l_b, 0, l_f,m_b, 0, m_f)
                print fac1, fac2, fac1*fac2
                print

    
    
    
    
    
    
    
    
    
    
    
    
    