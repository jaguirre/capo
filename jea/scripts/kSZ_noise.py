#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 14:57:50 2017

@author: jaguirre
"""

import numpy as np
from astropy import units as u

Tsys = 400.*u.K
dnu = 8.*u.MHz
t_transit = 15.*u.min
N_days = 100*u.day
t_coherent = t_transit/u.day*N_days # transit time of t_transit, over N_days
N_bl = 4200. # Number of 14.6 m baselines
dl = 20. # Based on C_l of beam * fringe for b = 14.6 m, lambda = 2 m
N_incoherent = 8.*u.hour/t_transit # 8 useful hours in t_transit chunks

Zaldarriaga57 = (np.power(2*np.pi*Tsys,2))/(dnu * t_coherent * np.power(dl,2))
# Assuming Zaldarriaga is for single baseline
Zaldarriaga57 /= (N_bl*np.sqrt(N_incoherent))

print Zaldarriaga57.to(u.uK**2)