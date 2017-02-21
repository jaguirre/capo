#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:15:27 2017

@author: jaguirre
"""

def toymodel(N,sigma_phi,lmbda):
    sigma = sigma_phi * np.power(lmbda,2)
    return 1./N + (1-1./N)*np.exp(-np.power(sigma,2))    