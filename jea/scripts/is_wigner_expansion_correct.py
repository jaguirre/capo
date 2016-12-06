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

def WignerFac(l1,l3,m1,m3):
    fac1 = wp.wigner3j(l1,l3,0, 0,0,0)
    fac2 = wp.wigner3j(l1,l3,0 ,m1,-m3,0)
    return fac1*fac2

def full_a_lm(a_lm,lmax):
    l,m = hp.Alm.getlm(lmax)
    assert (len(a_lm) == len(l))
    A_lm = {}
    for L in np.arange(lmax+1):
        Mk = {}
        for M in np.arange(-L,L+1):
            i = hp.Alm.getidx(lmax,L,np.abs(M))
            if M < 0:
                a = np.power(-1,M)*np.conj(a_lm[i])
            else:
                a = a_lm[i]
            Mk[M] = a
        A_lm[L] = Mk
    return A_lm

nside = 128
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
theta,phi = hp.pix2ang(nside,ipix)
lmax = 3*nside-1
l,m = hp.Alm.getlm(lmax)
n_lm = len(l)
ell = np.arange(0,lmax+1)

below_horizon = np.where(theta > np.pi/2.)[0]

nu = np.linspace(100,200,num=203)*u.MHz

# Define primary beam
#A = np.exp(-np.power(theta,2)/(2.*0.1))
#A = np.exp(-np.power(theta,2)/(2.*0.1))
A = np.ones_like(theta)
A[below_horizon] = 0.
a_lm = hp.map2alm(A,lmax=lmax)
Cl_A = hp.alm2cl(a_lm,lmax=lmax)
Cl_G = hp.gauss_beam(np.radians(90.),lmax)
A_nu = np.outer(A,np.ones_like(nu))

#apod = np.ones_like(theta)
#apod = np.exp(-np.power(theta,2)/(2.*0.5))
#disc = hp.query_disc(nside,[0,0,1],np.radians(80.))
#apod[disc] = 0.
#apod = 1-np.exp(-np.power(theta,2)/(2.*0.5))
#apod = 
#apod = np.power(np.sin(phi),2)

# Define fringe 
bmag = 30.
bvec = np.array([0,1,0])*bmag
b = np.outer(bvec,np.ones(npix))*u.m
s = np.array(hp.pix2vec(nside,ipix))
bdots = np.sum(b*s,axis=0)
F_nu = np.exp(-2.*np.pi*1j*np.outer(bdots.value,nu.to(u.Hz).value)/c.c.value)
#fringe *= np.outer(apod,np.ones_like(nu))
f_lm = hp.map2alm(F_nu[:,100],lmax=lmax)
Cl_F = hp.alm2cl(f_lm,lmax=lmax)

dOmega = 4.*np.pi/npix
V_0_map = np.sum(A_nu*F_nu*dOmega,axis=0)

W0 = np.zeros(n_lm)

# Brute force
V_0 = np.zeros_like(a_lm)
for i,(L,M) in enumerate(zip(l,m)):    
    W0[i] = WignerFac(L,L,M,M)
    #print i,L,M,W0[i],WignerFac(L,L,M,M)
    V_0[i] = a_lm[i] * f_lm[i] * W0[i] 

taunu = (2.*np.pi*bmag*u.m/c.c*nu).to(u.dimensionless_unscaled).value
plt.figure(1)
plt.clf()
# Can't get the normalization to quite work out
# Also, V_0_map pick up a mean
V_0_to_plot = V_0_map.real
V_0_to_plot -= V_0_to_plot.mean()
V_0_to_plot /= V_0_to_plot.real.max()
plt.plot(nu,V_0_to_plot,label='Simulation')
sinc = np.sin(taunu)/taunu
sinc /= sinc.max()
plt.plot(nu,sinc,label='Analytic l=0')
plt.legend()


if False:
    plt.figure(1)
    plt.clf()
    plt.plot(Cl_A/Cl_A.max())
    plt.plot(Cl_G)
    plt.xlim([0,30])    
    plt.plot(l,np.abs(V_0)/np.abs(V_0).max(),'o')
    
    
A_lm = full_a_lm(a_lm,lmax)
    
    
    
    
    
    
    
    