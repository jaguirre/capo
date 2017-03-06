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

def VisualizeAlm(alm,figno=1,max_l=None,annot=''):
    """ Visualize a healpy a_lm vector """
    lmax = hp.Alm.getlmax(f_lm.size)
    l,m = hp.Alm.getlm(lmax)
    mag = np.zeros([lmax+1,lmax+1])
    phs = np.zeros([lmax+1,lmax+1])
    a_lm = np.zeros([lmax+1,lmax+1],dtype='complex128')
    mag[m,l] = np.abs(alm)
    phs[m,l] = np.angle(alm)
    a_lm[m,l] = alm
    cl = hp.alm2cl(alm)
    # Decide the range of l to plot
    if max_l != None:
        max_l = (max_l if (max_l <= lmax) else lmax)
    else:
        max_l = lmax 
    print max_l
    plt.figure(figno,figsize=(4,4))
    plt.clf()
    plt.subplot(211)
    plt.imshow(mag[0:max_l,0:max_l],interpolation='nearest',origin='lower')
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$m$')
    plt.colorbar()
    plt.title(annot+'Magnitude')
    plt.subplot(212)
    plt.imshow(phs[0:max_l,0:max_l],interpolation='nearest',origin='lower')
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$m$')
    plt.colorbar()
    plt.title(annot+'Phase')
    # plt.subplot(313)
    #plt.semilogy(cl[0:max_l])
    return {'mag':mag,'phs':phs,'cl':cl,'a_lm':a_lm}
    
#%%

# re-inventing the wheel
def calc_delay(n,df):
    # Convert from MHz to Hz and tau back to ns
    tau = np.fft.fftshift(np.fft.fftfreq(n,d=df*1e6)*1e9)
    return tau

def hp2full_a_lm(a_lm,lmax):
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
lmax = 3*nside-1
ell = np.arange(0,lmax+1)
l,m = hp.Alm.getlm(lmax)
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
theta,phi = hp.pix2ang(nside,ipix)
below_horizon = np.where(theta > np.pi/2.)[0]

nu = np.linspace(100,200,num=203)*u.MHz
dnu = np.median(np.diff(nu))
tau = calc_delay(len(nu),df=dnu)

#---

runname='PAPER-like'
## Instrument parameters
lmbda0 = 2.
D = 2.
fwhm = 1.22*lmbda0/D
sigma = fwhm/2.35
b0 = 15.
tau0 = (b0*u.m/c.c).to(u.ns)
l_pk = 2.*np.pi*D/lmbda0

## Construct beam
sidelobe_level = 3e-4
beam = np.exp(-np.power(theta,2)/(2.*np.power(sigma,2)))
beam += sidelobe_level
#beam_nu += sidelobe_level
#beam = np.power(np.sin(theta),2)*np.power(np.sin(phi),16)
beam /= beam.max()
beam_nu = np.outer(beam,np.ones_like(nu))

## Construct fringe
bvec = np.array([0,1,0])*b0
b = np.outer(bvec,np.ones(npix))*u.m
s = np.array(hp.pix2vec(nside,ipix))
bdots = np.sum(b*s,axis=0)
fringe = np.exp(-2.*np.pi*1j*np.outer(bdots.value,nu.to(u.Hz).value)/c.c.value)
if True:
    fringe[below_horizon] = 0.
    beam[below_horizon] = 0.

## V_G(nu) = s_00(u) T(nu) 
# Average compensates for any nside changes, normalized to beam solid angle (?)
omega_nu = 4.*np.pi*beam_nu.sum(axis=0)/npix
T_nu = np.average(beam_nu*fringe,axis=0)
window = np.hanning(len(T_nu))
dtrans_T_nu = np.fft.fftshift(np.fft.fft(window*T_nu))
Ptau_T_nu = np.abs(dtrans_T_nu)

# Need to calculate this for all frequencies
f_lm = hp.map2alm(fringe[:,100])
C_f = hp.alm2cl(f_lm)
a_lm = hp.map2alm(beam_nu[:,100])
C_a = hp.alm2cl(a_lm)

C_fa = hp.alm2cl(f_lm*a_lm)



plt.figure(1,figsize=(12,7))
plt.clf()
plt.subplot(231)
plt.plot(nu,T_nu.real,'b',label=r'$Re[\Xi(\nu)]$')
plt.plot(nu,T_nu.real-T_nu.mean(),'b--',label=r'$Re[\Xi(\nu)-\bar{\Xi}(\nu)]$')
plt.plot(nu,window*T_nu.real,'g',label=r'$Re[\Xi(\nu)] W(\nu)$')
plt.plot(nu,T_nu.imag,'r',label=r'$Im[\Xi(\nu)]$')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Amplitude (arb. units)')
plt.legend(loc='upper right')
plt.title('Global signal transfer function')

plt.subplot(232)
plt.plot(nu,T_nu.real/omega_nu,label=r'$Re[\Xi(\nu)]/\Omega_\nu$')
plt.legend(loc='upper right')
plt.title('Normalized global signal transfer function')
plt.xlabel('Frequency [MHz]')

plt.subplot(233)
plt.semilogy(tau,Ptau_T_nu)

plt.semilogy(-tau0.value*np.array([1.,1.]),[Ptau_T_nu.min(),Ptau_T_nu.max()],'r--')
plt.semilogy(tau0.value*np.array([1.,1.]),[Ptau_T_nu.min(),Ptau_T_nu.max()],'r--')
plt.title('Delay transform of transfer function')
plt.xlim([-3*tau0.value,3*tau0.value])
plt.ylim([Ptau_T_nu.max()*1e-6,Ptau_T_nu.max()*2.])
plt.xlabel('Delay [ns]')

plt.subplot(234)
plt.plot(C_f/C_f.max(),label=r'$C_{\ell,fringe}$')
plt.plot(C_a/C_a.max(),label=r'$C_{\ell,beam}$')
plt.xlim([0,l_pk*1.5])
plt.xlabel(r'Multipole $\ell$')
plt.ylabel(r'$C_\ell/C_0$')
plt.legend()

plt.subplot(235)
plt.plot(ell,C_fa/C_fa.max(),label=r'$\frac{1}{2 \ell+1} \sum_m |f_{lm} a_{lm}|^2$')
plt.xlim([0,100])
plt.xlabel(r'Multipole $\ell$')
plt.ylabel(r'$C_\ell/C_0$')
plt.legend(loc='upper right')

hp.orthview((beam_nu*fringe)[:,100].real,rot=[0,90],half_sky=True,sub=236,title='')
hp.graticule()
plt.savefig(runname+'_spectrum_delay.png')

#%%

valm = VisualizeAlm(f_lm*a_lm,figno=2,max_l = int(l_pk*1.2),annot=r'$a_{lm} f_{lm}$ ')
plt.savefig(runname+'_a_lm_visualization.png')

plt.figure(3)
plt.plot(valm['a_lm'][0,:].real,'bo-',label=r'Re[$f_{\ell 0} a_{\ell 0}]')
plt.plot(valm['a_lm'][0,:].imag,'g',label=r'Re[$f_{\ell 0} a_{\ell 0}]')
plt.xlim([0,int(l_pk*1.2)])
plt.xlabel(r'Multipole $\ell$')
plt.ylabel(r'$f_{\ell 0} a_{\ell 0}(\nu=150 \, \rm{MHz} )$')
plt.subplots_adjust(left=0.2)
plt.savefig(runname+'_xi_of_ell.png')
    
    
    
    
    
    
    
    