import numpy as np
from aipy import dsp
from pylab import *

c = 0.3 #m/ns
twopi = 2.*np.pi

################
# UTILITY BELT #
################

def nu2lam2(nu): return (c/nu)**2
def lam22nu(lam2): return (c/np.sqrt(lam2))
def nu_lam_swap(x): return c/x

def better_guess_l2(nu):
    dnu = np.mean(np.diff(nu))
    return (c/nu)**2 * (1. + (3.*dnu/(4.*nu)))

####################
# GENERATE SPECTRA #
####################

def gen_rm_spec(nu,rm,amp=1.):
    return amp*np.exp(2.j*rm*nu2lam2(nu))

def dlam2(nu):
    dnu = nu[1]-nu[0]
    return -2. * (c/nu)**2 * (dnu/nu)

def gen_rm_samples(nu):
    N = float(len(nu))
    j = np.arange(N)
    j -= N/2.
    return np.pi * (j / N) / np.max(np.abs(np.diff(better_guess_l2(nu))))

def Lam2Measure(nu):
    return np.abs(dlam2(nu))

#################
# REBIN AND FFT #
#################

def rebin_nu2lam2(nu,f_nu,bin=100):
    l2 = nu2lam2(nu)
    hist1,bins = np.histogram(l2, bins=bin, weights=f_nu)
    hist2,bins = np.histogram(l2, bins=bin)
    l2 = .5 * (bins[1:] + bins[:-1])
    return l2, hist1 / np.where(hist2 == 0, 1., hist2)

#######
# DFT #
#######

def RMTmat(nu,window='hamming'):
    N = len(nu)
    W = np.indices(np.array((N,N)))
    RMs = gen_rm_samples(nu)
    L2 = better_guess_l2(nu)
    wgt = dsp.gen_window(N,window)
    W = np.exp(-2.j*RMs[W[1]]*L2[W[0]]) * wgt[W[0]] * Lam2Measure(nu[W[0]]) 
    return RMs,W.T

def RMT(spec,W):
    return np.dot(W,spec)
    

########
# iDFT #
########

def iRMTmat(nu,window='hamming'):
    N = len(nu)
    W = np.indices(np.array((N,N)))
    RMs = gen_rm_samples(nu)
    L2 = better_guess_l2(nu)
    wgt = dsp.gen_window(N,window)
    W = np.exp(2.j*RMs[W[0]]*L2[W[1]]) / wgt[W[1]]
    W *= (RMs[-1]-RMs[0])/(np.pi*N)
    return W.T

def iRMT(spec,W):
   return np.dot(W,spec)















