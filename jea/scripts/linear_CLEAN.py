# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 14:21:37 2016

@author: jaguirre
"""

import aipy as a
import sys
import capo
from glob import glob
import pylab as plt
import numpy as np
import time
import h5py
from astropy import units as u
from astropy import constants as c
import pandas as pd
#%%
def stime(message=''):
    print message
    t0 = time.time()
    return t0

def etime(t0):
    dt = time.time() - t0
    print 'Took',dt,'seconds'
    return dt

#%%
def flag_avg(data,flags,axis=0):
    sum_flags = flags.sum(axis=axis)
    sum_data = (data*flags).sum(axis=axis)
    avg = sum_data/sum_flags
    avg[sum_flags == 0] = 0
    avg_flag = np.zeros_like(avg.real)
    avg_flag[np.isfinite(avg)] = 1
    return (avg,avg_flag)
#%%

#%% From Saul, untested
def dly_xform(d,flags,tol = 1e-3):
    """ d is shape (Ntime, Nfreq) """
    w = a.dsp.gen_window(d.shape[-1], window='blackman-harris')
    _dw = np.fft.ifft(d*w)
    _ker= np.fft.ifft(flags*w)
    gain = a.img.beam_gain(_ker)
    for t in range(_dw.shape[0]):
        _dw[t,:],info = a.deconv.clean(_dw[t,:], _ker[t,:], tol=tol)
        _dw[t,:] += info['res']/gain
    # Dangerous, potentially
    _dw = np.fft.fftshift(_dw)
    return _dw

caldir='/Users/jaguirre/PyModules/capo/jea/cals/'
sys.path.append(caldir)
calfile = 'psa6622_v003'
hdf5dir = '/Users/jaguirre/Documents/PAPER/GlobalSignal/'
#%%
#filename = 'zen.2456680.xx.uvcRREO.hdf5'
filename='HERA.hdf5'
fileptr = h5py.File(hdf5dir+filename,'r')
# Gotta be a better way
def copy_dataset(dataset):
    datadict = {}
    keys = dataset.keys()
    for key in keys:
        datadict[key] = dataset[key].value
    return datadict
vis = copy_dataset(fileptr['waterfalls'])
flags = copy_dataset(fileptr['iflags'])
#%%
iflags = {}
bls = vis.keys()
for bl in bls:
    iflags[bl] = np.ones_like(vis[bl])
    iflags[bl][flags[bl]==True] = 0
#%%
# Argh. How to know this in general?
#bl = '072_096'
bl = '020_089'
b_len = 14.6*u.m
wf = fileptr['waterfalls'][bl].value
# How did this get messed up?  What is 1 is supposed to be 0.
ifl = np.abs(1.-fileptr['iflags'][bl].value) 
wfnan = wf/ifl
ntime,nfreq = wf.shape
#%%
n = wf.shape[1]
f,df = np.linspace(0.1+0.1/(2.*n),0.2,num=n,retstep=True) # GHz
tau = np.fft.fftfreq(n,d=df)


#%%

# Build a matrix
def tau_select(tau,b_len,tau_pad=160.):
    """ Select taus to fit based on baseline length """
    tau_max = (b_len/c.c).to(u.ns).value
    #tau_pad = 150. # ns.  That's CRAZY!  Why do I need so many modes?
    imax = np.where((tau>0)*(tau<=tau_max+tau_pad))[0].max()
    taus = tau[1:imax+1]
    return (taus,imax)

def buildA(taus,f):
    nfreq = len(f)
    #ntaus = len(taus)
    #A = np.zeros([nfreq,2*ntaus+1])
    #A[:,0] = np.ones(nfreq) # Constant term
    arg = 2.*np.pi*np.outer(f,taus)
    A = np.concatenate((np.ones([nfreq,1]),np.cos(arg),np.sin(arg)),axis=1)
    return A

#taus,imax = tau_select(tau,b_len,tau_pad=160) #20
#print 'imax',imax
#    
#A = buildA(taus,f_tofit)
#A_all = buildA(taus,f)
#
##
#x,residuals,rank,singular = np.linalg.lstsq(A,data_tofit.real)
#model = np.dot(A_all,x)

#%%
#if True:
def linCLEAN(wf,ifl,b_len,tau_pad=160.):
    #resid = np.zeros_like(wf)
    t0 = stime('One baseline')
    model = np.zeros_like(wf)
    ntime,nfreq = wf.shape
    f,df = np.linspace(0.1+0.1/(2.*nfreq),0.2,num=nfreq,retstep=True) # GHz
    tau = np.fft.fftfreq(nfreq,d=df)
    taus,imax = tau_select(tau,b_len,tau_pad)
    for i in np.arange(ntime):
        #print i
        whf=np.where(ifl[i,:])[0]
        f_tofit = f[whf]
        A = buildA(taus,f_tofit)
        xreal,residuals,rank,singular = np.linalg.lstsq(A,wf[i,whf].real)
        ximag,residuals,rank,singular = np.linalg.lstsq(A,wf[i,whf].imag)
        model[i,whf] = np.dot(A,xreal) + 1j*np.dot(A,ximag)
    resid = wf - model
    etime(t0)
    return (model,resid)

#%%
model,resid = linCLEAN(wf,ifl,b_len,tau_pad=200.)

#%% Calculate a running RMS of the residual
resid_nan = resid/ifl
resid_rms = np.zeros_like(resid.real)
resid_rms_poly = np.zeros_like(resid.real)
for i in np.arange(ntime):
    ser = pd.Series(resid_nan[i,:].real)
    resid_rms[i,:] = ser.rolling(window=10).std()
    wh = np.isfinite(resid_rms[i,:])
    poly_params = np.polyfit(f[wh],resid_rms[i,wh],5)
    resid_rms_poly[i,:] = (np.poly1d(poly_params))(f)

#%% Some stuff with globsig investigation
avg,avg_flag = flag_avg(wf,ifl)
#avg,avg_flag = flag_avg(model,np.ones_like(model.real))
avg = np.reshape(avg,[1,len(avg)])
avg_flag = np.reshape(avg_flag,[1,len(avg_flag)])
d_avg = dly_xform(avg,avg_flag)

plt.figure(2)
plt.clf()
d_spec = np.abs(d_avg[0])
plt.plot(np.fft.fftshift(tau),np.roll(d_spec,-4))
tau_hor = ((b_len/c.c).to(u.ns)).value
plt.axvline(x=0,color='r')
plt.axvline(x=-tau_hor,color='r',linestyle='--')
plt.axvline(x=tau_hor,color='r',linestyle='--')
plt.axvline(x=-tau_hor*2,color='g',linestyle='--')
plt.axvline(x=tau_hor*2,color='g',linestyle='--')
plt.xlim([-tau_hor*5,tau_hor*5])
#%%
plt.figure(1)
plt.clf()
plt.subplot(131)
plt.imshow((wf/ifl).real,aspect='auto',vmin=-0.3,vmax=0.3)
plt.colorbar()
plt.subplot(132)
plt.imshow((model/ifl).real,aspect='auto',vmin=-0.3,vmax=0.3)
plt.colorbar()
plt.subplot(133)
plt.imshow((resid/ifl).real,aspect='auto',vmax=0.01,vmin=-0.01)
plt.colorbar()


#%% That residual ... go there and back again
#resid_arp = (data - model)*flags


#%%
if False:
    plt.figure(1)
    plt.clf()
    plt.plot(f,datanan.real,label='Data')
    #plt.plot(f_tofit,data_tofit,'.',label='Flagged Data')
    plt.plot(f,model,label='Model')
    plt.plot(f,(data-model)/flags,label='Residual')
    plt.legend()

