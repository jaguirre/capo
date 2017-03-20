# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 14:21:37 2016

@author: jaguirre
"""

import aipy
import sys
import capo
from glob import glob
import pylab as plt
import numpy as np
import time
import h5py
#%%
# Useful functions
def stime(message=''):
    print message
    t0 = time.time()
    return t0

def etime(t0):
    dt = time.time() - t0
    print 'Took',dt,'seconds'
    return dt

def tup2str(tup):
    """ Take a typical two element tuple describing a baseline and convert 
    to the string name used in the hdf5 file """
    if tup[0] < tup[1]:
        i = tup[0]
        j = tup[1]
    else:
        i = tup[1]
        j = tup[0]
    blstr = (str(i).zfill(3)+'_'+str(j).zfill(3))
    return blstr
    
def flag_avg(data,flags):
    """ iflags here are floats """
    avg = (data*flags).sum(axis=0)/flags.sum(axis=0)
    return avg
    
#%%
# Paths    
caldir='/Users/jaguirre/PyModules/capo/jea/cals/'
sys.path.append(caldir)
calfile = 'psa6622_v003'
hdf5dir = '/Users/jaguirre/Documents/PAPER/GlobalSignal/'
#%%
# Setup the information describing the PAPER array
ants = np.arange(0,128)
nant = len(ants)
nbl = nant*(nant+1)/2.
polstr = 'xx'
nfreqs = 203

freqs = np.linspace(100.,200.,num=nfreqs)

aa = aipy.cal.get_aa(calfile,np.array([.15]))
info = capo.omni.aa_to_info(aa)
reds = info.get_reds()

#%%
# Find (10,64), which is known to be a XXX m E-W baseline
myblgrp = -1
for i,group in enumerate(reds):
    for bl in group:
        if (bl == (10,64) or bl == (64,10)):
            myblgrp = i
#%%
badants = [18,5,61,72,76,77,78,79,16,26,34,38,46,50,84,99]
# Open the HDF file
f = h5py.File(hdf5dir+'zen.2456680.xx.uvcRREO.hdf5','r')


#%%
t0 = stime('Starting read')
avgspecs = np.zeros([nant,nant,nfreqs],dtype='complex128')
ibl = 1
for i in range(nant):
    print ibl,' of ',nbl
    for j in range(i,nant):
        blstr = tup2str((i,j))
        avgspecs[i,j,:]=flag_avg(f['waterfalls'][blstr].value,
        f['iflags'][blstr].value)
        ibl += 1
etime(t0)
#%%        
_g# Edit out bad antennas
reds_good = []
for bl in reds[myblgrp]:
    i = bl[0]
    j = bl[1]
    if i in badants or j in badants:
        continue
    reds_good.append((i,j))
    
#%%    
plt.figure(1)
plt.clf()
reds_really_good = []
for bl in reds_good:
    i = bl[0]
    j = bl[1]
    dev = np.nanstd(np.abs(avgspecs[1,47,:]-avgspecs[i,j,:]))
    if dev < 0.019:
        reds_really_good.append(bl)
        print bl
        plt.plot(freqs,np.abs(avgspecs[i,j,:]))

print len(reds_good)
print len(reds_really_good)

#%%
plt.figure(2)
plt.clf()
for num,bl in enumerate(reds_good):
    i = bl[0]
    j = bl[1]  
    plt.subplot(10,9,num+1)
    print bl,num+1
    plt.plot(freqs,np.abs(avgspecs[i,j,:]),label=tup2str((i,j)))
    plt.xlim([100,200])
    plt.ylim([0,0.08])
    plt.legend()

plt.show()

    


    
