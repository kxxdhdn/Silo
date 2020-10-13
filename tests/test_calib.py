#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, logging
# logging.disable(sys.maxsize)

## Set dir
testdir = os.path.dirname(os.path.abspath(__file__))
datdir = testdir+'/dat/'
outdir = testdir+'/out/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

## Local
sys.path.insert(0, testdir+'/..') ## astylo path
from astylo.calib import (intercalib, photometry_profile, 
)

print('\n TEST intercalib ')
print('-----------------')
print('* via FITS file *')
phots = 'IRAC1', 'IRAC4'
c1 = intercalib(datdir+'M82')
sp1 = c1.synthetic_photometry(phots)
print('Fnu_filt = ', sp1.Fnu_filt.shape)
print('wcen = ', sp1.wcen)
print('smat = ', sp1.smat)

print('* via data array *')
with fits.open(datdir+'M82.fits') as hdul:
    spec = hdul[0].data[:,34,40]
    wave = hdul[1].data
c2 = intercalib()
sp2 = c2.synthetic_photometry('IRAC3', wave, spec, verbose=False)
print('Fnu_filt = ', sp2.Fnu_filt)
print('wcen = ', sp2.wcen)
print('smat = ', sp2.smat)

print('\n TEST specorrect ')
print('-----------------')
a1 = 1.5 # slope
b1 = -200. # offset
a2 = .5
b2 = 100.
wlim11=(None,5.2)
wlim12=(5.2,14.5)
wlim2=(5.2,7.6)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8,9))
plt.subplots_adjust(left=.1, bottom=.05, \
                    right=.99, top=.95, wspace=0., hspace=.2)
ax1, ax2 = axes

print('* via FITS file *')
new_spec11 = c1.specorrect(a1, b1, wlim=wlim11)
new_spec12 = c1.specorrect(a2, b2, wlim=wlim12)
ax1.plot(c1.wvl, c1.im[:,34,40], c='k', label='y')
ax1.plot(c1.wvl, new_spec11[:,34,40], c='y',
         label='{}*y+{}, x=({},{})'.format(a1,b1,wlim11[0],wlim11[1]))
ax1.plot(c1.wvl, new_spec12[:,34,40], c='g',
         label='{}*y+{}, x=({},{})'.format(a2,b2,wlim12[0],wlim12[1]))
ax1.legend(loc='upper left')

print('* via data array *')
new_spec2 = c2.specorrect(a1, b1, w_spec=wave, Fnu_spec=spec, wlim=wlim2)
ax2.plot(wave, spec, c='k', label='y')
ax2.plot(wave, new_spec2, c='y',
         label='{}*y+{}, x=({},{})'.format(a1,b1,wlim2[0],wlim2[1]))
ax2.legend(loc='upper left')

fig.savefig(outdir+'calib_specorrect.png')
print('See ./out/calib_specorrect.png [Done]')

print('\n TEST photometry_profile ')
print('-------------------------')
phot_lib = photometry_profile(None,
                              'IRAC1', 'IRAC2', 'IRAC3', 'IRAC4', 'MIPS1',
                              'WISE1', 'WISE2', 'WISE3', 'WISE4',)
phot_lib.save(outdir+'calib_phot')
print('See ./out/calib_phot.png [Done]')

# plt.show()
