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

## Local
sys.path.insert(0, testdir+'/..') ## astylo path
from astylo.calib import (intercalib, photometry_profile)

print('\n TEST intercalib ')
print('-----------------')
print('* via FITS file *')
phots = 'IRAC1', 'IRAC4'
ca = intercalib(datdir+'M82').synthetic_photometry(phots)
print('Fnu_filt = ', ca.Fnu_filt.shape)
print('wcen = ', ca.wcen)
print('smat = ', ca.smat)

print('* via data array *')
with fits.open(datdir+'M82.fits') as hdul:
    spec = hdul[0].data[:,30,43]
    wave = hdul[1].data
ca1 = intercalib().synthetic_photometry('IRAC3', wave, spec, verbose=False)
print('Fnu_filt = ', ca1.Fnu_filt)
print('wcen = ', ca1.wcen)
print('smat = ', ca1.smat)

print('\n TEST photometry_profile ')
print('-------------------------')
phot_lib = photometry_profile(None,
                              'IRAC1', 'IRAC2', 'IRAC3', 'IRAC4', 'MIPS1',
                              'WISE1', 'WISE2', 'WISE3', 'WISE4',)
phot_lib.save(outdir+'phot_profiles')
# phot_lib.show()
print('TEST photometry_profile [Done]')

