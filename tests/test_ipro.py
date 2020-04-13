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
from astropy.wcs import WCS

## Local
sys.path.insert(0, testdir+'/..') ## astylo path
from astylo.ipro import (improve, islice, icrop, iswarp, iconvolve,
                          sextract, wclean, interfill, hextract, hswarp)
"""
print('\n TEST improve ')
print('--------------')
imp = improve(datdir+'M83', verbose=True)
print('raw: ', imp.im[0,0,:])
imp.rand_norm(datdir+'M83_unc')
print('add sym rand: ', imp.im[0,0,:])
imp.rand_splitnorm([datdir+'M83_unc', datdir+'M83_unc'])
print('add split rand', imp.im[0,0,:])
imp.crop(outdir+'M83_imp_crop', sizval=(0.005,0.009), cenpix=(7,8))
print('See out/M83_imp_crop.fits [Done]')

print('\n TEST islice ')
print('-------------')
slc = islice(datdir+'M83', filSL=outdir+'M83_inv_sqrt',
             filUNC=datdir+'M83_unc', slicetype='inv_sqrt', postfix='_test_')
print(slc.filenames()[0])
if input('Clean slices ? (y/n) ')=='y':
    slc.clean()

print('\n TEST icrop ')
print('------------')
crp = icrop(datdir+'M83', filOUT=outdir+'M83_icrop',
            sizpix=(3,6), cenval=(204.2529675, -29.8656962),
            filUNC=[datdir+'M83_unc', datdir+'M83_unc'],
            dist='splitnorm', verbose=True)
print('See out/M83_icrop.fits [Done]')
"""
print('\n TEST iswarp ')
print('-------------')
print('Reproject M82 (towards north) to M83 orientation')
swp = iswarp(datdir+'M83',
             center=None, pixscale=5.,
             tmpdir=outdir+'swp/') # M82 center
swp.combine(datdir+'M82', uncpdf='norm',
            filOUT=outdir+'M82_rep_to_M83', tmpdir=outdir+'swp/')

print('\n TEST iconvolve ')
print('----------------')

print('\n TEST sextract ')
print('---------------')

print('\n TEST wclean ')
print('-------------')

print('\n TEST interfill ')
print('----------------')

print('\n TEST hextract ')
print('---------------')

print('\n TEST hswarp ')
print('-------------')
