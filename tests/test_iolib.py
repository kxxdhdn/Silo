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
from astylo.iolib import (fclean, read_fits, write_fits,
                         read_hdf5, write_hdf5, read_ascii,
                         read_csv, write_csv, )

print('\n TEST read_fits ')
print('----------------')
get_fdata = read_fits(datdir+'M83', datdir+'M83_unc')
print('header \n', get_fdata.header)
print('header of TAB \n', get_fdata.header_w)
print('data (shape) \n', get_fdata.data.shape)
print('wave (shape) \n', get_fdata.wave.shape)
print('uncertainty (shape) \n', get_fdata.unc.shape)

print('\n TEST write_fits ')
print('----------------')
write_fits(outdir+'M83_copy',
           get_fdata.header, get_fdata.data, get_fdata.wave)
write_fits(outdir+'M83_irs_data_format_copy',
           get_fdata.header, get_fdata.data, get_fdata.wave, wmod=1)

print('\n TEST write_hdf5 ')
print('-----------------')
label = [['alpha', 'beta', 'gamma'], [1,2,3]]
data1 = np.arange(5, 23, 2).reshape((3,3))
data2 = [1,2,3]
write_hdf5(outdir+'test_iolib', 'Label', label)
write_hdf5(outdir+'test_iolib', 'Data1', data1, append=True)
write_hdf5(outdir+'test_iolib', 'Data2', data2, append=True)
print('See ./out [Done]')

print('\n TEST read_hdf5 ')
print('----------------')
get_hdata = read_hdf5(outdir+'test_iolib', 'Label', 'Data1', 'Data2')
print('Label \n', get_hdata[0])
print('Data1 \n', get_hdata[1])
print('Data2 \n', get_hdata[2])

print('\n TEST read_ascii ')
print('-----------------')
get_adata = read_ascii(datdir+'IRAC1')
print('col1 - Wave (microns) \n', get_adata[:,0])
print('col2 - Spectral Response (electrons/photon) \n', get_adata[:,1])

print('\n TEST write_csv ')
print('----------------')
write_csv(outdir+'test_iolib', ['Wave', 'Spectral Response'], get_adata)
print('See ./out [Done]')

print('\n TEST read_csv ')
print('---------------')
get_cdata = read_csv(outdir+'test_iolib', 'Spectral Response', 'Wave')
print('Spectral Response \n', get_cdata[0]) # changed list order here
print('Wave \n', get_cdata[1])

m = input('Test fclean ? (y/n)')
if m=='y':
    print('\n TEST fclean ')
    print('-------------')
    info = './out is empty! [Done]'
    fclean('out', info)
