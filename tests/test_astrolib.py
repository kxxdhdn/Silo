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
from astylo.astrolib import (pix2sr, sr2arcsec2, rad2arcsec, hour2deg,
                             deg2hour, get_cd, get_pc, fixwcs, )

pc = np.array(((-0.928395134922, -0.371594501376), (0.371594501376, -0.928395134922)))
cdelt = np.array((-0.00141111109406, 0.00141111109406))
cd = np.array(((0.00035731, 0.00036933), (0.00036933, -0.00035731)))

print('\n TEST pix2sr ')
print('-------------')
print('1 pix = {} sr'.format(pix2sr(1, cdelt)))

print('\n TEST sr2arcsec2 ')
print('-----------------')
print('1 sr = {} arcsec^2'.format(sr2arcsec2(1)))

print('\n TEST rad2arcsec ')
print('-----------------')
print('1 rad = {} arcsec'.format(rad2arcsec(1)))

print('\n TEST hour2deg ')
print('---------------')
ra, dec = hour2deg(9, 55, 52.725, 69, 40, 45.78)
print('09h55m52.725s, +69d40m45.78s = RA: {0} DEC: {1}'.format(ra, dec))

print('\n TEST deg2hour ')
print('---------------')
h, m, s, deg, arcmin, arcsec = deg2hour(148.969687, 69.679383)
print('RA: 148.969687, DEC: 69.679383 = {}h{}m{}s, {}d{}m{}s'.format(h, m, s, deg, arcmin, arcsec))

print('\n TEST get_cd ')
print('-------------')
get_cd_0 = get_cd(pc, cdelt) # via pc+cdelt
print('cd (pc+cdelt)\n', get_cd_0.cd, '\n')

hdr = fits.open(datdir+'M83.fits')[0].header
get_cd_1 = get_cd(header=hdr) # via header
print('pc (header)\n', get_cd_1.pc)
print('cdelt (header)\n', get_cd_1.cdelt)
print('cd (header)\n', get_cd_1.cd, '\n')

w = fixwcs(datdir+'M83').wcs
get_cd_2 = get_cd(wcs=w) # via wcs
print('pc (wcs)\n', get_cd_2.pc)
print('cdelt (wcs)\n', get_cd_2.cdelt)
print('cd (wcs)\n', get_cd_2.cd)

print('\n TEST get_pc ')
print('-------------')
print('pc\n', get_pc(cd).pc)

'''
print('\n TEST fixwcs ')
print('-------------')
# print(WCS(datdir+'M82_09_LL2.fits')) # not working

# hdul = fits.open('data/M83_LL1.fits')
hdr = fits.open(datdir+'M82_09_SL2.fits')[0].header
header = hdr.copy()
# for kw in hdr.keys():
	# if '3' in kw:
		# del header[kw]
		# print(kw)

header['NAXIS'] = 2 # SIP kw sensible
del header['NAXIS3']
del header['PC3_3'] # SIP kw sensible
del header['CRPIX3'] # SIP kw sensible
del header['CRVAL3'] # SIP kw sensible
del header['CTYPE3'] # SIP kw sensible
del header['CUNIT3'] # SIP kw sensible
del header['PS3_0'] # SIP kw sensible
del header['PS3_1'] # SIP kw sensible
# del header['FLXCON03']
# del header['FLXERR03']
# del header['GAIN3']
# del header['BASECH3']
print(WCS(header))
# print(list(header.keys()))

# print(fixwcs(datdir+'M82_09_LL2').header)
print(fixwcs().wcs)
'''
