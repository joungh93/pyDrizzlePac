#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:43:27 2019

@author: jlee
"""


import numpy as np
import glob
import os
import time
from astropy.io import fits
from pyraf import iraf


start_time = time.time()


# ========== Code overview ========== #
# Startup settings for drizzling of HST ACS, WFC3/IR, WFC3/UVIS images
# =================================== #


# ----- Initialization ----- #
os.system('rm -rfv Phot_* tweak')


# ----- Raw images ----- #
img_dir = 'Images/'
img_name = glob.glob(img_dir+'*.fits')


# ----- Reading instrument and filters from the headers of images ----- #
inst, filt, inst_filt = [], [], []
for i in np.arange(len(img_name)):
    hdr = fits.getheader(img_name[i], ext=0)
    inst.append(hdr['INSTRUME'])
    if (hdr['INSTRUME'] == 'ACS'):
        if (hdr['FILTER1'][0:5] != 'CLEAR'):
            hfilt = hdr['FILTER1']
            filt.append(hfilt)
        if (hdr['FILTER2'][0:5] != 'CLEAR'):
            hfilt = hdr['FILTER2']
            filt.append(hfilt)
    if (hdr['INSTRUME'] == 'WFC3'):
        hfilt = hdr['FILTER']
        filt.append(hfilt)
    inst_filt.append(hdr['INSTRUME']+'/'+hfilt)


# ----- Creating directories ----- #
ufilt, uinst_filt = list(set(filt)), list(set(inst_filt))
nfilt = len(ufilt)

os.system('mkdir tweak')
for i in np.arange(nfilt):
    os.system('mkdir Phot_'+ufilt[i][1:4])    # F475W -> 475, F850LP -> 850
    os.system('mkdir tweak/drz_'+ufilt[i][1:4])


# ----- Copying images to the corresponding directories ----- #
for i in np.arange(len(img_name)):
    dir_filt = [filt[j][1:4] for j in np.arange(len(filt))]
    os.system('cp -rpv '+img_name[i]+' Phot_'+dir_filt[i])


# ----- Making lists & Imcopy tasks ----- #
current_dir = os.getcwd()

for i in np.arange(nfilt):
    os.chdir(current_dir+'/'+'Phot_'+ufilt[i][1:4])
    img = glob.glob('*_fl*.fits')    # *_flc.fits (ACS), *_flt.fits (WFC3/IR)
    date = np.array([])

    # Sorting with MJD
    for j in img:
        hdr = fits.getheader(j, ext=0)
        date = np.append(date, hdr['EXPSTART'])
    order = date.argsort()
    
    # Making input.list
    f1 = open('input.list', 'w')
    for j in np.arange(len(order)):
        f1.write(str(np.array(img)[order][j])+' \n')
    f1.close()
    os.system('cp -rpv input.list ../tweak/input_'+ufilt[i][1:4]+'.list')
    os.system('cp -rpv input.list ../tweak/'+'drz_'+ufilt[i][1:4]+'/input_'+ufilt[i][1:4]+'.list')

    # Making catalog.list ---> ACS 2 CCD chips & WFC3/IR 1 CCD chip
    f2 = open('catalog.list', 'w')
    for j in np.arange(len(order)):
        if (uinst_filt[i].split('/')[0] == 'ACS'):
            f2.write(str(np.array(img)[order][j])+' f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_sci1.mat.coo'+' f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_sci2.mat.coo'+' \n')
        if (uinst_filt[i].split('/')[0] == 'WFC3'):
            f2.write(str(np.array(img)[order][j])+' f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_sci1.mat.coo'+' \n')
    f2.close()
    os.system('cp -rpv catalog.list ../tweak/catalog_'+ufilt[i][1:4]+'.list')

    # Imcopy tasks (to separate extensions) ---> ACS 4 extensions & WFC3/IR 2 extensions
    iraf.chdir(current_dir+'/'+'Phot_'+ufilt[i][1:4])    # Check iraf directory 'iraf.pwd()'    
    for j in np.arange(len(order)):
        if (uinst_filt[i].split('/')[0] == 'ACS'):
            iraf.imcopy(input=str(np.array(img)[order][j])+'[1]',
                        output='f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_temp1.fits')    # [SCI]
            iraf.imcopy(input=str(np.array(img)[order][j])+'[3]',
                        output='f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_temp2.fits')    # [DQ]
            iraf.imcopy(input=str(np.array(img)[order][j])+'[4]',
                        output='f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_temp3.fits')    # [SCI]
            iraf.imcopy(input=str(np.array(img)[order][j])+'[6]',
                        output='f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_temp4.fits')    # [DQ]
        if (uinst_filt[i].split('/')[0] == 'WFC3'):
            iraf.imcopy(input=str(np.array(img)[order][j])+'[1]',
                        output='f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_temp1.fits')    # [SCI]
            iraf.imcopy(input=str(np.array(img)[order][j])+'[3]',
                        output='f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_temp2.fits')    # [DQ]


os.chdir(current_dir)
iraf.chdir(current_dir)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
