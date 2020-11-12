#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:43:27 2019

@author: jlee
"""


import numpy as np
import glob, os
import time
from astropy.io import fits
from astroscrappy import detect_cosmics
import init_param as ip


start_time = time.time()


# ========== Code overview ========== #
# Startup settings for drizzling of HST ACS, WFC3/IR, WFC3/UVIS images
# Creating cosmic ray removed images
# =================================== #


# ----- Initialization ----- #
os.system('rm -rfv Phot_* '+ip.dir_twk)
os.system('rm -rfv '+ip.dir_out)
os.system('mkdir '+ip.dir_out)


# ----- Reading instrument and filters from the headers of images ----- #
inst, filt, inst_filt = [], [], []
for i in np.arange(len(ip.img_name)):
    hdr = fits.getheader(ip.img_name[i], ext=0)
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
inst, filt, inst_filt = np.array(inst), np.array(filt), np.array(inst_filt)


# ----- Creating directories ----- #
ufilt, uidx = np.unique(filt, return_index=True)
uinst_filt = inst_filt[uidx]
nfilt = len(uidx)


os.system('mkdir '+ip.dir_twk)
for i in np.arange(nfilt):
    os.system('mkdir Phot_'+ufilt[i][1:4])    # i.e. F475W -> 475, F850LP -> 850
    os.system('mkdir '+ip.dir_twk+'drz_'+ufilt[i][1:4])


# ----- Copying images to the corresponding directories ----- #
for i in np.arange(len(ip.img_name)):
    dir_filt = [filt[j][1:4] for j in np.arange(len(filt))]
    os.system('cp -rpv '+ip.img_name[i]+' Phot_'+dir_filt[i])


# ----- Making lists & masking cosmic rays ----- #
current_dir = os.getcwd()

for i in np.arange(nfilt):
    os.chdir(current_dir+'/'+'Phot_'+ufilt[i][1:4])
    img = glob.glob('*_fl*.fits')    # *_flc.fits (ACS, WFC3/UVIS), *_flt.fits (WFC3/IR)
    date = np.array([])
    expt = np.array([])

    # Sorting with MJD
    for j in img:
        hdr = fits.getheader(j, ext=0)
        date = np.append(date, hdr['EXPSTART'])
        expt = np.append(expt, hdr['EXPTIME'])

    # Image order
    if ip.date_order:
    	order = date.argsort()
    else:
    	order = expt.argsort()[::-1]


    # Making input.list
    f1 = open('input.list', 'w')
    for j in np.arange(len(order)):
        f1.write(str(np.array(img)[order][j])+' \n')
    f1.close()
    os.system('cp -rpv input.list ../'+ip.dir_twk+'input_'+ufilt[i][1:4]+'.list')
    os.system('cp -rpv input.list ../'+ip.dir_twk+'drz_'+ufilt[i][1:4]+'/input_'+ufilt[i][1:4]+'.list')

    # Making catalog.list ---> ACS, WFC3/UVIS 2 CCD chips & WFC3/IR 1 CCD chip
    f2 = open('catalog.list', 'w')
    for j in np.arange(len(order)):
        if (uinst_filt[i].split('/')[0] == 'ACS'):
            f2.write(str(np.array(img)[order][j])+' f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_sci1.mat.coo'+' f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_sci2.mat.coo'+' \n')
        if (uinst_filt[i].split('/')[0] == 'WFC3'):
            f2.write(str(np.array(img)[order][j])+' f'+ufilt[i][1:4]+'_%02d' %(j+1)+'_sci1.mat.coo'+' \n')
    f2.close()
    os.system('cp -rpv catalog.list ../'+ip.dir_twk+'catalog_'+ufilt[i][1:4]+'.list')

    # Creating cosmic ray removed images ---> ACS, WFC3/UVIS 4 extensions & WFC3/IR 2 extensions
    os.chdir(current_dir+'/'+'Phot_'+ufilt[i][1:4])
    f = open('image_names.log','w')    
    for j in np.arange(len(order)):
        hd0 = fits.getheader(str(np.array(img)[order][j]), ext=0)
        if (uinst_filt[i].split('/')[0] == 'ACS'):
            ext_list = [1, 3, 4, 6]
        if (uinst_filt[i].split('/')[0] == 'WFC3'):
            ext_list = [1, 3]
        for k in np.arange(len(ext_list)):
            if (k % 2 == 0):
                print("Writing "+'f'+ufilt[i][1:4]+'_%02d_sci%1d.fits' %(j+1, 1+k//2)+"...")
                sci, hdr = fits.getdata(str(np.array(img)[order][j]), ext=ext_list[k], header=True)
                # dq = fits.getdata(str(np.array(img)[order][j]), ext=ext_list[k+1], header=False)
                # cr = (dq >= ip.cr_thre)
                # if (np.sum(cr) >= 1):
                #     sci[cr] = ip.cr_mask
                epadu = hd0['CCDGAIN']
                crmask, cleanarr = detect_cosmics(sci, gain=epadu, cleantype='medmask')
                fits.writeto('f'+ufilt[i][1:4]+'_%02d_sci%1d.fits' %(j+1, 1+k//2),
                             cleanarr/epadu, hdr, overwrite=True)
                f.write(str(np.array(img)[order][j]))
                f.write('\t')
                f.write('f'+ufilt[i][1:4]+'_%02d_sci%1d.fits' %(j+1, 1+k//2))
                f.write('\n')
            else:
                continue
    f.close()

os.chdir(current_dir)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
