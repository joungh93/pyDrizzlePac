#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:49:33 2019

@author: jlee
"""


import numpy as np
import glob
import os
import time
from astropy.io import fits
from astropy.stats import SigmaClip
from photutils import MMMBackground
from photutils import StdBackgroundRMS


start_time = time.time()


# ========== Code overview ========== #
# 1. Making cosmic raw removed images
# 2. Sky & sky sigma estimation
# 3. SExtractor
# =================================== #


current_dir = os.getcwd()


# ----- SExtractor initial parameters ----- #
d_aper = "3.0,4.0,6.0"
zmag = "25.0"
pixscl = "0.05"
seeing = "0.10"


for w_dir in glob.glob('Phot_*'):

    w_dir = w_dir+'/'
    
    print('#---------- Working on '+w_dir+' ----------#')

    # ----- Temporary images for photometry ----- #
    phot_dir = w_dir

    tmp_name = glob.glob(phot_dir+'f*temp*.fits')

    N_tmp1 = np.array([10*int(tmp_name[i].split('/')[1].split('_')[1]) for i in np.arange(len(tmp_name))])
    N_tmp2 = np.array([int(tmp_name[i].split('.')[0][-1]) for i in np.arange(len(tmp_name))])
    N_tmp = N_tmp1 + N_tmp2

    tmp_order = np.argsort(N_tmp)

    tmp_name = np.array(tmp_name)[tmp_order]


    # ----- 1. Making cosmic ray removed images ----- #
    os.system('rm -rfv '+phot_dir+'f*sci*.fits')

    nele = int(len(tmp_name)/4)

    for i in np.arange(nele):
        for j in [0, 1]:
            img1, hdr1 = fits.getdata(tmp_name[4*i+2*j], header=True)
            img2, hdr2 = fits.getdata(tmp_name[4*i+2*j+1], header=True)

            var = (img2 >= 4000.0)
            if (np.sum(var) >= 1):
                img1[var] = -20000

            fits.writeto(tmp_name[4*i+2*j].split('temp')[0]+'sci'+'%1d' %(j+1)+'.fits',
                         img1, hdr1, overwrite=True)


    # ----- Science images for photometry ----- #
    sci_name = glob.glob(phot_dir+'f*sci*.fits')

    N_sci1 = np.array([10*int(sci_name[i].split('/')[1].split('_')[1]) for i in np.arange(len(sci_name))])
    N_sci2 = np.array([int(sci_name[i].split('.')[0][-1]) for i in np.arange(len(sci_name))])
    N_sci = N_sci1 + N_sci2

    sci_order = np.argsort(N_sci)

    sci_name = np.array(sci_name)[sci_order]


    # ----- 2. Sky & sky sigma estimation ----- #
    sigma_clip = SigmaClip(sigma = 3.0)
    bkg = MMMBackground(sigma_clip = sigma_clip)
    bkgrms = StdBackgroundRMS(sigma_clip)

    nshot = int(len(sci_name)/2)

    os.system('rm -rfv '+phot_dir+'opt_1.dat')
    f = open(phot_dir+'opt_1.dat','w')
    f.write('#  shot  sky  sig \n')

    for i in np.arange(nshot):
        img1, hdr1 = fits.getdata(sci_name[2*i]), fits.getheader(sci_name[2*i])
        var1 = (img1 > -20000.0)
        sky1, sig1 = bkg(img1[var1]), bkgrms(img1[var1])
        print(sci_name[2*i].split('/')[1]+' : %.2f  %.2f' %(sky1, sig1))
        f.write('%d  %.2f  %.2f \n' %(i+1, sky1, sig1))

        img2, hdr2 = fits.getdata(sci_name[2*i+1]), fits.getheader(sci_name[2*i+1])
        var2 = (img2 > -20000.0)
        sky2, sig2 = bkg(img2[var2]), bkgrms(img2[var2])
        print(sci_name[2*i+1].split('/')[1]+' : %.2f  %.2f' %(sky2, sig2))
        f.write('%d  %.2f  %.2f \n' %(i+1, sky2, sig2))

    f.close()


    # ----- 3. SExtractor ----- #
    f = open(phot_dir+'sephot.sh', 'w')
    for i in np.arange(len(sci_name)):
        sci_img = sci_name[i].split('/')[1]
        txt1 = 'sex '+sci_img+' -c '+current_dir+'/npsf.sex'
        txt2 = ' -CATALOG_NAME '+sci_img.split('.fits')[0]+'.cat -DETECT_MINAREA 4 -PHOT_APERTURES '+d_aper
        txt3 = ' -MAG_ZEROPOINT '+zmag+' -PIXEL_SCALE '+pixscl+' -SEEING_FWHM '+seeing
        f.write(txt1+txt2+txt3+'\n')
    f.close()

    os.chdir(phot_dir)
    os.system('sh sephot.sh')
    os.chdir(current_dir)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
