#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:49:33 2019

@author: jlee
"""


import numpy as np
import glob, os
import time
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clip
from photutils import MMMBackground
from photutils import StdBackgroundRMS
import sep
# from photutils import DAOStarFinder as DSF
# from photutils import IRAFStarFinder as ISF
# from photutils.aperture import CircularAperture as CAp
# from photutils.aperture import CircularAnnulus as CAn
# from photutils import aperture_photometry as apphot
from stwcs.wcsutil import hstwcs
import init_param as ip


start_time = time.time()


# ========== Code overview ========== #
# Sky & sky sigma estimation
# Aperture photometry with Python/sep
# =================================== #


# ----- Basic settings ----- #
current_dir = os.getcwd()

for w_dir in glob.glob('Phot_*'):
    w_dir = w_dir+'/'
    phot_dir = w_dir

    print('#---------- Working on '+phot_dir+' ----------#')

    # ----- Reading image lists ----- #
    imglist = np.genfromtxt(phot_dir+'image_names.log', dtype=None,
                            encoding='ascii', names=('raw','sci'))
    nshot = len(np.unique(imglist['raw']))
    nchip = len(imglist) // nshot


    # ----- The main loop ----- #
    os.system('rm -rfv '+phot_dir+'opt_1.dat')
    f = open(phot_dir+'opt_1.dat','w')
    f.write('# shot  sky  sig \n')  
    for i in np.arange(nshot):
        for j in np.arange(nchip):
            # Reading images
            imgidx = nchip*i+j
            sci_name = imglist['sci'][imgidx].split('.fits')[0]

            dat = fits.getdata(phot_dir+imglist['sci'][imgidx], header=False)
            var = (dat > ip.cr_mask)
            msk = np.zeros(dat.shape, dtype=bool)
            msk[~var] = True
            nsky = np.sum(msk)

            h0 = fits.getheader(phot_dir+imglist['raw'][imgidx], ext=0)
            w = hstwcs.HSTWCS(phot_dir+imglist['raw'][imgidx]+'[SCI,%d]' %(j+1))

            # Sky & sky sigma estimation
            dat = dat.byteswap().newbyteorder()
            bkg = sep.Background(dat, mask=None, bw=64, bh=64, fw=3, fh=3)
            sky, sig = bkg.globalback, bkg.globalrms
            print(sci_name+' : %.2f  %.2f' %(sky, sig))
            f.write('%d  %.2f  %.2f \n' %(i+1, sky, sig))

            # Extracting sources
            dat_sub = dat - bkg
            src = sep.extract(dat_sub, ip.detect_thresh, err=bkg.globalrms)
            n_src = len(src)
            print("{0:d} sources are found.".format(n_src))

            g = open(phot_dir+'src_'+sci_name+'.reg','w')
            for k in np.arange(n_src):
                g.write('{0:.4f}  {1:.4f}\n'.format(src['x'][k]+1, src['y'][k]+1))
            g.close()

            # Aperture photmetry
            flx1, e_flx1, flag1 = sep.sum_circle(dat_sub, src['x'], src['y'],
                                                 ip.r_ap1, err=bkg.globalrms,
                                                 gain=ip.gain, subpix=0)
            flx2, e_flx2, flag2 = sep.sum_circle(dat_sub, src['x'], src['y'],
                                                 ip.r_ap2, err=bkg.globalrms,
                                                 gain=ip.gain, subpix=0)
            itime, exptime = h0['EXPTIME'], h0['EXPTIME']
            nflx1, nflx2 = flx1/itime, flx2/itime
            tflx1, tflx2 = nflx1*exptime, nflx2*exptime
            mag1 = ip.zmag - 2.5*np.log10(nflx1)
            mag2 = ip.zmag - 2.5*np.log10(nflx2)
            e_mag1 = (2.5/np.log(10.0)) * (e_flx1/tflx1)
            e_mag2 = (2.5/np.log(10.0)) * (e_flx2/tflx2)

            # XY to WCS
            ra, dec = w.all_pix2world(src['x']+1, src['y']+1, 1)

            # Saving the results
            df = pd.DataFrame(data = {'x' : src['x']+1,
                                      'y' : src['y']+1,
                                      'mag1' : mag1,
                                      'mag2' : mag2,
                                      'e_mag1' : e_mag1,
                                      'e_mag2' : e_mag2,
                                      'flx1' : nflx1,
                                      'flx2' : nflx2,
                                      'e_flx1' : e_flx1/itime,
                                      'e_flx2' : e_flx2/itime,
                                      'ra' : ra,
                                      'dec' : dec})
            df = df.fillna(99.0)
            df.to_pickle(phot_dir+sci_name+'.pkl')

    f.close()


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
