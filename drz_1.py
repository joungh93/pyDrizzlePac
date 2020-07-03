#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 16:23:45 2019

@author: jlee
"""


import numpy as np
import glob, os
import time
from drizzlepac import astrodrizzle
import init_param as ip

start_time = time.time()


# ----- Running astrodrizzle task ----- #
current_dir = os.getcwd()

os.chdir(ip.dir_twk)

dir_drz = []
dir_drz.append('drz_'+ip.ref_flt)
glob_drz = glob.glob('drz_*')
glob_drz.remove('drz_'+ip.ref_flt)
for i in np.arange(len(glob_drz)):
	dir_drz.append(glob_drz[i])

for di in dir_drz:
	flt = di.split('_')[1]
	os.chdir(di)
	os.system('rm -rfv final*.fits')
	
	if (di == dir_drz[0]):
		astrodrizzle.AstroDrizzle('@input_'+flt+'.list', preserve=False, driz_sep_kernel='gaussian', driz_sep_pixfrac=1.0,
			                      combine_type='minmed', combine_nlow=0, combine_nhigh=1, final_kernel='gaussian', final_scale=ip.pixscl1,
			                      skymethod='globalmin+match', driz_sep_bits=32, driz_cr_scale='1.5 1.2', final_bits=352, final_wcs=True, final_rot=360)
	else:
		os.system('cp -rpv ../../'+ip.dir_out+ip.ref_flt+'.fits .')
		astrodrizzle.AstroDrizzle('@input_'+flt+'.list', preserve=False, driz_sep_kernel='gaussian', driz_sep_pixfrac=1.0,
			                      combine_type='minmed', combine_nlow=0, combine_nhigh=1, final_kernel='gaussian', final_scale=ip.pixscl1,
			                      skymethod='globalmin+match', driz_sep_bits=32, driz_cr_scale='1.5 1.2', final_bits=352, final_wcs=True, final_rot=360,
			                      final_refimage=ip.ref_flt+'.fits')

	os.system('rm -rfv final_dr*_ctx.fits final_med.fits *single* *mask*.fits *blt.fits *Mask.fits')
	os.system('cp -rpv final_dr*_sci.fits ../../'+ip.dir_out+flt+'.fits')
	os.chdir('../')


# ----- After astrodrizzle task ----- #
os.chdir(current_dir)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
