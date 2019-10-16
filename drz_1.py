#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 16:23:45 2019

@author: jlee
"""


import numpy as np
import glob
import os
import time
import drizzlepac
from drizzlepac import astrodrizzle

start_time = time.time()


### Note that you should set astroconda environment before running the code!
### Note that this code is for drizzling of only two bands.


# ----- Running astrodrizzle task ----- #
current_dir = os.getcwd()
dir_twk = 'tweak/'

os.chdir(dir_twk)

ref_flt = str(input("Enter reference filter name (i.e. F814W -> type '814'): "))

dir_drz = []
dir_drz.append('drz_'+ref_flt)
glob_drz = glob.glob('drz_*')
glob_drz.remove('drz_'+ref_flt)
for i in np.arange(len(glob_drz)):
	dir_drz.append(glob_drz[i])

for di in dir_drz:
	flt = di.split('_')[1]
	os.chdir(di)
	os.system('rm -rfv final*.fits')
	
	if (di == dir_drz[0]):
		astrodrizzle.AstroDrizzle('@input_'+flt+'.list', preserve=False, driz_sep_kernel='gaussian', driz_sep_pixfrac=1.0,
			                      combine_type='minmed', combine_nlow=0, combine_nhigh=1, final_kernel='gaussian', final_scale=0.05,
			                      skymethod='globalmin+match', driz_sep_bits=32, driz_cr_scale='1.5 1.2', final_bits=352, final_wcs=True, final_rot=360)
		os.system('rm -rfv final_drc_ctx.fits final_med.fits *single* *mask*.fits *blt.fits *Mask.fits')
		os.system('cp -rpv final_drc_sci.fits ../'+flt+'.fits')
		os.chdir('../')
	else:
		os.system('cp -rpv ../'+ref_flt+'.fits .')
		astrodrizzle.AstroDrizzle('@input_'+flt+'.list', preserve=False, driz_sep_kernel='gaussian', driz_sep_pixfrac=1.0,
			                      combine_type='minmed', combine_nlow=0, combine_nhigh=1, final_kernel='gaussian', final_scale=0.05,
			                      skymethod='globalmin+match', driz_sep_bits=32, driz_cr_scale='1.5 1.2', final_bits=352, final_wcs=True, final_rot=360,
			                      final_refimage=ref_flt+'.fits')
		os.system('rm -rfv final_drc_ctx.fits final_med.fits *single* *mask*.fits *blt.fits *Mask.fits')
		os.system('cp -rpv final_drc_sci.fits ../'+flt+'.fits')
		os.chdir('../')


# ----- After astrodrizzle task ----- #
os.chdir(current_dir)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))