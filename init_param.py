#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 23:23:03 2020

@author: jlee
"""


import numpy as np
import glob, os


# ----- Raw images ----- #
dir_img = 'Images/'
img_name = glob.glob(dir_img+'*.fits')


# ----- Cosmic ray masking parameters ----- #
cr_thre = 4000.0
cr_mask = -20000.0


# ----- Initial photometric parameters ----- #
pixscl0 = 0.05    # arcsec/pix (initial)
pixscl1 = 0.05    # arcsec/pix (final)
zmag = 25.0
gain = 2.0
fwhm = 2.0
r_ap1, r_ap2 = 3.0, 6.0
detect_thresh = 1.5


# ----- Parameter cut ----- #
err_cut = 0.20    # magnitude error cut
cidx_lo_cut = 0.05    # lower C index cut
cidx_hi_cut = 0.30    # higher C index cut
mag_lo_cut = 16.5    # brighter magnitude cut
mag_hi_cut = 23.0    # fainter magnitude cut


# ----- Initial setting for running drizzlepac ----- #
ref_flt = '850'
dir_twk = 'tweak/'
dir_out = 'Results/'