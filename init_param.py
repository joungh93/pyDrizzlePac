#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 23:23:03 2020

@author: jlee
"""


import numpy as np
import glob, os


# ----- Raw images ----- #
dir_img = '../Raw/'    # Raw data directory
img_name = glob.glob(dir_img+'*.fits')    # Names of raw images
date_order = True    # Order of image (the first image will be the reference.)


# # ----- Cosmic ray masking parameters ----- #
# cr_thre = 4000.0
# cr_mask = -20000.0


# ----- Initial photometric parameters ----- #
r_ap1, r_ap2 = 3.0, 6.0    # Radii of apertures for photometry
detect_thresh = 1.5    # Detection threshold (sigma)
tolerance = 1.5    # Matching tolerance & search radius (arcsec)


# ----- Parameter cut (these should be revised interactively!) ----- #
err_cut = 0.40    # Magnitude error cut
cidx_lo_cut = 0.05    # Lower C index cut
cidx_hi_cut = 0.60    # Higher C index cut
mag_lo_cut = 17.0    # Rrighter magnitude cut
mag_hi_cut = 27.0    # Fainter magnitude cut


# ----- Initial setting for running drizzlepac ----- #
ref_flt = '814'    # Reference filter
pixscl1 = 0.05    # arcsec/pix (final)
dir_twk = 'tweak/'    # The directory for running TWEAKREG
dir_out = 'Results/'    # Output image directory
n_minobj = 15    # Minimum number of objects for running TWEAKREG task

