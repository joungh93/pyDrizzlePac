#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 16:23:45 2019

@author: jlee
"""


import numpy as np
import glob, os
import time
from drizzlepac import tweakreg
import init_param as ip


start_time = time.time()


### Note that you should set astroconda environment before running the code!
### Note that this code is for drizzling of only two bands.


# ----- Initial structure setting ----- #
current_dir = os.getcwd()

comb_name = "com"
comb_lst = "input_"+comb_name+".list"
comb_cat = "catalog_"+comb_name+".list"

os.system("rm -rfv "+ip.dir_twk+comb_lst)
os.system("rm -rfv "+ip.dir_twk+comb_cat)
os.system("rm -rfv "+ip.dir_twk+"tweakreg*.log")
os.system("rm -rfv "+ip.dir_twk+"*_flc.fits")

print('# ----- Total ----- #')
os.system('cp -rpv '+ip.dir_img+'*.fits '+ip.dir_twk)
os.chdir(ip.dir_twk)

for inp in glob.glob('input_*.list'):
	flt = inp.split('input')[1].split('.')[0]
	print('# ----- filter '+flt[1:]+' ----- #')
	flc = np.genfromtxt(inp, dtype=None, encoding='ascii')
	os.system('cp -rpv ../Phot'+flt+'/*.mat.coo .')
	if (flt[1:] == ip.ref_flt):
		ref_img = flc[0]
	else:
		nref_flt = flt[1:]
		nref_img = flc

ref_line = []
nref_line = []
for cat in glob.glob('catalog_*.list'):
	flt = cat.split('catalog')[1].split('.')[0]
	g = open(cat, 'r')
	if (flt[1:] == ip.ref_flt):
		ref_line.append(g.readline())
	else:
		for i in np.arange(len(nref_img)):
			nref_line.append(g.readline())


# ----- Writing combined list & catalog files ----- #
f = open(comb_lst, "w")
f.write(ref_img+'\n')
for i in np.arange(len(nref_img)):
	f.write(nref_img[i]+'\n')
f.close()

f = open(comb_cat, "w")
f.write(ref_line[0])
for i in np.arange(len(nref_img)):
	f.write(nref_line[i])
f.close()


# ----- Running tweakreg task ----- #
twk_order = [ip.ref_flt, comb_lst.split('_')[1].split('.')[0], nref_flt]
for twk in twk_order:
	tweakreg.TweakReg('@input_'+twk+'.list', interactive=False, updatehdr=True,
		              runfile='tweakreg_'+twk+'.log', wcsname='TWEAK_'+twk,
		              catfile='catalog_'+twk+'.list',
		              xcol=1, ycol=2, fluxcol=3, minobj=15, searchrad=1.5,
					  tolerance=2.0)


# ----- After tweakreg task ----- #
os.system('rm -rfv *catalog.coo *.match *xy*.list')

# Copying header-updated images to drizzle path
for inp in glob.glob('input_*.list'):
	flt = inp.split('input')[1].split('.')[0]
	if not (flt == '_'+comb_name):
		print('# ----- filter '+flt[1:]+' ----- #')
		flc = np.genfromtxt(inp, dtype=None, encoding='ascii')
		for file in flc:
			os.system('cp -rpv '+file+' drz'+flt+'/')
			
os.chdir(current_dir)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
