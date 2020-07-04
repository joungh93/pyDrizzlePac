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


# ----- Initial structure setting ----- #
current_dir = os.getcwd()

comb_name = "com"
comb_lst = "input_"+comb_name+".list"
comb_cat = "catalog_"+comb_name+".list"

os.system("rm -rfv "+ip.dir_twk+comb_lst)
os.system("rm -rfv "+ip.dir_twk+comb_cat)
os.system("rm -rfv "+ip.dir_twk+"tweakreg*.log")
os.system("rm -rfv "+ip.dir_twk+"*_fl*.fits")
os.system("rm -rfv "+ip.dir_twk+"*.png")

print('# ----- Total ----- #')
os.system('cp -rpv '+ip.dir_img+'*.fits '+ip.dir_twk)
os.chdir(ip.dir_twk)

lst_input = glob.glob('input_*.list')
lst_input = sorted(lst_input)

# Reference filter & Non-reference filters
flt = [b.split('input_')[1].split('.')[0] for b in lst_input]
flt = np.array(flt)
non_ref = flt[flt != ip.ref_flt]

# Copying *.mat.coo to this directory
for i in np.arange(len(flt)):
	os.system('cp -rpv ../Phot_'+flt[i]+'/*.mat.coo .')

# Reference data from input/catalog file
ref_imgs = np.genfromtxt('input_'+ip.ref_flt+'.list',
	                     dtype=None, encoding='ascii')
g = open('catalog_'+ip.ref_flt+'.list','r')
ref_lines = g.readlines()
g.close()


# ----- Writing combined list & catalog files ----- #
for i in np.arange(len(non_ref)):
	f = open(comb_lst, "w")
	f.write(ref_imgs[0]+'\n')
	nref_imgs = np.genfromtxt('input_'+non_ref[i]+'.list',
	                          dtype=None, encoding='ascii')
	for j in np.arange(len(nref_imgs)):
		f.write(nref_imgs[j]+'\n')
	f.close()

	f = open(comb_cat, "w")
	f.write(ref_lines[0])
	g = open('catalog_'+non_ref[i]+'.list')
	nref_lines = g.readlines()
	g.close()
	for j in np.arange(len(nref_lines)):
		f.write(nref_lines[j])
	f.close()

	print("... TWEAKREG task for "+ip.ref_flt+"(ref), "+non_ref[i]+" ...")

	# ----- Running tweakreg task ----- #
	twk_order = [ip.ref_flt, comb_name, non_ref[i]]
	for twk in twk_order:
		tweakreg.TweakReg('@input_'+twk+'.list', interactive=False, updatehdr=True,
			              runfile='tweakreg_'+twk+'.log', wcsname='TWEAK_'+twk+'_'+str(i+1),
			              catfile='catalog_'+twk+'.list',
			              xcol=1, ycol=2, fluxcol=3, minobj=ip.n_minobj, searchrad=ip.tolerance,
						  tolerance=2.0)		


# ----- After tweakreg task ----- #
os.system('rm -rfv *catalog.coo *.match *xy*.list')

# Copying header-updated images to drizzle path
for inp in lst_input:
	flt = inp.split('input')[1].split('.')[0]
	if not (flt == '_'+comb_name):
		print('# ----- filter '+flt[1:]+' ----- #')
		flc = np.genfromtxt(inp, dtype=None, encoding='ascii')
		for file in flc:
			os.system('cp -rpv '+file+' drz'+flt+'/')
		
os.chdir(current_dir)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
