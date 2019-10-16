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
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy import units as u
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


start_time = time.time()


# ========== Code overview ========== #
# Step 1 : Plotting & selecting point sources with cut conditions
# Step 2 : Writing catalogs and regions of selected point sources
### *.cut.dat : x, y, flux
### *.cut.coo : ra, dec (distortion not corrected)
### *.cut.reg : x, y (DS9)
# Step 3 : Matching point sources from each field and creating *.coo
### *.mat.coo : x, y, flux
### *.mat.reg : x, y (DS9)
# =================================== #


for w_dir in glob.glob('Phot_*'):

    w_dir = w_dir+'/'
   
    print(' \n')
    print('#---------- Working on '+w_dir+' ----------#')

    
    # ========== Parameter cut ========== #
    err_cut = float(input("Enter magnitude error cut (default : 0.20) : ") or "0.20")
    cidx_cut1 = float(input("Enter higher C index cut (default : 0.60) : ") or "0.60")
    cidx_cut2 = float(input("Enter lower C index cut (default : 0.30) : ") or "0.30")
    mag_cut1 = float(input("Enter brighter magnitude cut (default : 18.5) : ") or "18.5")
    mag_cut2 = float(input("Enter fainter magnitude cut (default : 22.5) : ") or "22.5")
    # =================================== #
    
    os.system('rm -rfv '+w_dir+'mat_1.log')
    l = open(w_dir+'mat_1.log','w')
    l.write('Magnitude error cut : %.2f \n' %(err_cut))
    l.write('Higher concentration index cut : %.2f \n' %(cidx_cut1))
    l.write('Lower concentration index cut : %.2f \n' %(cidx_cut2))
    l.write('Brighter magnitude cut : %.2f \n' %(mag_cut1))
    l.write('Fainter magnitude cut : %.2f \n' %(mag_cut2))
    l.close()
    
    # ----- SExtractor photometry catalog ----- #
    phot_dir = w_dir

    cat_name = glob.glob(phot_dir+'f*sci*.cat')

    N_cat1 = np.array([10*int(cat_name[i].split('/')[1].split('_')[1]) for i in np.arange(len(cat_name))])
    N_cat2 = np.array([int(cat_name[i].split('.')[0][-1]) for i in np.arange(len(cat_name))])
    N_cat = N_cat1 + N_cat2

    cat_order = np.argsort(N_cat)

    cat_name = np.array(cat_name)[cat_order]

    nshot = int(len(cat_name) / 2)
    
    
    # ====================================== #
    # =============== Step 1 =============== #
    # ====================================== #
            
    # ----- Plotting (PDF) ----- #

    os.system('rm -rfv '+phot_dir+'mat_1.pdf')

    with PdfPages(phot_dir+'mat_1.pdf') as pdf:
        for i in np.arange(nshot):

            # ---- Figure ---- #
            fig1 = plt.figure(1, figsize=(10,10))
            fig1.suptitle(cat_name[2*i].split('/')[1].split('.')[0][:-1], x=0.5, y=0.95, 
                          fontsize=20.0, fontweight='bold', ha='center', va='center')

            # ----- Reading catalog : SCI1 ----- # 
            print('... Reading '+cat_name[2*i].split('/')[1]+' ...')
            dat_cat1 = np.genfromtxt(cat_name[2*i], dtype=None, encoding='ascii',
                                     names=('x','y','num','mag','merr',
                                            'mag1','mag2','mag3','merr1','merr2','merr3',
                                            'kron','backgr','ra','dec','cxx','cyy','cxy','a','b','theta',
                                            'mu0','flags','fwhm','flxrad','cl'))

            hdr = fits.getheader(cat_name[2*i].split('.cat')[0]+'.fits')
            exp1 = hdr['EXPTIME']

            dat_cat1['mag'] = dat_cat1['mag'] + 2.5*np.log10(exp1)

            flx1 = 10.0 ** ((25.0-dat_cat1['mag'])/2.5)

            mag_cnd1 = ((dat_cat1['mag'] > mag_cut1) & (dat_cat1['mag'] < mag_cut2))

             
            # ----- Figure 1 : mag-merr diagram ----- #
            ax1 = fig1.add_subplot(2,2,1)
            ax1.set_position([0.10,0.52,0.42,0.40])
            ax1.set_xticks([15.0,20.0,25.0,30.0])
            ax1.set_yticks([0.0,0.2,0.4,0.6])
            ax1.tick_params(labelbottom=False)
            plt.ylabel('magnitude error',fontsize=15)
            plt.xlim([15.0,30.0]) ; plt.xticks(fontsize=15)
            plt.ylim([0.0,0.6]) ; plt.yticks(fontsize=15)
            plt.tick_params(width=1.5, length=8.0)
            plt.minorticks_on()
            plt.tick_params(width=1.5,length=5.0,which='minor')
            for axis in ['top','bottom','left','right']:
                ax1.spines[axis].set_linewidth(1.5)
            # ---------------------------------------- #
            xx = mag_cut1 + (mag_cut2-mag_cut1)*np.arange(1000)/1000.
            ax1.fill_between(xx, np.repeat(0.0, 1000), np.repeat(err_cut, 1000),
                             color='khaki', alpha=0.6)
            err_cnd1 = ((dat_cat1['merr'] > 0.0) & (dat_cat1['merr'] < err_cut))

            plt.plot(dat_cat1['mag'], dat_cat1['merr'], 'o', ms=3.5, color='gray',
                     mec='black', mew=0.8, alpha=0.7)
            plt.text(0.05, 0.95, 'SCI-1', fontsize=17.5, fontweight='bold',
                     ha='left', va='top', transform=ax1.transAxes)

            # ----- Figure 2 : mag-concentration diagram ----- #
            ax1 = fig1.add_subplot(2,2,2)
            ax1.set_position([0.10,0.10,0.42,0.40])
            ax1.set_xticks([15.0,20.0,25.0,30.0])
            ax1.set_yticks([-1.0,-0.5,0.0,0.5,1.0])
            plt.xlabel('magnitude',fontsize=15)
            plt.ylabel('C index',fontsize=15)
            plt.xlim([15.0,30.0]) ; plt.xticks(fontsize=15)
            plt.ylim([0.0,1.0]) ; plt.yticks(fontsize=15)
            plt.tick_params(width=1.5, length=8.0)
            plt.minorticks_on()
            plt.tick_params(width=1.5,length=5.0,which='minor')
            for axis in ['top','bottom','left','right']:
                ax1.spines[axis].set_linewidth(1.5)
            # -------------------------------------------- #
            xx = mag_cut1 + (mag_cut2-mag_cut1)*np.arange(1000)/1000.
            ax1.fill_between(xx, np.repeat(cidx_cut2, 1000), np.repeat(cidx_cut1, 1000),
                             color='khaki', alpha=0.6)
            
            cidx1 = dat_cat1['mag1'] - dat_cat1['mag3']
            cidx_cnd1 = ((cidx1 > cidx_cut2) & (cidx1 < cidx_cut1))

            plt.plot(dat_cat1['mag'], cidx1, 'o', ms=3.5, color='gray',
                     mec='black', mew=0.8, alpha=0.7)

            # ----- Reading catalog : SCI2 ----- # 
            print('... Reading '+cat_name[2*i+1].split('/')[1]+' ...')
            dat_cat2 = np.genfromtxt(cat_name[2*i+1], dtype=None, encoding='ascii',
                                     names=('x','y','num','mag','merr',
                                            'mag1','mag2','mag3','merr1','merr2','merr3',
                                            'kron','backgr','ra','dec','cxx','cyy','cxy','a','b','theta',
                                            'mu0','flags','fwhm','flxrad','cl'))

            hdr = fits.getheader(cat_name[2*i+1].split('.cat')[0]+'.fits')
            exp2 = hdr['EXPTIME']

            dat_cat2['mag'] = dat_cat2['mag'] + 2.5*np.log10(exp2)

            flx2 = 10.0 ** ((25.0-dat_cat2['mag'])/2.5)

            mag_cnd2 = ((dat_cat2['mag'] > mag_cut1) & (dat_cat2['mag'] < mag_cut2))

            # ----- Figure 3 : mag-merr diagram ----- #
            ax1 = fig1.add_subplot(2,2,3)
            ax1.set_position([0.54,0.52,0.42,0.40])
            ax1.set_xticks([15.0,20.0,25.0,30.0])
            ax1.set_yticks([0.0,0.2,0.4,0.6])
            ax1.tick_params(labelbottom=False)
            ax1.tick_params(labelleft=False)
            plt.xlim([15.0,30.0]) ; plt.xticks(fontsize=15)
            plt.ylim([0.0,0.6]) ; plt.yticks(fontsize=15)
            plt.tick_params(width=1.5, length=8.0)
            plt.minorticks_on()
            plt.tick_params(width=1.5,length=5.0,which='minor')
            for axis in ['top','bottom','left','right']:
                ax1.spines[axis].set_linewidth(1.5)
            # ---------------------------------------- #
            xx = mag_cut1 + (mag_cut2-mag_cut1)*np.arange(1000)/1000.
            ax1.fill_between(xx, np.repeat(0.0, 1000), np.repeat(err_cut, 1000),
                             color='khaki', alpha=0.6)

            err_cnd2 = ((dat_cat2['merr'] > 0.0) & (dat_cat2['merr'] < err_cut))

            plt.plot(dat_cat2['mag'], dat_cat2['merr'], 'o', ms=3.5, color='gray',
                     mec='black', mew=0.8, alpha=0.7)
            plt.text(0.05, 0.95, 'SCI-2', fontsize=17.5, fontweight='bold',
                     ha='left', va='top', transform=ax1.transAxes)

            # ----- Figure 4 : mag-concentration diagram ----- #
            ax1 = fig1.add_subplot(2,2,4)
            ax1.set_position([0.54,0.10,0.42,0.40])
            ax1.set_xticks([15.0,20.0,25.0,30.0])
            ax1.set_yticks([-1.0,-0.5,0.0,0.5,1.0])
            plt.xlabel('magnitude',fontsize=15)
            ax1.tick_params(labelleft=False)
            plt.xlim([15.0,30.0]) ; plt.xticks(fontsize=15)
            plt.ylim([0.0,1.0]) ; plt.yticks(fontsize=15)
            plt.tick_params(width=1.5, length=8.0)
            plt.minorticks_on()
            plt.tick_params(width=1.5,length=5.0,which='minor')
            for axis in ['top','bottom','left','right']:
                ax1.spines[axis].set_linewidth(1.5)
            # -------------------------------------------- #
            xx = mag_cut1 + (mag_cut2-mag_cut1)*np.arange(1000)/1000.
            ax1.fill_between(xx, np.repeat(cidx_cut2, 1000), np.repeat(cidx_cut1, 1000),
                             color='khaki', alpha=0.6)

            cidx2 = dat_cat2['mag1'] - dat_cat2['mag3']
            cidx_cnd2 = ((cidx2 > cidx_cut2) & (cidx2 < cidx_cut1))

            plt.plot(dat_cat2['mag'], cidx2, 'o', ms=3.5, color='gray',
                     mec='black', mew=0.8, alpha=0.7)


            pdf.savefig()
            plt.close()

            
            # ====================================== #
            # =============== Step 2 =============== #
            # ====================================== #

            # ----- Writing catalogs & regions (SCI-1) ----- # 
            os.system('rm -rfv '+phot_dir+cat_name[2*i].split('/')[1].split('.')[0]+'.cut.dat')
            os.system('rm -rfv '+phot_dir+cat_name[2*i].split('/')[1].split('.')[0]+'.cut.coo')
            os.system('rm -rfv '+phot_dir+cat_name[2*i].split('/')[1].split('.')[0]+'.cut.reg')
            
            poi1 = (1*mag_cnd1 + 1*err_cnd1 + 1*cidx_cnd1 == 3)
            print('\n # of point sources (SCI-1) : %d \n' %(np.sum(poi1)))
            
            print('... Writing '+cat_name[2*i].split('/')[1].split('.')[0]+'.cut.dat ...')
            f1 = open(phot_dir+cat_name[2*i].split('/')[1].split('.')[0]+'.cut.dat','w')
            for k in np.arange(np.sum(poi1)):
                f1.write('%.2f   %.2f   %.2f \n' \
                         %(dat_cat1['x'][poi1][k], dat_cat1['y'][poi1][k], flx1[poi1][k]))
            f1.close()

            print('... Writing '+cat_name[2*i].split('/')[1].split('.')[0]+'.cut.coo ...')
            g1 = open(phot_dir+cat_name[2*i].split('/')[1].split('.')[0]+'.cut.coo','w')
            for k in np.arange(np.sum(poi1)):
                g1.write('%.6f   %.6f \n' \
                         %(dat_cat1['ra'][poi1][k], dat_cat1['dec'][poi1][k]))
            g1.close()
            
            print('... Writing '+cat_name[2*i].split('/')[1].split('.')[0]+'.cut.reg ...')
            h1 = open(phot_dir+cat_name[2*i].split('/')[1].split('.')[0]+'.cut.reg','w')
            h1.write('global color=green font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source width=2 \n')
            for k in np.arange(np.sum(poi1)):
                h1.write('image;circle(%.2f, %.2f, 20) \n' \
                         %(dat_cat1['x'][poi1][k], dat_cat1['y'][poi1][k]))
            h1.close()

            # ----- Writing catalogs & regions (SCI-2) ----- # 
            os.system('rm -rfv '+phot_dir+cat_name[2*i+1].split('/')[1].split('.')[0]+'.cut.dat')
            os.system('rm -rfv '+phot_dir+cat_name[2*i+1].split('/')[1].split('.')[0]+'.cut.coo')
            os.system('rm -rfv '+phot_dir+cat_name[2*i+1].split('/')[1].split('.')[0]+'.cut.reg')

            poi2 = (1*mag_cnd2 + 1*err_cnd2 + 1*cidx_cnd2 == 3)
            print('\n # of point sources (SCI-2) : %d \n' %(np.sum(poi2)))
           
            print('... Writing '+cat_name[2*i+1].split('/')[1].split('.')[0]+'.cut.dat ...')
            f2 = open(phot_dir+cat_name[2*i+1].split('/')[1].split('.')[0]+'.cut.dat','w')
            for k in np.arange(np.sum(poi2)):
                f2.write('%.2f   %.2f   %.2f \n' \
                         %(dat_cat2['x'][poi2][k], dat_cat2['y'][poi2][k], flx2[poi2][k]))            
            f2.close()

            print('... Writing '+cat_name[2*i+1].split('/')[1].split('.')[0]+'.cut.coo ...')
            g2 = open(phot_dir+cat_name[2*i+1].split('/')[1].split('.')[0]+'.cut.coo','w')
            for k in np.arange(np.sum(poi2)):
                g2.write('%.6f   %.6f \n' \
                         %(dat_cat2['ra'][poi2][k], dat_cat2['dec'][poi2][k]))            
            g2.close()

            print('... Writing '+cat_name[2*i+1].split('/')[1].split('.')[0]+'.cut.reg ...')
            h2 = open(phot_dir+cat_name[2*i+1].split('/')[1].split('.')[0]+'.cut.reg','w')
            h2.write('global color=green font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source width=2 \n')
            for k in np.arange(np.sum(poi2)):
                h2.write('image;circle(%.2f, %.2f, 20) \n' \
                         %(dat_cat2['x'][poi2][k], dat_cat2['y'][poi2][k]))            
            h2.close()


    # ====================================== #
    # =============== Step 3 =============== #
    # ====================================== #
    for i in np.arange(len(cat_name)-2):

        # Reading the point source coordinate (WCS)
        ### Note that the WCS coordinate from SExtractor is the one before the distortion correction.
        dat_name = cat_name[i].split('/')[1].split('.')[0]+'.cut.dat'

        dcut1 = np.genfromtxt(phot_dir+cat_name[i].split('/')[1].split('.')[0]+'.cut.coo',
                              dtype=None, encoding='ascii', names=('ra','dec'))
        dcut2 = np.genfromtxt(phot_dir+cat_name[i+2].split('/')[1].split('.')[0]+'.cut.coo',
                              dtype=None, encoding='ascii', names=('ra','dec'))

        ra1, dec1 = dcut1['ra'], dcut1['dec']
        ra2, dec2 = dcut2['ra'], dcut2['dec']

        dcut1 = np.genfromtxt(phot_dir+cat_name[i].split('/')[1].split('.')[0]+'.cut.dat',
                              dtype=None, encoding='ascii', names=('x','y','flx'))
        dcut2 = np.genfromtxt(phot_dir+cat_name[i+2].split('/')[1].split('.')[0]+'.cut.dat',
                              dtype=None, encoding='ascii', names=('x','y','flx'))


        # Matching point sources
        tol = 0.2/3600.0

        src1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
        src2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)

        idx, d2d, d3d = src1.match_to_catalog_sky(src2)
        matched = d2d.value < tol
        midx1 = np.where(matched)[0]
        midx2 = idx[matched]

        print("# Point sources --- {0:d},{1:d} : {2:d} matched \n".format(len(ra1), len(ra2), np.sum(matched)))


        # Writing matched catalogs & regions
        f1 = open(phot_dir+cat_name[i].split('/')[1].split('.')[0]+'.mat.coo','w')
        for k in np.arange(len(midx1)):
            f1.write('%.2f  %.2f  %.2f \n' %(dcut1['x'][midx1][k], dcut1['y'][midx1][k], dcut1['flx'][midx1][k]))
        f1.close()

        g1 = open(phot_dir+cat_name[i].split('/')[1].split('.')[0]+'.mat.reg','w')
        g1.write('global color=green font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source width=2 \n')
        for k in np.arange(len(midx1)):
            g1.write('image;circle(%.2f, %.2f, 20 \n' \
                     %(dcut1['x'][midx1][k], dcut1['y'][midx1][k]))
        g1.close()

        f2 = open(phot_dir+cat_name[i+2].split('/')[1].split('.')[0]+'.mat.coo','w')
        for k in np.arange(len(midx2)):
            f2.write('%.2f  %.2f  %.2f \n' %(dcut2['x'][midx2][k], dcut2['y'][midx2][k], dcut2['flx'][midx2][k]))
        f2.close()

        g2 = open(phot_dir+cat_name[i+2].split('/')[1].split('.')[0]+'.mat.reg','w')
        g2.write('global color=green font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source width=2 \n')
        for k in np.arange(len(midx2)):
            g2.write('image;circle(%.2f, %.2f, 20 \n' \
                     %(dcut2['x'][midx2][k], dcut2['y'][midx2][k]))
        g2.close()


# Printing the running time
print('\n')
print('--- %s seconds ---' %(time.time()-start_time))
