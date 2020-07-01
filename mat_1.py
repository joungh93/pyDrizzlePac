#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:49:33 2019

@author: jlee
"""


import numpy as np
import glob, os
import time
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy import units as u
import init_param as ip


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


# ====================================== #
# =============== Step 1 =============== #
# ====================================== #
for w_dir in glob.glob('Phot_*'):

    w_dir = w_dir+'/'
   
    print(' \n')
    print('#---------- Working on '+w_dir+' ----------#')
    
    os.system('rm -rfv '+w_dir+'mat_1.log')
    l = open(w_dir+'mat_1.log','w')
    l.write('Magnitude error cut : %.2f \n' %(ip.err_cut))
    l.write('Lower concentration index cut : %.2f \n' %(ip.cidx_lo_cut))
    l.write('Higher concentration index cut : %.2f \n' %(ip.cidx_hi_cut))
    l.write('Brighter magnitude cut : %.2f \n' %(ip.mag_lo_cut))
    l.write('Fainter magnitude cut : %.2f \n' %(ip.mag_hi_cut))
    l.close()
    
    # ----- Reading SEP photometry results ----- #
    phot_dir = w_dir

    imglist = np.genfromtxt(phot_dir+'image_names.log', dtype=None,
                            encoding='ascii', names=('raw','sci'))
    nshot = len(np.unique(imglist['raw']))
    nchip = len(imglist) // nshot

    pkl_name = glob.glob(phot_dir+'*.pkl')
    pkl_name = sorted(pkl_name)
            
    # ====================================== #
    # =============== Step 1 =============== #
    # ====================================== #
    os.system('rm -rfv '+phot_dir+'mat_1.pdf')

    with PdfPages(phot_dir+'mat_1.pdf') as pdf:
        for i in np.arange(nshot):

            # ---- Figure setting ---- #
            fig = plt.figure(1, figsize=(10,10))
            fig.suptitle(pkl_name[nchip*i].split('/')[1].split('.')[0][:-1], x=0.5, y=0.95, 
                         fontsize=20.0, fontweight='bold', ha='center', va='center')
            gs = GridSpec(2, 2, left=0.10, bottom=0.10, right=0.92, top=0.92,
                          width_ratios=[1., 1.], height_ratios=[1., 1.], wspace=0.05, hspace=0.05)

            for j in np.arange(nchip):
                # ----- Reading data ----- #
                print('... Reading '+pkl_name[nchip*i+j].split('/')[1]+' ...')
                sci_name = pkl_name[nchip*i+j].split('/')[1].split('.pkl')[0]

                df = pd.read_pickle(pkl_name[nchip*i+j])

                mag_cnd = ((df['mag2'].values > ip.mag_lo_cut) & \
                           (df['mag2'].values < ip.mag_hi_cut))

                # Figure 1 : mag-merr diagram
                ax = fig.add_subplot(gs[0,j])
                xt, xl = [15.0, 20.0, 25.0, 30.0], [14.0, 29.0]
                yt, yl = [0.0, 0.2, 0.4, 0.6], [0.0, 0.5]
                ax.set_xticks(xt)
                ax.set_xticklabels(xt, fontsize=15.0)
                ax.set_yticks(yt)
                ax.set_yticklabels(yt, fontsize=15.0)
                ax.tick_params(labelbottom=False)
                if (j == 0):
                    ax.set_ylabel('magnitude error', fontsize=15.0)
                else:
                    ax.tick_params(labelleft=False)
                ax.set_xlim(xl)
                ax.set_ylim(yl)
                ax.tick_params(width=1.5, length=8.0)
                ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
                ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=4))
                ax.tick_params(width=1.5, length=5.0, which='minor')
                for axis in ['top', 'bottom', 'left', 'right']:
                    ax.spines[axis].set_linewidth(1.5)
                # -------------------- #
                xx = ip.mag_lo_cut + (ip.mag_hi_cut-ip.mag_lo_cut)*np.arange(1000)/1000.
                ax.fill_between(xx, np.repeat(0.0, 1000), np.repeat(ip.err_cut, 1000),
                                color='khaki', alpha=0.6)
                X = df['mag2'].values
                Y = df['e_mag1'].values
                ax.plot(X, Y, 'o', ms=3.5, color='gray', mec='black', mew=0.7, alpha=0.7)
                ax.text(0.05, 0.95, 'SCI-{0:d}'.format(j+1), fontsize=17.5, fontweight='bold',
                        ha='left', va='top', transform=ax.transAxes)

                err_cnd = ((Y > 0.0) & (Y < ip.err_cut))

                # Figure 2 : mag-concentration diagram
                ax = fig.add_subplot(gs[1,j])
                xt, xl = [15.0, 20.0, 25.0, 30.0], [14.0, 29.0]
                yt, yl = [-1.0, -0.5, 0.0, 0.5, 1.0], [-1.0, 1.0]
                ax.set_xticks(xt)
                ax.set_xticklabels(xt, fontsize=15.0)
                ax.set_yticks(yt)
                ax.set_yticklabels(yt, fontsize=15.0)
                ax.set_xlabel('magniutude', fontsize=15.0)
                if (j == 0):
                    ax.set_ylabel('C index', fontsize=15.0)
                else:
                    ax.tick_params(labelleft=False)
                ax.set_xlim(xl)
                ax.set_ylim(yl)
                ax.tick_params(width=1.5, length=8.0)
                ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
                ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
                ax.tick_params(width=1.5, length=5.0, which='minor')
                for axis in ['top', 'bottom', 'left', 'right']:
                    ax.spines[axis].set_linewidth(1.5)
                # -------------------- #
                xx = ip.mag_lo_cut + (ip.mag_hi_cut-ip.mag_lo_cut)*np.arange(1000)/1000.
                ax.fill_between(xx, np.repeat(ip.cidx_lo_cut, 1000), np.repeat(ip.cidx_hi_cut, 1000),
                                color='khaki', alpha=0.6)
                X = df['mag2'].values
                Y = df['mag1'].values - df['mag2'].values
                ax.plot(X, Y, 'o', ms=3.5, color='gray', mec='black', mew=0.7, alpha=0.7)

                cidx_cnd = ((Y > ip.cidx_lo_cut) & (Y < ip.cidx_hi_cut))

    # ====================================== #
    # =============== Step 2 =============== #
    # ====================================== #

                # Point source selection
                poi = (1*mag_cnd + 1*err_cnd + 1*cidx_cnd == 3)
                df_poi = df[poi]
                n_poi = np.sum(poi)
                print('\n # of point sources (SCI-{0:d}) : {1:d} \n'.format(j+1, n_poi))

                f = open(phot_dir+sci_name+'.cut.dat', 'w')
                g = open(phot_dir+sci_name+'.cut.coo', 'w')
                h = open(phot_dir+sci_name+'.cut.reg', 'w')
                for k in np.arange(n_poi):
                    f.write('%.2f   %.2f   %.2f \n' \
                            %(df_poi['x'].values[k], df_poi['y'].values[k], df_poi['flx2'].values[k]))
                    g.write('%.6f   %.6f \n' \
                            %(df_poi['ra'].values[k], df_poi['dec'].values[k]))
                    h.write('image;circle(%.2f, %.2f, 20) \n' \
                            %(df_poi['x'].values[k], df_poi['y'].values[k]))
                f.close()
                g.close()
                h.close()

            pdf.savefig()
            plt.close()

    # os.system('evince '+phot_dir+'mat_1.pdf &')

    # ====================================== #
    # =============== Step 3 =============== #
    # ====================================== #

    for i in np.arange(nchip*(nshot-1)):
        sci1_name = pkl_name[i].split('/')[1].split('.pkl')[0]
        sci2_name = pkl_name[i+nchip].split('/')[1].split('.pkl')[0]

        # Reading the point source coordinate (WCS)
        dcut1 = np.genfromtxt(phot_dir+sci1_name+'.cut.coo',
                              dtype=None, encoding='ascii', names=('ra','dec'))
        dcut2 = np.genfromtxt(phot_dir+sci2_name+'.cut.coo',
                              dtype=None, encoding='ascii', names=('ra','dec'))
        ra1, dec1 = dcut1['ra'], dcut1['dec']
        ra2, dec2 = dcut2['ra'], dcut2['dec']

        # Reading the point source data
        dcut1 = np.genfromtxt(phot_dir+sci1_name+'.cut.dat',
                              dtype=None, encoding='ascii', names=('x','y','flx'))
        dcut2 = np.genfromtxt(phot_dir+sci2_name+'.cut.dat',
                              dtype=None, encoding='ascii', names=('x','y','flx'))

        # Matching point sources
        tol = ip.tolerance / 3600.0

        src1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
        src2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)

        idx, d2d, d3d = src1.match_to_catalog_sky(src2)
        matched = d2d.value < tol
        n_mch = np.sum(matched)
        midx1 = np.where(matched)[0]
        midx2 = idx[matched]

        print("# Point sources --- {0:d},{1:d} : {2:d} matched \n".format(len(ra1), len(ra2), n_mch))

        # Writing matched catalogs & regions
        f1 = open(phot_dir+sci1_name+'.mat.coo','w')
        f2 = open(phot_dir+sci2_name+'.mat.coo','w')
        g1 = open(phot_dir+sci1_name+'.mat.reg','w')
        g1.write('global color=green font="helvetica 10 normal" ')
        g1.write('select=1 edit=1 move=1 delete=1 include=1 fixed=0 source width=2 \n')
        g2 = open(phot_dir+sci2_name+'.mat.reg','w')
        g2.write('global color=green font="helvetica 10 normal" ')
        g2.write('select=1 edit=1 move=1 delete=1 include=1 fixed=0 source width=2 \n')
        for k in np.arange(n_mch):
            f1.write('%.2f  %.2f  %.2f \n' \
                     %(dcut1['x'][midx1][k], dcut1['y'][midx1][k], dcut1['flx'][midx1][k]))
            f2.write('%.2f  %.2f  %.2f \n' \
                     %(dcut2['x'][midx2][k], dcut2['y'][midx2][k], dcut2['flx'][midx2][k]))
            g1.write('image;circle(%.2f, %.2f, 20 \n' \
                     %(dcut1['x'][midx1][k], dcut1['y'][midx1][k]))
            g2.write('image;circle(%.2f, %.2f, 20 \n' \
                     %(dcut2['x'][midx2][k], dcut2['y'][midx2][k]))
        f1.close()
        f2.close()
        g1.close()
        g2.close()


# Printing the running time
print('\n')
print('--- %s seconds ---' %(time.time()-start_time))
