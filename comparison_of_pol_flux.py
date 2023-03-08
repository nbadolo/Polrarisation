#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:40:42 2022

@author: nbadolo
"""

"""
#
# Pour la comparaison entre la polarisation linéaire standart  et la polarisation azimutale de Engler
# 
"""


import numpy as np
import os 
from os.path import exists
from astropy.io import fits
from scipy import optimize
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.pyplot import Figure, subplot
    
#%% 
star_name = 'R_Aqr'
  
fdir= '/home/nbadolo/Bureau/Aymard/Donnees_sph/log/' + star_name +'/'
   #fdir_star = fdir +'star/alone/'
   
#fdir_star_ = fdir+'star/both/V_N_R/'
fdir_star_ = fdir+'star/alone/V/'
fdir_psf = fdir+'psf'
fname1 ='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2 ='-zpl_science_p23_REDUCED'
file_I_star = fdir_star_ + fname1+'_I'+ fname2 +'_I.fits'
file_PI_star = fdir_star_+ fname1+'_PI'+ fname2 +'_PI.fits'
file_DOLP_star = fdir_star_ + fname1 +'_DOLP' + fname2 +'_DOLP.fits'
file_AOLP_star = fdir_star_ + fname1 +'_AOLP'+ fname2 +'_AOLP.fits' 
file_Q_PHI_star =  fdir_star_ + fname1 +'_Q_PHI'+ fname2 +'_Q_PHI.fits'
file_I_psf = fdir_psf+ fname1 + '_I'+ fname2 +'_I.fits'
file_PI_psf = fdir_psf+fname1 + '_PI'+ fname2 +'_PI.fits'
file_DOLP_psf = fdir_psf+ fname1 + '_DOLP'+ fname2 +'_DOLP.fits'
file_AOLP_psf = fdir_psf + fname1 + '_AOLP' + fname2 + '_AOLP.fits'
 
file_lst = [file_PI_star,file_Q_PHI_star]
lst_Frame_name = ['standar polarized flux, based on Q and U', 'Q_phi in log scale', 'Engler polrised flux,  based on Q_phi']
nFrames = len(file_lst)


nDim = 1024
nSubDim = 200 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
nDimfigj = [9,10,11]
nDimfigk = [0,1,2]
vmin0 = 3.5
vmax0 = 15
pix2mas = 6.8  #en mas/pix
x_min = -pix2mas*nSubDim//2
x_max = pix2mas*(nSubDim//2-1)
y_min = -pix2mas*nSubDim//2
y_max = pix2mas*(nSubDim//2-1)


sub_v_arr =np.empty((nFrames,nSubDim,nSubDim))
im_name_lst = ['Mira I','Mira PI','Mira DOLP','Mira AOLP',
'HD204971 I','HD204971 PI','HD204971 DOLP','HD204971 AOLP']
Vmin=np.empty((nFrames))
Vmax=np.empty((nFrames))
Vmin2=np.empty((nFrames))
Vmax2=np.empty((nFrames))

position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)

x, y = np.meshgrid(np.arange(nSubDim),np.arange(nSubDim)) #cree un tableau 

R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1)

r_mas=pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde

# """
# Filtre utilisé: I_PRIM 
# """
for i in range(nFrames):
   hdu = fits.open(file_lst[i])   
   data = hdu[0].data   
   i_v = data[0,:,:]
  
   cutout = Cutout2D(i_v, position=position, size=size)
   zoom_hdu = hdu.copy()
   sub_v = cutout.data
   sub_v_arr[i] = sub_v
  
   Vmin[i] = np.min(np.log10(sub_v + np.abs(np.min(sub_v))+10))
   Vmax [i] = np.max(np.log10(sub_v + np.abs(np.min(sub_v))+10))
   
   Vmin2[i]= np.min(np.log10(np.abs(sub_v) + 10))
   Vmax2[i]= np.min(np.log10(np.abs(sub_v) + 10))
   
  
plt.figure('flux polarized comparison')
plt.clf()
for i in range(nFrames):
    plt.subplot(1, 2, i+1)
    plt.imshow(np.log10(sub_v_arr[i] + np.abs(np.min(sub_v_arr[i])) + 10), cmap='inferno', origin='lower',vmin=Vmin[i], vmax=Vmax[i], extent = [x_min , x_max, y_min , y_max])
    plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)
    plt.title(star_name + ' ' + f'{lst_Frame_name[i]}', fontsize=10)
    plt.xlabel('Relative R.A.(mas)', size=10) 
    plt.ylabel('Relative Dec.(mas)', size=10)
plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/Comp/Q_phi.pdf', 
            dpi=100, bbox_inches ='tight')
plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/Comp/Q_phi.png', 
            dpi=100, bbox_inches ='tight')

#%%
plt.figure('flux polarized comparison_')
plt.clf()
plt.subplot(1, 2, 1)
plt.imshow(np.log10(sub_v_arr[0] + np.abs(np.min(sub_v_arr[0])) + 10), cmap='inferno', origin='lower',vmin=Vmin[0], vmax=Vmax[0], extent = [x_min , x_max, y_min , y_max])
plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)
plt.title(star_name + ' ' + f'{lst_Frame_name[0]}', fontsize=10)
plt.xlabel('Relative R.A.(mas)', size=10) 
plt.ylabel('Relative Dec.(mas)', size=10)
plt.subplot(1, 2, 2)
plt.imshow(np.log10(np.abs(sub_v_arr[1]) +  10), cmap='inferno', origin='lower',vmin = 1.25, vmax = 1.15, extent = [x_min , x_max, y_min , y_max])
plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)
plt.title(star_name + ' ' + f'{lst_Frame_name[2]}', fontsize=10)
plt.xlabel('Relative R.A.(mas)', size=10) 
plt.ylabel('Relative Dec.(mas)', size=10)
plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/Comp/pol_flux.pdf', 
            dpi=100, bbox_inches ='tight')
plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/Comp/pol_flux.png', 
            dpi=100, bbox_inches ='tight')