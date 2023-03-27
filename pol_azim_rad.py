#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 09:38:47 2023

@author: nbadolo
"""

"""
bouts de codes pour les calcules des polarisations azimutales et radiales
"""

#packages
import numpy as np
import os
import scipy 
from os.path import exists
from astropy.io import fits
from scipy import optimize
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.pyplot import Figure, subplot
import webbrowser
import Julien_image_tool as imtool

#%%

# bout de code de plotDir à comparer
x = image.x / nc.au
y = image.y / nc.au
iix = [int(np.floor(i)) for i in np.arange(nx) * float(x.shape[0]) / nx]
iiy = [int(np.floor(i)) for i in np.arange(ny) * float(x.shape[0]) / ny]
xr = x[iix]
yr = y[iiy]
xxr, yyr = np.meshgrid(x, y, indexing='ij')
qqr = (np.squeeze(dum_image.image[:, :, 1, ifreq])
        / np.squeeze(dum_image.image[:, :, 0, ifreq]).clip(1e-60))[np.ix_(iix, iiy)]
uur = (np.squeeze(dum_image.image[:, :, 2, ifreq])
        / np.squeeze(dum_image.image[:, :, 0, ifreq]).clip(1e-60))[np.ix_(iix, iiy)]
lpol = np.sqrt(qqr**2 + uur**2).clip(1e-60)
qqr /= lpol
uur /= lpol
ang = np.arccos(qqr) / 2.0
ii = (uur < 0)
if True in ii:
    ang[ii] = np.pi - ang[ii]
vx = np.cos(ang)
vy = np.sin(ang)
ii = (lpol < 1e-6)
vx[ii] = 0.001
vy[ii] = 0.001
plt.quiver(xxr, yyr, vx, vy, color=color, pivot='mid', scale=2. * np.max([nx, ny]), headwidth=1e-10,
            headlength=1e-10, headaxislength=1e-10)

ang_arr2 = imtool.angle_array((nSubDim,nSubDim),centerx=None,centery=None)  #calcul of the positional angle selon de Boer
ang_arr1 = imtool.angleN_array((nSubDim,nSubDim),centerx=None,centery=None)  #calcul of the positional angle selon Schmid

#%%
# Q_file = fits.getdata(Q_path)
# U_file = fits.getdata(U_path)
# Q = np.reshape(Q,(nSubDim,nSubDim))
# U = np.reshape(U,(nSubDim, nSubDim))

#%%
#données d'observtion de SW_Col

star_name = 'SW_Col'  
obsmod = 'both'     
fdir= '/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+ '/'
fdir_star = fdir + 'star/'+obsmod+ '/V_N_R/' 

#fdir_star_fltr = fdir + fdir_star
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_Q_star = fdir_star + fname1+'_Q'+fname2+'_Q.fits'
file_U_star = fdir_star +fname1+'_U'+fname2+'_U.fits'
file_Q_phi_star = fdir_star +fname1+'_Q_PHI'+fname2+'_Q_PHI.fits'
file_DoLP_star = fdir_star +fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AoLP_star = fdir_star +fname1+'_AOLP'+fname2+'_AOLP.fits'

file_lst_obs = [file_Q_star, file_U_star, file_Q_phi_star, file_DoLP_star, file_AoLP_star]
nFrames_obs = len(file_lst_obs)

nDim = 1024 
nSubDim = 150 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)

dpc = 100 
scale = dpc*np.pi /(36*18 ) # pc/arcs

pix2mas = 3.4  # en mas/pix
x_min = -pix2mas*nSubDim//2
x_max = pix2mas*(nSubDim//2-1)
y_min = -pix2mas*nSubDim//2
y_max = pix2mas*(nSubDim//2-1)
X, Y= np.meshgrid(np.linspace(-nSubDim/2,nSubDim/2-1,nSubDim), np.linspace(-nSubDim/2,nSubDim/2-1,nSubDim))
X_, Y_= np.meshgrid(np.linspace(-nDim/2,nDim/2-1,nDim), np.linspace(-nDim/2,nDim/2-1,nDim))

X *= pix2mas
Y *= pix2mas
X_ *= pix2mas
Y_ *= pix2mas

X_step = 10
X_step_ = 50


i_v_arr_obs = np.empty((nFrames_obs, nDim, nDim))   
sub_v_arr_obs = np.empty((nFrames_obs, nSubDim, nSubDim))


im_name_lst = ['I','PI','DOLP', 'AOLP']
Vmin_obs = np.empty((nFrames_obs))
Vmax_obs = np.empty((nFrames_obs))

position = (nDim//2,nDim//2)


x, y = np.meshgrid(np.arange(nSubDim), np.arange(nSubDim)) #cree un tableau 

R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1) # creation d'un tableau de distance radiale

r_mas = pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde

ang_arr2 = imtool.angle_array((nSubDim,nSubDim),centerx=None,centery=None)  #calcul of the positional angle selon de Boer
ang_arr1 = imtool.angleN_array((nSubDim,nSubDim),centerx=None,centery=None)  #calcul of the positional angle selon Schmid

# opening Q and U of observations
for i in range (nFrames_obs):
    hdu_obs = fits.open(file_lst_obs[i])[0]   
    Data_obs = hdu_obs.data   
    i_v_obs = Data_obs[0,:,:]
    i_v_arr_obs[i] = i_v_obs.data
   
    #fltr = hdu.header.get('HIERARCH ESO INS3 OPTI5 NAME')     
    #print(fltr)            
    #Q = sub_v_arr[0]
    #U = sub_v_arr[1]       
    cutout = Cutout2D(i_v_obs, position=position, size=size)
    #zoom_hdu = hdu.copy()
    sub_v_obs = cutout.data
    sub_v_arr_obs[i] = sub_v_obs

# calcul des angles

AOLP_obs =np.rad2deg((0.5*np.arctan(sub_v_arr_obs[1]/sub_v_arr_obs[0])))

AOLP_obs2 =np.rad2deg(0.5*np.arctan(sub_v_arr_obs[0]/sub_v_arr_obs[1]))

AOLP_obs3 = sub_v_arr_obs[4]
plt.figure('angl2')
plt.clf()
plt.subplot(1,3,1) # AOLP selon Julien selon la litterature
plt.imshow(AOLP_obs, cmap = 'inferno', vmin=np.min(AOLP_obs), vmax= np.max(AOLP_obs))
plt.title('AOLP de SW_Col litt', size = 18)
plt.subplot(1,3,2) # AOLP selon Julien et calculer par moi
plt.imshow(AOLP_obs2, cmap = 'inferno',vmin=np.min(AOLP_obs2), vmax= np.max(AOLP_obs2))
plt.title('AOLP de SW_Col zpl_polarmap', size = 18)
plt.subplot(1,3,3) # AOLP selon Julien et deja calculer par le DC
plt.imshow(AOLP_obs3, cmap = 'inferno',vmin=np.min(AOLP_obs3), vmax= np.max(AOLP_obs3))
plt.title('AOLP de SW_Col zpl_polarmap2', size = 18)
plt.savefig('/home/nbadolo/SIM_CODES/RADMC3D/newest_version/radmc3d-2.0/AymardModels/PolData/SW_ColAOLP_compar.pdf', dpi=100, bbox_inches ='tight')
plt.tight_layout()


phi1 = ang_arr1 #  le postional angle selon Schmid, calcul inspiré par le code de Julien
phi2 = ang_arr2 #  le postional angle selon de Boer et calculer par Julien

# Schmid
Co1 = np.cos(2*phi1) 
Si1 = np.sin(2*phi1)
Q_r = Qpol*Co1 + Upol*Si1  # polarisation radiale de Schmid 2006 
U_r = -Qpol*Si1 + Upol*Co1

# De Boer
Co2 = np.cos(2*phi2) 
Si2 = np.sin(2*phi2)
Q_phi = -Qpol*Co2 - Upol*Si2 # polarisation azimutale de De Boer 2020
U_phi = Qpol*Si2 -Upol*Co2

Vmin = [np.min(Q_r), np.min(Q_phi)]
Vmax = [np.max(Q_r), np.max(Q_phi)]
#%%
fig, ax = plt.subplots() 
plt.figure('Q_r')
plt.clf()
plt.imshow(Q_r, cmap = 'inferno', origin='lower', vmin=Vmin[1], vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])
plt.colorbar(label='ADU', shrink = 0.6)
plt.xlabel('pix', size=20)
plt.ylabel('pix', size=20)
plt.title('plarisation radiale de  Schmid 2006', size = 18)
plt.savefig('/home/nbadolo/SIM_CODES/RADMC3D/newest_version/radmc3d-2.0/AymardModels/PolData/pol_rad_pyr_mg100.pdf', dpi=100, bbox_inches ='tight')
plt.tight_layout()
#%%
plt.figure('angl')
plt.clf()
plt.subplot(1,2,1)
plt.imshow(AOLP, cmap = 'inferno')
plt.title('AOLP according to literature', size = 20)
plt.subplot(1,2,2)
plt.imshow(AOLP2, cmap = 'inferno')
plt.title('AOLP used in the DC', size = 20)
plt.savefig('/home/nbadolo/SIM_CODES/RADMC3D/newest_version/radmc3d-2.0/AymardModels/PolData/radmcAOLP_compar.pdf', dpi=100, bbox_inches ='tight')
plt.tight_layout()