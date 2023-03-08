#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:21:29 2022

@author: nbadolo
"""


"""
Code simplifié pour la comparaison des deux methodes de calcul du Qphi et l'évaluation de leur impact 
sur les vecteurs de polaristaion. Resultats à envoyer à Miguel :
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

file_lst = [file_Q_star, file_U_star, file_Q_phi_star, file_DoLP_star, file_AoLP_star]
nFrames = len(file_lst)


nDim = 1024 
nSubDim = 150 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)

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



i_v_arr = np.empty((nFrames,nDim,nDim))   
sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))


#im_name_lst = ['I','PI','DOLP', 'AOLP']
Vmin = np.empty((nFrames))
Vmax = np.empty((nFrames))

position = (nDim//2,nDim//2)


x, y = np.meshgrid(np.arange(nSubDim), np.arange(nSubDim)) #cree un tableau 

R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1) # creation d'un tableau de distance radiale

r_mas = pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde

ang_arr2 = imtool.angle_array((nSubDim,nSubDim),centerx=None,centery=None)  #calcul of the positional angle selon de Boer
ang_arr1 = imtool.angleN_array((nSubDim,nSubDim),centerx=None,centery=None)  #calcul of the positional angle selon Schmid


# opening Q and U
for i in range (nFrames):
    hdu = fits.open(file_lst[i])[0]   
    Data = hdu.data   
    i_v = Data[0,:,:]
    i_v_arr[i] = i_v.data
    
    #fltr = hdu.header.get('HIERARCH ESO INS3 OPTI5 NAME')     
    #print(fltr)            Q = sub_v_arr[0]
    U = sub_v_arr[1]       
    cutout = Cutout2D(i_v, position=position, size=size)
    #zoom_hdu = hdu.copy()
    sub_v = cutout.data
    sub_v_arr[i] = sub_v
    
    


Qpol = sub_v_arr[0]
Upol = sub_v_arr[1]
Q_phi_sphr = sub_v_arr[2]
phi1 = ang_arr1 #  le postional angle selon Schmid, calcul inspiré par le code de Julien
phi2 = ang_arr2 #  le postional angle selon de Boer et calculer par Julien


Co1 = np.cos(2*phi1) 
Si1 = np.sin(2*phi1)
Q_r = Qpol*Co1 + Upol*Si1  # polarisation radiale de Schmid 2006 
U_r = -Qpol*Si1 + Upol*Co1

Co2 = np.cos(2*phi2) 
Si2 = np.sin(2*phi2)
Q_phi = -Qpol*Co2 - Upol*Si2 # polarisation azimutale de De Boer 2020
U_phi = Qpol*Si2 -Upol*Co2

Vmin =[np.min(Q_phi),np.min(Q_phi_sphr), np.min(Q_r)]
Vmax = [np.max(Q_phi),np.max(Q_phi_sphr), np.max(Q_r)] 

U = sub_v_arr[3]*np.cos(np.pi*sub_v_arr[4]/180)
V = sub_v_arr[3]*np.sin(np.pi*sub_v_arr[4]/180)

# U1 = sub_v_arr[2]*np.cos(np.pi*sub_v_arr[3]/180)
# V1 = sub_v_arr[2]*np.sin(np.pi*sub_v_arr[3]/180)

# U2 = sub_v_arr[2]*np.cos(np.pi*sub_v_arr[3]/180)
# V2 = sub_v_arr[2]*np.sin(np.pi*sub_v_arr[3]/180)



#%%
fig = plt.figure('azm_rad_param_' +star_name)
plt.clf()
#fig.set_size_inches(12,6, forward = True)
plt.subplot(2,2,1)
plt.imshow((Q_phi), cmap = 'inferno', origin='lower', vmin=Vmin[0], 
                vmax=Vmax[0], extent = [x_min , x_max, y_min , y_max])
plt.colorbar(label='ADU', shrink = 0.6)
q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U[::X_step,::X_step], V[::X_step,::X_step], color='w', pivot='mid', width=0.004, headwidth=0.8, headlength=1.5)
plt.quiverkey(q, X = 0.9, Y = 1.03, U = 0.02, label='Pol Dir ('+ str(0.02) +' without units)', labelpos='W', color = 'r')
plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
           ' Q_phi,' +' De Boer 2020', color='w',
       fontsize='small', ha='center')
plt.xlabel('mas', size=14)
plt.ylabel('mas', size=14)



plt.subplot(2,2,2)
plt.imshow((Q_r), cmap = 'inferno', origin='lower', vmin=Vmin[2], 
                vmax=Vmax[2], extent = [x_min , x_max, y_min , y_max])
plt.colorbar(label='ADU', shrink = 0.6)
q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U[::X_step,::X_step], V[::X_step,::X_step], color='w', pivot='mid', width=0.004, headwidth=0.8, headlength=1.5)
plt.quiverkey(q, X = 0.9, Y = 1.03, U = 0.02, label='Pol Dir ('+ str(0.02) +' without units)', labelpos='W', color = 'r')
plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
           ' Q_r,' +' Schmid 2006', color='w',
    fontsize='small', ha='center')

plt.xlabel('mas', size=14)
plt.ylabel('mas', size=14)
#plt.title('polarisation radiale de  Schmid 2006', size = 12)

plt.subplot(2,2,3)
plt.imshow((U_phi), cmap = 'inferno', origin='lower', vmin=Vmin[2], 
                vmax=Vmax[2], extent = [x_min , x_max, y_min , y_max])
plt.colorbar(label='ADU', shrink = 0.6)
q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U[::X_step,::X_step], V[::X_step,::X_step], color='w', pivot='mid', width=0.004, headwidth=0.8, headlength=1.5)
plt.quiverkey(q, X = 0.9, Y = 1.03, U = 0.02, label='', color = 'r')
plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
           ' U_phi,' +' De Boer 2020', color='w',
    fontsize='small', ha='center')

plt.xlabel('mas', size=14)
plt.ylabel('mas', size=14)

plt.subplot(2,2,4)
plt.imshow((U_r), cmap = 'inferno', origin='lower', vmin=Vmin[2], 
                vmax=Vmax[2], extent = [x_min , x_max, y_min , y_max])
plt.colorbar(label='ADU', shrink = 0.6)
q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U[::X_step,::X_step], V[::X_step,::X_step], color='w', pivot='mid', width=0.004, headwidth=0.8, headlength=1.5)
plt.quiverkey(q, X = 0.9, Y = 1.03, U = 0.02, label= '', color = 'r')
plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
           ' U_r,' +' Schmid 2006', color='w',
    fontsize='small', ha='center')

plt.xlabel('mas', size=14)
plt.ylabel('mas', size=14)

plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/Donnees_Miguel/azm_rad_param_' +star_name+'.pdf', dpi=100, bbox_inches ='tight')
plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/Donnees_Miguel/azm_rad_param_' +star_name+'.png', dpi=100, bbox_inches ='tight')
plt.tight_layout()


#plt.title('polarisation azimutale de Boer 2020', size = 12)

# plt.subplot(2,2,2)
# plt.imshow((Q_phi_sphr), cmap = 'inferno', origin='lower', vmin=Vmin[1], 
#                 vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])
# plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)
# plt.xlabel('mas', size=14)
# plt.ylabel('mas', size=14)
# plt.title('polarisation azimutale de Sphere', size = 12)
#%%

plt.figure('PIL')
plt.clf()
plt.imshow(np.sqrt(Qpol**2+Upol**2), cmap = 'inferno', origin='lower', vmin=np.min(np.sqrt(Qpol**2+Upol**2)), 
                vmax=np.max(np.sqrt(Qpol**2+Upol**2)), extent = [x_min , x_max, y_min , y_max])
plt.xlabel('mas', size=14)
plt.ylabel('mas', size=14)
plt.colorbar(label='ADU', shrink = 0.6)
plt.title('Linear polarisation', size = 18)
plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/Donnees_Miguel/linear_pol'+star_name+'.pdf', dpi=100, bbox_inches ='tight')

stop 

plt.figure('Q_r')
plt.clf()
plt.imshow((Q_r), cmap = 'inferno', origin='lower')
plt.xlabel('pix', size=20)
plt.ylabel('pix', size=20)
plt.title('plarisation radiale de  Schmid 2006', size = 18)
#%%
           
#%%        
