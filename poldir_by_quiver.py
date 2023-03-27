#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 10:59:59 2023

@author: nbadolo
"""

"""
Determination de la direction de la polaristion des données de radmc3d à l'aide de
quiver en vue de la comparer à celle de PolDir de radmc3dPy. Travail préliminaire 
à celui à envoyer à Miguel. C'est bon. Il marche maintenant. Je dois juste comprendre 
la difference entre ce que radmc3dPy et notre routine à ce propos. Au pire, je pourrais
appliquer la methode de radmc3dPy à mes données sphere
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
# Parameters
nDim=1024
nSubDim = 100
pix2mas = 1.9  #en mas/pix

# mas2au = 
x_min = -pix2mas*nSubDim//2
x_max = pix2mas*(nSubDim//2-1)
y_min = -pix2mas*nSubDim//2
y_max = pix2mas*(nSubDim//2-1)

size = (nSubDim, nSubDim)
X, Y= np.meshgrid(np.linspace(-nSubDim/2, nSubDim/2-1, nSubDim), np.linspace(-nSubDim/2,nSubDim/2-1,nSubDim))
X_, Y_= np.meshgrid(np.linspace(-nDim/2,nDim/2-1,nDim), np.linspace(-nDim/2,nDim/2-1,nDim))

X *= pix2mas
Y *= pix2mas
X_ *= pix2mas
Y_ *= pix2mas

X_step = 12
X_step_ = 50

# radmc3d file paths

file_I = '/home/nbadolo/SIM_CODES/RADMC3D/newest_version/radmc3d-2.0/AymardModels/PolData/I_data.fits'
file_Q = '/home/nbadolo/SIM_CODES/RADMC3D/newest_version/radmc3d-2.0/AymardModels/PolData/Q_data.fits'
file_U = '/home/nbadolo/SIM_CODES/RADMC3D/newest_version/radmc3d-2.0/AymardModels/PolData/U_data.fits'

file_lst = [file_I, file_Q, file_U]
nFrames = len(file_lst)

# lists

i_v_arr = np.empty((nFrames,nDim,nDim))
sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))

position = (nDim//2, nDim//2)

# opening Q and U from radmc3d data
for i in range (nFrames):
    hdu = fits.open(file_lst[i])[0]   
    Data = hdu.data
    i_v = Data
    i_v_arr[i] = i_v
                       
    cutout = Cutout2D(i_v, position = position, size = size)
    zoom_hdu = hdu.copy()
    sub_v = cutout.data
    sub_v_arr[i] = sub_v
    
#stop : probleme de dimension. Je dois russir à rogner les cartes  à l'aide cutout    

# les paramètres de stokes
I =  sub_v_arr[0]
Qpol = sub_v_arr[1]
Upol = sub_v_arr[2]
#AOLP = 0.5*np.arctan((Upol/Qpol))
AOLP = 0.5*np.arctan2(Upol,Qpol)
DOLP = np.sqrt((Qpol**2 + Upol**2))/I

# les conventions et notations proposées par radmc3dPy
qr = Qpol/I
ur = Upol/I
lpol = DOLP.clip(1e-60)
qqr = qr/lpol
uur = ur/lpol

# determination des vecteurs de pol par notre appoche
#AOLP = -(AOLP + np.pi/2)
# ii = (uur < 0)
# if True in ii:
#     AOLP[ii] = np.pi - AOLP[ii]
U = np.cos(-(AOLP + np.pi/2))
V = np.sin(-(AOLP + np.pi/2))
# ii = (lpol < 1e-6)
# U = 0.001
# V = 0.001
#%%
#stop : attention, DS9 retourne les image de 90 degre. terminer les plot des orienations 
#et les ajouter aux planches pour les envoyer à Eric pour Julien
#%%
# determination des vecteurs de pol par radmc3dPy
nx = ny = 20
qr = Qpol/I
ur = Upol/I
lpol = DOLP.clip(1e-60)
qqr = qr/lpol
uur = ur/lpol
ang = np.arccos(qqr)/2. # equivaut à 0.5*arccos(Q/PIL) ce qui semble équival à AOLP = 0.5.sqrt(U/Q)
ang = -(ang + np.pi/2)  # Convention d'angle utiliser dans radmc3dPy. A comprendre !!!!
ii = (uur < 0)
if True in ii:
    ang[ii] = np.pi - ang[ii]
vx = np.cos(ang)
vy = np.sin(ang)

ii = (lpol < 1e-6)
vx[ii] = 0.001
vy[ii] = 0.001
#%%
plt.figure('Pol vect by my code')
plt.clf()
plt.subplot(1,2,1)# Ce que mon code fait
plt.imshow(np.sqrt(Qpol**2+Upol**2), cmap = 'inferno', origin='lower', vmin=np.min(np.sqrt(Qpol**2+Upol**2)), 
                vmax=np.max(np.sqrt(Qpol**2+Upol**2)), extent = [x_min , x_max, y_min , y_max])
plt.xlabel('mas', size=14)
plt.ylabel('mas', size=14)
plt.colorbar(label='ADU', shrink = 0.6)
q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U[::X_step,::X_step], V[::X_step,::X_step], color='w', pivot='mid', scale=2. * np.max([nx, ny]), headwidth=1e-10,
          headlength=1e-10, headaxislength=1e-10)

#plt.quiverkey(q, X = 0.9, Y = 1.03, U = 0.06, label='Pol Dir ('+ str(0.06) +' without units)', labelpos='W', color = 'r')
plt.title('Lin. Pol. and vect. with psi =-[0.5*arctan2(u,q) + pi/2]', size = 14)

plt.subplot(1,2,2) # ce que radmc3dPy fait
plt.imshow(np.sqrt(Qpol**2 + Upol**2), cmap = 'inferno', origin='lower', vmin=np.min(np.sqrt(Qpol**2 + Upol**2)), 
                vmax=np.max(np.sqrt(Qpol**2 + Upol**2)), extent = [x_min , x_max, y_min , y_max])
plt.xlabel('mas', size=14)
plt.ylabel('mas', size=14)
plt.colorbar(label='ADU', shrink = 0.6)
q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],vx[::X_step,::X_step], vy[::X_step,::X_step], color='w', pivot='mid', scale=2. * np.max([nx, ny]), headwidth=1e-10,
            headlength=1e-10, headaxislength=1e-10)
#plt.quiverkey(q, X = 0.9, Y = 1.03, U = 0.06, label='Pol Dir ('+ str(0.06) +' without units)', labelpos='W', color = 'r')
plt.title('Lin. Pol. and vect. by radmc3dPy', size = 14)
plt.savefig('/home/nbadolo/SIM_CODES/RADMC3D/newest_version/radmc3d-2.0/AymardModels/PolData/pythonPolDir.pdf', dpi=100, bbox_inches ='tight')
plt.savefig('/home/nbadolo/SIM_CODES/RADMC3D/newest_version/radmc3d-2.0/AymardModels/PolData/pythonPolDir.png', dpi=100, bbox_inches ='tight')


