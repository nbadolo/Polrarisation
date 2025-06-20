#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 12:52:29 2023

@author: nbadolo
"""


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
la difference entre ce que fait radmc3dPy et notre routine à ce propos. Au pire, je pourrais
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

star_name = '17_Lep'
obs_mode = 'alone'
filtr = 'I_PRIM'

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
fname1 = '/home/nbadolo/Bureau/Aymard/Donnees_sph/large_log/'
fname2= 'zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED_'
fname3 = '-zpl_science_p23_REDUCED_'
file_Q = fname1+star_name+'/star/'+ obs_mode+ '/'+ filtr + '/' + fname2 + 'Q' + fname3 + 'Q' + '.fits'
file_U =  fname1+star_name+'/star/'+ obs_mode+ '/' + filtr +'/' + fname2 + 'U' + fname3 + 'U' + '.fits'
file_AOLP = fname1+star_name+'/star/'+ obs_mode+ '/'+ filtr +'/'+  fname2 + 'AOLP' + fname3 + 'AOLP' + '.fits'

file_lst = [file_Q, file_U, file_AOLP]
nFrames = len(file_lst)

# lists

i_v_arr = np.empty((nFrames,nDim,nDim))
sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))

position = (nDim//2, nDim//2)

# opening Q and U from radmc3d data
for i in range (nFrames):
    hdu = fits.open(file_lst[i])[0]   
    Data = hdu.data
    i_v = Data[0,:,:]
    i_v_arr[i] = i_v
                       
    cutout = Cutout2D(i_v, position = position, size = size)
    zoom_hdu = hdu.copy()
    sub_v = cutout.data
    sub_v_arr[i] = sub_v

Q =  sub_v_arr[0]  
U = sub_v_arr[1]
AOLP_1 = sub_v_arr[2]
AOLP_2 = 0.5*np.arctan(U/Q)
AOLP_3 = 0.5*np.arctan2(U,Q)

plt.figure(0)
plt.clf()
plt.subplot(2,2,1)
plt.imshow(AOLP_1, cmap='inferno', origin='lower', extent = [x_min , x_max, y_min , y_max])
#plt.text(size[0]//10, 2*pix2mas*size[1]//6.,'AOLP_1', color='w',fontsize='large', ha='center')
plt.title('AOLP_1')
plt.subplot(2,2,2)
plt.imshow(AOLP_2, cmap='inferno', origin='lower', extent = [x_min , x_max, y_min , y_max])
#plt.text(size[0]//10, 2*pix2mas*size[1]//6.,'AOLP_2', color='w',fontsize='large', ha='center')
plt.title('AOLP_2')
plt.subplot(2,2,3)
plt.imshow(AOLP_3, cmap='inferno', origin='lower', extent = [x_min , x_max, y_min , y_max])
plt.text(size[0]//10, 2*pix2mas*size[1]//6.,'AOLP_3', color='w',fontsize='large', ha='center')
#plt.title('AOLP_3')