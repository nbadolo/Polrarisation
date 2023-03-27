#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 14:02:24 2020

@author: emily
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from SPHERE_binning import cutsphere
from scipy.optimize import curve_fit
from gauss_2d import gaussian2d
import matplotlib.ticker as ticker
from scipycb import *

def angle_array(shape,centerx=None,centery=None,verbose=True,fullOutput=False):
    """
    Creates a 2d array with the angle in rad (0 is North, then the array is positive going East
    and becomes negative after we pass South).
    Input: 
        - shape: a tuple indicating the desired shape of the output array, e.g. (100,100)
                The 1st dim refers to the y dimension and the 2nd to the x dimension
        - centerx: the center of the frame from which to compute the distance from
                    by default shape[1]/2 (integer division). Accepts numerical value
        - centery: same for the y dimension
        - verbose: to print a warning for even dimensions
        - fullOutput: if True returns the angle array and in addition the 2d array
                    of x values and y values in 2nd and 3rd ouptuts.
    """
    if len(shape) != 2 :
        raise ValueError('The shape must be a tuple of 2 elements for the y and x dimension!')
    if centerx == None:
        centerx = shape[1]//2
        if np.mod(shape[1],2) == 0 and verbose:
            print('The X dimension is even ({0:d}), the center is assumed to be in {1:d}'.format(shape[1],centerx))
    if centery == None:
        centery = shape[0]//2
        if np.mod(shape[0],2) == 0 and verbose:
            print('The Y dimension is even ({0:d}), the center is assumed to be in {1:d}'.format(shape[0],centery))
    x_array = np.arange(shape[1])-centerx
    y_array = np.arange(shape[0])-centery
    xx_array,yy_array=np.meshgrid(x_array,y_array)
    theta = -np.arctan2(xx_array,yy_array)
    if fullOutput:
        return theta,xx_array,yy_array
    return theta

def sphereprep(afile,  radius_px, nogauss = False):
    
    """
    A function to take a cut out of SPHERE data centred on the star
    
    Parameters
    ----------

    **afile** : Fits file of observations

    **new_FOV** : FOV in mas

    Returns
    ----------
    Clips of observations and limx for plotting
    """
    # open file
    with fits.open(afile) as hdu:
        data = hdu[0].data
        I1 = data[0,:,:]
        I2 = data[2,:,:]
        Q = data[1,:,:]
        U = data[3,:,:]
        hdr = hdu[0].header
    intens = (I1**2 + I2**2)**0.5
    polflux = (Q**2 + U**2)**0.5
    poldeg = ((Q**2 + U**2)**0.5)/intens
    polangle = np.arctan2(U, Q)*(180/np.pi)*(0.5)
    plt.figure()
    plt.imshow(poldeg, norm = PowerNorm(1, vmax =0.1))


    #Calculate fov in pixels
#    radius_px_dec =  new_FOV / (hdr['PIXSCAL'])  #should be divided by 2 but pixscal seems to be twice as large
#    radius_px = round(radius_px_dec)
    limx = - radius_px*hdr['PIXSCAL']/2 ### Miguel ou should change this!! this was ue to a header keyword problem


    # Fit gaussian to find centre
    if nogauss == True:
        loc_max = np.unravel_index(intens.argmax(), intens.shape)
        centre_x = loc_max[0]
        centre_y = loc_max[1]
    else:
        params_1 = fitgaussian(intens)
        centre_x = int(round(params_1[1]))
        centre_y = int(round(params_1[2]))  

    # Cut out central clip
    intens_clip = intens[centre_x - radius_px:centre_x + radius_px, centre_y - radius_px:centre_y + radius_px]
    poldeg_clip = poldeg[centre_x - radius_px:centre_x + radius_px, centre_y - radius_px:centre_y + radius_px]
    polangle_clip = polangle[centre_x - radius_px:centre_x + radius_px, centre_y - radius_px:centre_y + radius_px]
    polflux_clip = polflux[centre_x - radius_px:centre_x + radius_px, centre_y - radius_px:centre_y + radius_px]
    Q_clip =  Q[centre_x - radius_px:centre_x + radius_px, centre_y - radius_px:centre_y + radius_px]
    U_clip =  U[centre_x - radius_px:centre_x + radius_px, centre_y - radius_px:centre_y + radius_px]
#    print(hdr)

    return intens_clip, polflux_clip, poldeg_clip, limx, hdr, polangle_clip, radius_px,Q_clip, U_clip
plt.close('all')
BHa_dc = '/lhome/emily/SPHERE_DC_DATA/HIP 107516_B_Ha,CntHa_2018-08-16_zpl_science_p23_142169/zpl_science_p23-ZPL_SCIENCE_P23_REDUCED-sci1.fits'

plt.close('all')


clips = sphereprep(BHa_dc, 300)
I = clips[0]
Q = clips[7]
U =clips[8]
poldeg = clips[2]

phi = angle_array(Q.shape)

Q_phi =- (Q *np.cos(2*phi) + U*np.sin(2*phi))
U_phi = -Q* np.sin(2*phi) + U*np.cos(2*phi)
plt.figure()


plt.subplot(121)
plt.imshow(phi)
plt.subplot(122)
plt.imshow(poldeg, cmap = 'CMRmap', vmax = 0.1)

plt.figure()
plt.subplot(121)
plt.imshow(Q_phi,  cmap='seismic',vmin = -600, vmax = 600)
plt.xlim(200,400)
plt.ylim(200,400)

plt.colorbar()
plt.subplot(122)
plt.imshow(U_phi,  cmap='seismic',  vmin = -600, vmax = 600)
plt.colorbar()
plt.xlim(200,400)
plt.ylim(200,400)

plt.figure()
plt.subplot(121)
plt.imshow(Q_phi/I,  cmap='CMRmap', vmin = 0, vmax = 0.1)
plt.xlim(200,400)
plt.ylim(200,400)
plt.colorbar()
plt.subplot(122)
plt.imshow(U_phi/I,  cmap='CMRmap', vmin = 0, vmax = 0.1)
plt.xlim(200,400)
plt.ylim(200,400)

plt.colorbar()
