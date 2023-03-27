#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 12:16:50 2020

@author: emily
"""

import numpy as np
from astropy.io import fits
from scipy import optimize

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def sphereprep(afile, new_FOV, outfits, nogauss = False):

    """
    A function to take a cut out of SPHERE data centred on the star

    Parameters
    ----------

    **afile** : Fits file of observations

    **new_FOV** : FOV in mas

    **outfits** Name of the output FITS file

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



    #Calculate fov in pixels
    radius_px_dec =  new_FOV / (hdr['PIXSCAL'])  #should be divided by 2 but pixscal seems to be twice as large
    radius_px = round(radius_px_dec)
    limx = - radius_px*hdr['PIXSCAL']/2


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

    #Export to FITS
    hdr['EXTNAME'] = 'INTENSITY'
    hdu_int = fits.PrimaryHDU(data=intens_clip, header=hdr)
    hdr_polflux = fits.Header()
    hdr_polflux['EXTNAME'] = 'POLFLUX'
    hdu_polflux = fits.PrimaryHDU(data=polflux_clip, header=hdr_polflux)
    hdr_poldeg = fits.Header()
    hdr_poldeg['EXTNAME'] = 'POLDEG'
    hdu_poldeg = fits.PrimaryHDU(data=poldeg_clip, header=hdr_poldeg)
    hdr_polangle = fits.Header()
    hdr_polangle['EXTNAME'] = 'POLANGLE'
    hdu_polangle = fits.PrimaryHDU(data=polangle_clip, header=hdr_polangle)
    hdul = fits.HDUList([hdu_int, hdu_polflux, hdu_poldeg, hdu_polangle])
    hdul.writeto(outfits, overwrite=True)
    print('\nData written in '+outfits+'\n')


    return intens_clip, polflux_clip, poldeg_clip, limx, hdr, polangle_clip, radius_px

