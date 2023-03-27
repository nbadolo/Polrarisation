#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 14:13:21 2017

Library to remove bad RAW files of ZIMPOL (opened AO loop).

@author: miguel
"""

#pylint: disable-msg=E1101,C0103

import os
import pickle
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from astropy.io import fits
from astropy.time import Time
from astropy.time import TimeDelta


def archive_raw_info(input_dir, output_file, window_size):
    """
    Creating the pickle file to gather the required information from the RAW data to analyze.

    Parameters
    ----------
    **input_dir** : Directory containing the RAW ZIMPOL files

    **output_file** : Name of the output .pickle file. It will be written in input_dir

    **window_size** : Size (in pixels) of the squared window, centered on the maximum flux location,
    on which the flux of each frame is derived.
    """

    #List of RAW files
    list_files = glob(os.path.join(input_dir, "SPH*.fits"))

       #For each file
    if len(list_files) == 0:

        print("\nNo RAW ZIMPOL files in the input directory")
        print("\nCheck your path ! Aborting...\n")

        return 1

    else:

        #Initializing dictionaries
        dic_filt_b = {}
        dic_filt_c = {}
        dic_flux_b = {}
        dic_flux_c = {}
        dic_dit = {}
        dic_mjd = {}
        dic_targ = {}

        #Window half size
        whs = int(window_size/2)
        list_files.sort()

        #Foe each RAW file
        for n, cur_file in enumerate(list_files):

            filename = os.path.basename(cur_file)
#            print("File %i of %i: " %(n, len(list_files)), filename)
            all_b = np.array([])
            all_c = all_b.copy()

            #Getting images
            try:

                with(fits.open(cur_file)) as hdu:
                    img_c = hdu['Callas'].data
                    img_b = hdu['Bartoli'].data
                    hdr = hdu[0].header

                #Forgetting BIASES
                if hdr['ESO DPR TYPE'] == 'OBJECT':

                    #Exploring the frames for polarization
                    if len(img_b.shape) == 3:
                        for frames_b, frames_c in zip(img_b, img_c):

                            #Locating maxima
                            max_b = np.unravel_index(frames_b.argmax(), frames_b.shape)
                            max_c = np.unravel_index(frames_c.argmax(), frames_c.shape)

                            #Getting flues
                            flux_b = frames_b[max_b[0]-whs:max_b[0]+whs, max_b[1]-whs:max_b[1]+whs]
                            flux_c = frames_c[max_c[0]-whs:max_c[0]+whs, max_c[1]-whs:max_c[1]+whs]
                            flux_b = flux_b.sum()
                            flux_c = flux_c.sum()

                            #Storing
                            all_b = np.append(all_b, flux_b)
                            all_c = np.append(all_c, flux_c)

                    #For classical imaging
                    elif len(img_b.shape) == 2:
                        #Locating maxima
                        max_b = np.unravel_index(img_b.argmax(), img_b.shape)
                        max_c = np.unravel_index(img_c.argmax(), img_c.shape)

                        #Getting flues
                        flux_b = img_b[max_b[0]-whs:max_b[0]+whs, max_b[1]-whs:max_b[1]+whs]
                        flux_c = img_c[max_c[0]-whs:max_c[0]+whs, max_c[1]-whs:max_c[1]+whs]
                        flux_b = flux_b.sum()
                        flux_c = flux_c.sum()

                        #Storing
                        all_b = np.append(all_b, flux_b)
                        all_c = np.append(all_c, flux_c)

                    else:
                        raise TypeError("Shape of the image not supported.")

                    #Filling dictionnaries for each files
                    dic_filt_c[filename] = hdr['ESO INS3 OPTI5 NAME']
                    dic_filt_b[filename] = hdr['ESO INS3 OPTI6 NAME']
                    dic_mjd[filename] = hdr['MJD-OBS']
                    dic_dit[filename] = hdr['ESO DET DIT1']
                    dic_targ[filename] = hdr['OBJECT']
                    dic_flux_b[filename] = all_b
                    dic_flux_c[filename] = all_c

            except KeyError:
#                    bad_dir = os.path.join(os.path.dirname(cur_file), \
#                    os.path.join(os.path.pardir, 'bad'))
#                    os.rename(cur_file, os.path.join(bad_dir, filename))
                pass


        #Filling log file
        raw_dict = {}
        raw_dict['MJD'] = dic_mjd
        raw_dict['DIT'] = dic_dit
        raw_dict['TARGET'] = dic_targ
        raw_dict['FILT_C'] = dic_filt_c
        raw_dict['FILT_B'] = dic_filt_b
        raw_dict['FLUX_C'] = dic_flux_c
        raw_dict['FLUX_B'] = dic_flux_b
        raw_dict['SIZE'] = window_size
        with (open(output_file, "wb")) as f:
            pickle.dump(raw_dict, f, protocol=-1)

def plot_and_thres(cur_mjd, cur_flux, cam_name, cur_targ, cur_filt):
    """
    Given a time range and the corresponding derived central flux, this function plots them and
    derives the threshold used for frame selection.

    Parameters
    ----------
    **cur_mjd** : MJD np.array

    **cur_flux** : flux np.array

    **cam_name** : Name of the ZIMPOL camera

    **cur_targ** : Current target name

    **cur_filt** : Current filter name

    Returns
    ----------
    **fig** : The figure instance that was generated

    **thres** : the derived threshold

    """

    fig, ax = plt.subplots(1, 1, sharex=True)

    #Plotting flux + date format
    all_times = Time(cur_mjd, format="mjd")
    all_times.format = "iso"
    ax.plot_date(all_times.plot_date, cur_flux, 'ok', ms=2, mec='k')
    fig.autofmt_xdate()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d-%b-%Y\n%H:%M"))

    #Cosmetics
    ax.grid()
    dt = TimeDelta(60.0, format='sec')
    ax.set_xlim((all_times.min()-dt).iso, (all_times.max()+dt).iso)
    ax.set_xlabel('UT')
    ax.set_ylabel('Stellar flux/DIT (ADU/s)')
    ax.set_title(cam_name+' '+cur_targ+" "+cur_filt)

    #Threshold
    if cur_flux.std() > cur_flux.mean()/10.:
        if cur_flux.std() > cur_flux.mean()/2.:
            thres = cur_flux.mean() - 0.5*cur_flux.std()
        else:
            thres = cur_flux.mean() - 2*cur_flux.std()
    else:
        thres = cur_flux.mean() - 5*cur_flux.std()
    #Plotting
    ax.plot(np.array([ax.get_xlim()[0], ax.get_xlim()[1]]), [thres, thres], 'r')

    #Saving figure and closing
    fig.set_size_inches(16, 8)

    return fig, thres


def select_raw_files(raw_dir, do_it, update=True, window_size=10):
    """
    This function actually performs the selection on the  ZIMPOL RAW files.

    For each frame it computes the flux in a 5px*5px square around the maximum flux location. The
    DIT is taken into account. The selection is then performed by deriving the mean and the standard
    deviation for each target and each of its associated filters. From that, each frame with a
    derived flux below a threshold is moved to a "bad" directory that will not be reduced by
    ESOREFLEX.

    The threshold is derived by:\n
    thres = mean - 0.5*std   (if std > mean/10 and std > mean/2, i.e. high dispersion)\n
    thres = mean - 2*std     (if std > mean/10 and std < mean/2, i.e. intermediate dispersion)\n
    thres = mean - 5*sdt     (if std < mean/10, i.e. low dispersion, meaning good data)

    If the same setting is used on several nights, each is considered separately.

    Parameters
    ----------
    **raw_dir** : Directory containing the raw data

    **do_it** : If set to *True*, the degraded files are put in *raw_dir/../bad*

    Keywords
    ----------
    **(update|True)** : At each run of the program, a .pickle file is created with all the
    necesssary information to make the selection. When *update* is True (default), the .pickle file
    is destroyed if it exists and recreated. As opening lots of FITS file is time consumming, when
    doing tests on the same sample of file, one can put *update* to False. This way, the pickle is
    not destroyed and is read. No FITS file is open.

    **(window_size|10)** : Window size on which the flux is derived around the maximum intensity.
    By default, the value is 10, meaning the flux is derived on 10*10px square centered on the
    maximum.

    Side effect
    ----------
    Whether *do_it* is set to True or False, the programs creates plots for each target and each
    filter representing the measured flux as a function of time. The computed threshold is shown as
    a red horizontal line.
    """

    plt.ioff()

    #Input dir
    bad_dir = os.path.join(os.path.join(raw_dir, os.path.pardir), "bad")

    #Creating bad directory
    if not os.path.isdir(bad_dir):
        os.mkdir(bad_dir)
    bad_files = 0

    #Log files
    log_raw = os.path.join(raw_dir, "log_raw.pickle")

    #For each raw file
    if update and os.path.exists(log_raw):
        os.remove(log_raw)
    if not os.path.exists(log_raw):
        archive_raw_info(raw_dir, log_raw, window_size)

    #Reading RAW log
    try:
        with (open(log_raw, "rb")) as f:
            raw_dict = pickle.load(f)
    except IOError:
        return

    #If the pickle file contains flux values with a different window size
    if window_size != raw_dict['SIZE']:
        os.remove(log_raw)
        archive_raw_info(raw_dir, log_raw, window_size)

    #Reading RAW log
    with (open(log_raw, "rb")) as f:
        raw_dict = pickle.load(f)

    #Getting RAW informations
    dic_mjd = raw_dict['MJD']
    dic_dit = raw_dict['DIT']
    dic_filt_c = raw_dict['FILT_C']
    dic_filt_b = raw_dict['FILT_B']
    dic_flux_c = raw_dict['FLUX_C']
    dic_flux_b = raw_dict['FLUX_B']
    dic_targ = raw_dict['TARGET']
    print("\nFound "+str(len(dic_mjd))+" ZIMPOL RAW FITS files")

    #Getting all the files associated to one target
    files_by_targ = {}
    for k, v in dic_targ.items():
        files_by_targ.setdefault(v, []).append(k)

    #Getting all the files associated to one a filter
    files_by_filt_c = {}
    files_by_filt_b = {}
    for cc, bb in zip(dic_filt_c.items(), dic_filt_b.items()):
        files_by_filt_c.setdefault(cc[1], []).append(cc[0])
        files_by_filt_b.setdefault(bb[1], []).append(bb[0])

    #For each target
    for cur_targ in files_by_targ:

        #Camera sets
        both_cams = [[files_by_filt_c, dic_flux_c, 'Callas'], \
        [files_by_filt_b, dic_flux_b, 'Bertoli']]

        #For each camera
        for cam in both_cams:

            #Camera data
            files_by_filt = cam[0]
            dic_flux = cam[1]
            cam_name = cam[2]

            #For each filter
            for cur_filt in files_by_filt:
                #Getting file list
                cur_list = files_by_filt[cur_filt]

                #Initializing
                all_mjd = np.array([])
                all_flux = all_mjd.copy()

                #Exploring files
                for cur_file in cur_list:

                    #If target matches
                    if dic_targ[cur_file] == cur_targ:

                        all_mjd = np.append(all_mjd, dic_mjd[cur_file]*np.ones(\
                        dic_flux[cur_file].shape))
                        all_flux = np.append(all_flux, dic_flux[cur_file]/dic_dit[cur_file])

                #Sorting chronologically
                all_flux = all_flux[all_mjd.argsort()]
                all_mjd.sort()

                #If we have observations
                if all_mjd.size > 0:

                    #Only one night : easy case
                    if (all_mjd.max() - all_mjd.min()) < 0.5:
                        fig, thres = plot_and_thres(all_mjd, all_flux, cam_name, cur_targ, cur_filt)
                        #Saving/closing figure
                        fig.savefig(os.path.join(raw_dir, cam_name+' '+cur_targ+" "+cur_filt+".png"),\
                        bbox_inches='tight', dpi=300)
                        plt.close(fig)

                    #Else, several night, we treat each night separately
                    else:

                        #Preparation for the exploration of each nights
                        remain_mjd = all_mjd
                        remain_flux = all_flux
                        print(cam_name+' '+cur_targ+" "+cur_filt)
                        nn = 0

                        #While we have not explored all nights
                        while remain_mjd.size > 0:

                            #1st of the remaining nights
                            tag_night = (remain_mjd < remain_mjd.min()+0.5)
                            nn += 1

                            #Dealing with it
                            fig, thres = plot_and_thres(remain_mjd[tag_night], remain_flux[tag_night], \
                            cam_name, cur_targ, cur_filt)

                            #Eliminating the night that was just dealt with
                            remain_mjd = remain_mjd[~tag_night]
                            remain_flux = remain_flux[~tag_night]

                            #Saving/closing figure
                            fig.savefig(os.path.join(raw_dir, cam_name+' '+cur_targ+" "+\
                            cur_filt+"_"+str(nn)+".png"),\
                            bbox_inches='tight', dpi=300)
                            plt.close(fig)

                    #Selection of good/bad files
                    for cur_file in cur_list:

                        if dic_targ[cur_file] == cur_targ and \
                        (dic_flux[cur_file]/dic_dit[cur_file] <= thres).any():

                            try:
                                #Only if asked
                                if do_it:
                                    os.rename(os.path.join(raw_dir, cur_file), \
                                    os.path.join(bad_dir, os.path.basename(cur_file)))
                                bad_files += 1
                            #File already moved
                            except OSError:
                                continue

    #Result of the selection
    if do_it:
        print("Moved "+str(bad_files)+" bad files")
    else:
        print("Found "+str(bad_files/2)+" bad files, not moved")

    plt.ion()
