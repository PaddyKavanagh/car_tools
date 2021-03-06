#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# @author: Patrick Kavanagh (DIAS)
#
"""

Tools to create a MIRISim scene object from an input catalogue.

Can either create a randomly generated catalogue or use Alistair's data on
the JWST Astrometric Field (exported from xlxs to csv after deleting the empty
column - possible to read xlsx file directly but it is significantly slower).

Can also be used as a command line script to generate a scene.ini file for
the JWST Astronometric catalogue:

python scene_from_cat.py inputfile ra dec

    inputfile   --  the path+name of the catalogue file

    instrument  --  MIRI instrument being simulated (IMA, MRS or LRS)

    ra          --  target RA in decimal degrees

    dec         --  target Dec in decimal degrees

Intended for simulations for CAR007 (Distortion and plate scale)

"""

from __future__ import absolute_import, division, print_function
import optparse
import os, sys, time

import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.optimize import curve_fit

from mirisim.skysim import *
from mirisim.config_parser import *


def generate_random_position_flux(ra, dec):
    """
    Generate random RA and DEC within 1 arcmin of supplied RA and DEC. Generate
    a flux in uJy.

    Parameters:
    ra, dec         --  centre position (in deg) around which random position
                        generated

    Returns:
    src_ra, src_dec --  random source position within 1 arcmin of ra, dec

    src_flux        --  random source flux
    """
    # generate random RAs and Decs in a 20x20 arcmin region
    src_ra = random.uniform(ra-0.017, ra+0.017)
    src_dec = random.uniform(dec-0.017, dec+0.017)

    # for fluxes, do randomisation in log space for even distribution
    log_src_flux = random.uniform(2, 7.)
    src_flux = np.exp(log_src_flux)

    return src_ra, src_dec, src_flux


def generate_random_offset_flux():
    """
    Generate random offsets within 1 arcmin 0,0. Generate a flux in uJy.

    Returns:
    x_off, y_off    --  random source position within 1 arcmin of ra, dec

    src_flux        --  random source flux
    """
    x_off = random.uniform(-60., 60.)
    y_off = random.uniform(-60., 60.)

    # for fluxes, do randomisation in log space for even distribution
    log_src_flux = random.uniform(1, 6.)
    src_flux = np.exp(log_src_flux)

    return x_off, y_off, src_flux


def generate_test_cat(source_num=None, centre_coords=[0.0, 0.0], save_file=False, out_file='test_catalogue.txt'):
    """
    Function to generate a test catalogue with the format:

    RA      DEC     Flux
    ===     ===     ===

    Can either save out to file or return array. Fluxes are in the 10 uJy to 1 Jy range

    Parameters:
    source_num      --  optionally set the number of sources to be included
                        in the catalogue
    centre coords   --  optionally pass coordinates as a list around which
                        sources will be created in the simulated catalogue.
                        Assumed that coords are in decimal degrees
                        (default = [0.0,0.0])
    save_file       --  if true, saves the catalogue to file
    out_file        --  if save_file is True, sets file name of catalogue
                        (default = 'test_catalogue.txt')

    Returns:
    cat             --  array containing test catalogue data

    """
    assert 0. <= centre_coords[0] <= 360., 'RA is out of allowable range, check units'
    assert -90. <= centre_coords[1] <= 90., 'DEC is out of allowable range, check units'

    # generate random data for n sources and write to list
    cat = np.zeros((source_num,3))
    for n in range(source_num):
        cat[n,0], cat[n,1], cat[n,2] = generate_random_position_flux(centre_coords[0], centre_coords[1])

    if save_file == True:
        # save a file if requested (delete old file if exists)
        try:
            os.remove(out_file)
        except OSError:
            pass

        np.savetxt(out_file, cat, fmt='%0.4e', header='RA\tDec\tFlux\ndeg\tdeg\tuJy')

    return cat


def read_cat_from_txt(txt_cat):
    """
    Simple function to read an input catalogue and return as a numpy array.

    Parameters:
    txt_cat     --      txt file containing catalogue data. Assumed to contain
                        RA, DEC, and flux columns (space seperated)

    Returns:
    catalogue   --      numpy array with catalogue data

    """
    # open file and read contents
    with open(txt_cat, 'rb') as catfile:
        data = np.loadtxt(catfile, delimiter='', comments='#', unpack=True)

        # check for the correct number of columns (should be 3: RA, DEC, flux)
        assert data.shape[0] == 3, 'Incorrect number of columns in catalogue file'

    return data


def get_nearby_sources(cat, target_ra, target_dec, instrument):
    """
    Simple function to read an input catalogue from an txt file, select
    sources within 2 arcmin of target (given by ra and dec)
    and return as a numpy array.

    Parameters:
    xls_cat     --      xls file containing catalogue data. For now, this is
                        taylored to read from Alistair's xls file.

    target_ra   --      ra of target source in decimal degrees

    target_dec  --      dec of target source in decimal degrees

    instrument  --      MIRI instrument being simulated (IMA, LRS or MRS)

    Returns:
    catalogue   --      numpy array with catalogue data

    """
    # open file and read contents
    with open(cat, 'rb') as catfile:
        data = np.loadtxt(catfile, delimiter=',', comments='#', unpack=True, skiprows=4)

    # check for the correct number of columns (should be 3: RA, DEC, flux)
    assert data.shape[0] == 14, 'Incorrect number of columns in catalogue file, should be 14'

    nearby_cat = []
    for n in range(0, data.shape[1]-1):

        #print(data[:,n])

        # determine distance of source to target
        dist = np.sqrt( (data[1,n] - target_ra)**2 + (data[2,n] - target_dec)**2 )

        # if distance <2' for IMA and <10" for MRS, write to output array
        if instrument == 'IMA':
            if dist < 0.033:
                nearby_cat.append([data[:,n]])

        if instrument == 'MRS':
            if dist < 0.003:
                nearby_cat.append([data[:,n]])

        if instrument == 'LRS':
            if dist < 0.01:
                nearby_cat.append([data[:,n]])

    return np.squeeze(np.asarray(nearby_cat))


def radec_to_offset(target_ra, target_dec, src_coords, test=None):
    """
    Function to take catalogue RA and DEC values and convert these to offsets
    from a specified coordinate (all in dec degrees)

    Parameters:
    target_ra, target_dec   --  coordinates of the target

    src_coords              --  coordinates of the sources in the catalogue

    test                    --  boolean switch to create test plots of source
                                coordinates and offsets

    Returns:
    src_offsets             --  array containing src offsets from target coords

    """
    # check the target coords are in the correct range for dec degrees
    assert 0. <= target_ra <= 360., 'Target RA is out of allowable range, check units'
    assert -90. <= target_dec <= 90., 'Target DEC is out of allowable range, check units'

    # assume sources are on a Cartesian grid an calculate offsets in arcsec
    src_offsets = np.zeros((src_coords.shape))
    for n in range(src_coords.shape[0]):
        src_offsets[n,0] = (src_coords[n,0] - target_ra) * 3600.
        src_offsets[n,1] = (src_coords[n,1] - target_dec) * 3600.

    # for testing, create a plot of the input RAs and DECs and the offsets
    if test == True:
            # Setup plot
            fig, axs = plt.subplots(1,2, figsize=(10,5))

            # ra, dec
            axs[0].plot(src_coords[:,0],src_coords[:,1], marker='.',
                        markersize=5, color='b', linewidth=0)
            axs[0].plot(target_ra,target_dec, marker='o', markersize=6,
                        color='r', linewidth=0)
            axs[0].set_xlabel('RA (degrees)')
            axs[0].set_ylabel('DEC (degrees)', rotation=90)

            # offsets
            axs[1].plot(src_offsets[:,0],src_offsets[:,1], marker='.',
                        markersize=5, color='b', linewidth=0)
            axs[1].plot([0],[0], marker='o', markersize=6, color='r',
                        linewidth=0)
            axs[1].set_xlabel('RA offset (arcsec)')
            axs[1].set_ylabel('DEC offset (arcsec)', rotation=90)

            plt.show()

    return src_offsets


def make_scene_obj(offsets,fluxes,wref,temp):
    """
    Function to take catalogue and construct a MIRISim scene object

    Note that this assumes all sources are points and described by
    blackbody SED.

    Parameters:
    offsets                 --  array containing source offsets

    fluxes                  --  array containing source fluxes (in mJy)

    wref                    --  reference wavelength of flux

    temp                    --  blackbody temperature

    Returns:
    scene_obj               --  catalogue scene object

    """
    # set background
    background = Background(level='low', gradient=5., pa=15.0, centreFOV=(0., 0.))

    # process the sources in the catalogue
    targets = []
    for n in range(offsets.shape[0]):
        sed = BBSed(Temp=temp, wref=wref, flux=fluxes[n]*1.e3)
        point = Point(Cen=(offsets[n,0],offsets[n,1]), vel=0.0)
        point.set_SED(sed)
        targets.append(point)

    # create the object
    scene_obj = SceneConfig.makeScene(loglevel=0, background=background, targets=targets)

    return scene_obj


def make_simple_cat_scene_obj(cat_file=None,target_coords=[1.0,1.0],random=True,\
                              source_num=100, centre_coords=[1.0,1.0],save_file=False, \
                              out_file='my_catalogue.txt', save_ini=False):
    """
    Generate the scene object from a simple input catalogue or a randomly generated
    catalogue, containing RA, DEC and flux for each source.

    Note that this assumes all sources are points and described by
    blackbody SED with T=300 K at wref=10 micron

    Parameters:
    cat_file                --  catalogue file. Assumed to be txt file with
                                three columns: RA,DEC,flux
                                (default = None)

    target_coords           --  set the coordinates of the observation target
                                in decimal degrees
                                (default = [1.0,1.0])

    random                  --  if True, generates a random catalogue
                                (default = True)

    source_num              --  if random=True, number of sources to generate
                                (default = 100)

    centre_coords           --  if random=True, centre point in RA,DEC (decimal
                                degrees) of the catalogue field
                                (default = [1.0,1.0])

    save_file               --  if random=True, switch to save the catalogue file
                                (default = False)

    out_file                --  if random=True, and save_file=True, name of file
                                to save the catalogue
                                (default = 'my_catalogue.txt')

    save_ini                --  if True, saves the scene object to an ini file
                                called 'scene.ini'
                                (default = False)

    Returns:
    scene_obj               --  catalogue scene object

    """
    # if random = True, generate random data
    if random == True:
        cat = generate_test_cat(source_num=source_num, centre_coords=centre_coords,
                                save_file=save_file, out_file=out_file)

    # convert source coordinates to offsets required for MIRISim
    offsets = radec_to_offset(target_coords[0], target_coords[1], cat[:,:2])

    # generate the scene object
    cat_scene_obj = make_scene_obj(offsets,cat[:,2], wref=10., temp=300.)

    # save ini file if requested (delete old one of exists)
    try:
        os.remove('scene.ini')
    except OSError:
        pass
    if save_ini == True: cat_scene_obj.write('scene.ini')

    return cat_scene_obj


def make_cat_scene_obj(cat_file=None,target_coords=[80.5,-69.5], instrument='IMA',save_ini=False):
    """
    Generate the scene object from an input csv catalogue exported from
    Alistair's Excel file

    Parameters:
    cat_file                --  catalogue csv file
                                (default = None)

    target_coords           --  set the coordinates of the observation target
                                in decimal degrees
                                (default = [80.5,-69.5])

    instrument              --  MIRI instrument being simulated. Used to determine
                                how to spatially filter catalogue (IMA,LRS or MRS)

    save_ini                --  if True, saves the scene object to an ini file
                                called 'scene.ini'
                                (default = False)

    Returns:
    scene_obj               --  catalogue scene object

    """
    # search the full astrometric catalogue for sources near the target
    cat = get_nearby_sources(cat_file, target_coords[0], target_coords[1], instrument)

    # convert source coordinates to offsets required for MIRISim
    offsets = radec_to_offset(target_coords[0], target_coords[1], cat[:,1:3])

    # generate the scene object using fluxes at 10 micron and assuming hot stars
    cat_scene_obj = make_scene_obj(offsets,cat[:,7],wref=10.,temp=10000.)

    # save ini file if requested (delete old one of exists)
    if save_ini == True:
        ini_file_name_stem = os.path.splitext(os.path.basename(cat_file))[0]
        ini_file_name = ini_file_name_stem + '_' + instrument + '_' + str(target_coords[0]) \
                        + '_' + str(target_coords[1]) + '.ini'
        try:
            os.remove(ini_file_name)
        except OSError:
            pass
        cat_scene_obj.write(ini_file_name)

    return cat_scene_obj


if __name__ == "__main__":
        # Parse arguments
        help_text = __doc__
        usage = "\n\n%prog inputfile ra -- dec\n"
        usage += "\nGenerate a scene.ini file from a csv file exported from "
        usage += "Alistair's JWST Astrometric Field Excel spreadsheet."

        parser = optparse.OptionParser(usage)
        (options,args) = parser.parse_args()

        # check for correct number of arguments
        try:
            input_file = args[0]
            instrument = args[1]
            target_ra = float(args[2])
            target_dec = float(args[3])

        except IndexError:
            print(help_text)
            time.sleep(1) # Ensure help text appears before error messages.
            parser.error("Not enough arguments provided")
            sys.exit(1)

        # check the catalogue file exists
        try:
            file = open(input_file, 'r')

        except IOError:
            parser.error("Catalogue file not found")
            sys.exit(1)

        # check the instrument
        try:
            assert instrument in ['IMA','MRS','LRS']

        except AssertionError:
            parser.error("Instrument not recognised, must be IMA, MRS or LRS")
            sys.exit(1)

        # check RA is in astrometric field range
        try:
            assert (80. <= target_ra) and (target_ra <= 81.)

        except AssertionError:
            parser.error("RA is not in astrometric field range of 80 .. 81 degrees")
            sys.exit(1)

        # check Dec is in astrometric field range
        try:
            assert (-70. <= target_dec) and (target_dec <=-69.)

        except AssertionError:
            parser.error("Dec is not in astrometric field range of -69 .. -70 degrees")
            sys.exit(1)

        #make_cat_scene_obj('/Users/patrickkavanagh/CARs/JA_LMC_MIRI.csv', target_coords=[80.508612, -69.510760], save_ini=True)
        make_cat_scene_obj(cat_file=input_file, target_coords=[target_ra, target_dec],
                            instrument=instrument, save_ini=True)
