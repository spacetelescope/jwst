#! /usr/bin/env python
# rotate_utils.py
# utility functions for creating a rotated dataset

import os
import shutil
import numpy as np
import Rotatefitscube
from astropy.io import fits

def write_file(data, bitpix, output_fname):
    """
    Short Summary
    -------------
    Write data and bitpix keyword to specified file

    Parameters
    ----------
    data: float
       array

    bitpix: integer
      value of bitpix keyword

    output_fname: string
      name of output file

    Returns
    _______

    """
    fitsobj = fits.HDUList()
    hdu = fits.PrimaryHDU()
    prihdr = hdu.header
    prihdr['BITPIX'] = bitpix
    hdu.data = data
    fitsobj.append(hdu)
    fitsobj.writeto(output_fname, overwrite=True)
    fitsobj.close()


def extract_ext_and_write(fh_in, ext):
    """
    Short Summary
    -------------
    Extract given extension for all integrations from input file and
       write to output file

    Parameters
    ----------
    fh_in: file handle
      handle of input file

    extension: string
      extension to extract

    Returns
    _______

    """
    hdr = fh_in[ext].header
    bitpix = hdr['BITPIX']
    data_ext = fh_in[ext].data
    ninteg = data_ext.shape[0]

    for ii_int in range(ninteg):
        data_integ = data_ext[ii_int, :, :, :]
        outfile_name = 'temp_' + ext + str(ii_int) + '.fits'
        write_file(data_integ, bitpix, outfile_name)


def extract_all_ext(input_file):
    """
    Short Summary
    -------------
    Extract all extensions from input file and write to extension-specific files

    Parameters
    ----------
    input_fname: string
       name of input file

    Returns
    _______
    """
    fh = fits.open(input_file, do_not_scale_image_data=True)
    try:
        extract_ext_and_write(fh, 'SCI')
        extract_ext_and_write(fh, 'ERR')
        extract_ext_and_write(fh, 'DQ')

    except Exception as errmess:
        print('Fatal ERROR: ', errmess)
    fh.close()


def copy_input_hdr(hdr_file, out_file):
    """
    Short Summary
    -------------
    Create outputfile with same header

    Parameters
    ----------
    hdr_file: string
      file supplying hdr info

    out_file: string
      output file

    Returns
    _______
    """
    try:
        shutil.copy(hdr_file, out_file)
        copy_data('temp_all_ext.fits', out_file)
    except Exception as errmess:
        print('Fatal ERROR: ', errmess)


def copy_data(data_file, hdr_file):
    """
    Short Summary
    -------------
    Update hdr_file (which has the desired header) with data from data_file

    Parameters
    ----------
    data_file: string
      file supplying data

    hdr_file: string
       output file, already with correct header

    Returns
    _______

    """
    fh_data = fits.open(data_file, do_not_scale_image_data=True)
    sci_data = fh_data[1].data
    err_data = fh_data[2].data
    dq_data = fh_data[3].data

    fh_hdr = fits.open(hdr_file, mode='update', do_not_scale_image_data=True)
    fh_hdr[1] = fits.ImageHDU(sci_data, fh_hdr[1].header)
    fh_hdr[2] = fits.ImageHDU(err_data, fh_hdr[2].header)
    fh_hdr[3] = fits.ImageHDU(dq_data, fh_hdr[3].header)

    fh_hdr.flush(output_verify='fix')
    fh_data.close()
    fh_hdr.close()


def glue_an_ext(prefix, nfiles, out_file):
    """
    Short Summary
    -------------
    Add data cube from extension of input to hypercube of output

    Parameters
    ----------
    prefix: string
       type of extension

    nfiles: integer
       number of files

    out_file: string
       output file

    Returns
    _______
    """

    nfiles = int(nfiles)
    in_file = prefix + '0.fits'
    fh_0 = fits.open(in_file, do_not_scale_image_data=True)
    cube_shape = fh_0[0].data.shape
    hypercube_shape = (nfiles,) + cube_shape # for output
    fh_0.close()

    if (prefix[-2:] == 'DQ'): # DQ extensions will have int16
        hyper_cube_ext = np.zeros(hypercube_shape, dtype=np.int16)
    else: # SCI and ERR extensions will have float32
        hyper_cube_ext = np.zeros(hypercube_shape, dtype=np.float32)

    for ii_int in range(nfiles):
        in_file = prefix + str(ii_int) + '.fits'
        fh_in = fits.open(in_file, do_not_scale_image_data=True)
        hyper_cube_ext[ii_int, :, :, :] = fh_in[0].data
        fh_in.close()

    append_data(hyper_cube_ext, out_file)


def glue_all_ext(nfiles):
    """
    Short Summary
    -------------
    Combine extracted and rotated extensions from all integrations into a file for each extension,
      then glue each of these into a single file

    Parameters
    ----------
    nfiles: integer
       number of files (one per integration)

    Returns
    _______
    """
    try:
        glue_an_ext('temp_r_SCI', nfiles, 'temp_all_SCI.fits')
        glue_an_ext('temp_r_ERR', nfiles, 'temp_all_ERR.fits')
        glue_an_ext('temp_r_DQ', nfiles, 'temp_all_DQ.fits')
        glue_all()

    except Exception as errmess:
        print('Fatal ERROR: ', errmess)


def glue_all():
    """
    Short Summary
    -------------
    Append all extension types to a single file 'temp_all_ext.fits'

    Parameters
    ----------

    Returns
    _______
    """

    fimg = fits.HDUList()  # output file ('temp_all_ext.fits')
    fimghdu = fits.PrimaryHDU()
    fimg.append(fimghdu)

    try:
        append_ext('temp_all_SCI.fits', fimg)
        append_ext('temp_all_ERR.fits', fimg)
        append_ext('temp_all_DQ.fits', fimg)
    except Exception as errmess:
        print('Fatal ERROR appending extension: ', errmess)

    try:
        os.remove("temp_all_ext.fits")
    except FileNotFoundError:
        pass

    fimg.writeto("temp_all_ext.fits")
    fimg.close()

def append_ext(in_file, fimg_out):
    """
    Short Summary
    -------------
    Append an extension's data to input file

    Parameters
    ----------
    in_file: string
       file to append

    fimg_out: hdulist
       hdulist for output file

    Returns
    _______

    """
    fh_in = fits.open(in_file, do_not_scale_image_data=True)
    fimg_out_hdu = fits.ImageHDU()
    fimg_out_hdu.data = fh_in[0].data
    fimg_out.append(fimg_out_hdu)
    fh_in.close()


def append_data(data, out_file):
    """
    Short Summary
    -------------
    Append data array to file

    Parameters
    ----------
    data: float
      data to append to output file

    out_file: string
       output file name

    Returns
    _______

    """
    fitsobj = fits.HDUList()
    hdu = fits.PrimaryHDU()
    hdu.data = data
    fitsobj.append(hdu)
    fitsobj.writeto(out_file, overwrite=True)


def rot_all_ext(nfiles):
    """
    Short Summary
    -------------
    For all extensions, rotate each of the input files and write to the appropriate
    output file

    Parameters
    ----------
    nfiles: integer
       number of each extension type

    Returns
    _______

    """

    try:
        rotate_ext('temp_SCI', 'temp_r_SCI', nfiles)
        rotate_ext('temp_ERR', 'temp_r_ERR', nfiles)
        rotate_ext('temp_DQ', 'temp_r_DQ', nfiles)

    except Exception as errmess:
        print('Fatal ERROR: ', errmess)

def rotate_ext(in_prefix, out_prefix, nfiles):
    """
    Short Summary
    -------------
    Rotate the only extension in each of the input files and write to each output file


    Parameters
    ----------
    in_prefix : string
       extension type, prefix of input file

    out_prefix : string
       extension type, prefix of output file

    nfiles: integer
       number of each extension type

    Returns
    _______
    """

    for ii in range(int(nfiles)):
        in_name = in_prefix + str(ii) + '.fits'
        out_name = out_prefix + str(ii) + '.fits'
        Rotatefitscube.Rotatefitscube(in_name, out_name)
