#! /usr/bin/env python
#
# rotate_dataset.py
# rotate arrays in all integrations for SCI,DQ, and ERR extensions
# in level 1b file (which has CRs added), and write rotated arrays
# to file for input into ramp_fit_arr.py
#

import sys
import os
from . import rotate_utils
from astropy.io import fits

def get_num_ints(in_file):
    """
    Short Summary
    -------------
    Read number of integrations from input file.

    Parameters
    ----------
    in_file: string
      input file

    Returns
    _______
    num_ints: integer
       number of integrations

    """
    fh = fits.open(in_file, do_not_scale_image_data=True)
    try:
        num_ints = fh['SCI'].header['NAXIS4'] # integrations
    except Exception as errmess:
        print('Fatal ERROR: ', errmess)
    fh.close()
    return int(num_ints)


def write_keys(in_file, out_file):
    """
    Short Summary
    -------------
    Write input file name and ROTATED=True in output file

    Parameters
    ----------
    in_file: string
      input file

    out_file: string
      output file

    Returns
    _______

    """
    try:
        fh_out = fits.open(out_file, mode='update', do_not_scale_image_data=True)
        fh_out[0].header['INFILE'] = (in_file, 'prerotated file')
        fh_out[0].header['ROTATED'] = (True, 'dataset is rotated')
        fh_out.close()
    except Exception as errmess:
        print('Fatal ERROR: ', errmess)


def delete_temp_file(fname, index):
    """
    Short Summary
    -------------
    Delete a temp file created

    Parameters
    ----------
    fname: string
      file name

    index: integer
      designation of integration

    Returns
    _______

    """
    fname = fname + str(index) + '.fits'
    os.remove(fname)


def delete_all_temp_files(num_ints):
    """
    Short Summary
    -------------
    Delete all temporary files created

    Parameters
    ----------
    num_ints: integer
      number of integrations

    Returns
    _______

    """
    for ii in range(num_ints):
        delete_temp_file('temp_DQ', ii)
        delete_temp_file('temp_SCI', ii)
        delete_temp_file('temp_ERR', ii)
        delete_temp_file('temp_r_DQ', ii)
        delete_temp_file('temp_r_SCI', ii)
        delete_temp_file('temp_r_ERR', ii)

    os.remove('temp_all_SCI.fits')
    os.remove('temp_all_ERR.fits')
    os.remove('temp_all_DQ.fits')
    os.remove('temp_all_ext.fits')


def do_all_rot_tasks(in_file, num_ints, out_file):
    """
    Short Summary
    -------------

    Extract and rotate all extensions from input file and
    write to output file.

    Parameters
    ----------
    in_file: string
      input file name

    num_ints: integer
      number of integrations

    out_file: string
      output file name

    Returns
    _______
    """

    print('Starting data rotation and writing.')
    print(' Input parameters are:')
    print('  in_file:', in_file)
    print('  out_file:', out_file)
    print('  num_ints:', num_ints)

    rotate_utils.extract_all_ext(in_file)
    rotate_utils.rot_all_ext(num_ints)
    rotate_utils.glue_all_ext(num_ints)
    rotate_utils.copy_input_hdr(in_file, out_file)

    print(' Done rotating input and writing output dataset.')


if __name__ == "__main__":
    """ get input and output files and do rotation
    """
    if (len(sys.argv) < 3):
        print('Fatal ERROR: ')
        print('usage: ./rotate_dataset.py in_file out_file')

    in_file = sys.argv[1]  # file with data to rotate, and header to copy
    out_file = sys.argv[2] # output file containing input header and rotated data

    try:
        num_ints = get_num_ints(in_file)
        do_all_rot_tasks(in_file, num_ints, out_file)
        write_keys(in_file, out_file)
        delete_all_temp_files(num_ints)
    except Exception as errmess:
        print('Fatal ERROR: ', errmess)
