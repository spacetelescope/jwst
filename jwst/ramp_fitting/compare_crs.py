#! /usr/bin/env python
#
# compare_crs.py - compare true and found cosmic rays
# pragma: no cover
import sys
import time
import numpy as np
from astropy.io import fits


def do_comparison(found_file, created_file):
    """
    Brief Summary
    ---------
    Open and read data from found and created files, then compare the locations
    of the created to the found cosmic rays and output arrays for the created
    only, found only, neither created nor found, and both created and found
    cosmic rays for both the full 3d or 4d input datasets and the 2d arrays.
    Data in the created file is in the 0th extension.  Data in the found file
    is in the 1st extension; the found file is output by ramp_fit.py.  For MIRI,
    copies are concatenated of the created data (except for the last frame) -
    this is done because the MC routines that produced the cosmic rays in the
    created file creates the same cosmic rays for each integration for simplicity.

    Parameters
    ---------
    found_file: string
       name of file with found crs

    created_file: string
       name of file with created crs

    Returns
    ---------
    """

    fh_f, fh_c, data_f, data_c = get_data(found_file, created_file)

    print('Initial found data shape ', data_f.shape)
    print(' and created data shape= ', data_c.shape)

    # Compare slice i of created to slice i+1 in found
    if (data_f.shape[0] == 1):  #  NIRCAM
        data_f = data_f[0, :, :, :]
        if (data_c.shape[0] == 1):  # to accept output of mc_4d
            data_c = data_c[0, :, :, :]
        data_c_start = data_c[:-1, :, :]
        data_f_end = data_f[1:, :, :]
    elif (fh_f['SCI'].header['NAXIS'] == 3):  #  NIRSPEC
        data_c_start = data_c[:-1, :, :]
        data_f_end = data_f[1:, :, :]
    elif (data_f.shape[0] > 1 and fh_f['SCI'].header['NAXIS'] == 4):  # MIRI
        # concatenate copies of created data (except for the last frame)
        num_ints = int(fh_f[1].data.shape[0]) # number of integrations
        data_c_start = (np.repeat(data_c[:-1, :, :], num_ints, axis=0))
        data_f_end = data_f[:, 1:, :, :]
        data_c_start = data_c_start.reshape(data_f_end.shape)
    else:
        print(' FATAL ERROR - unsupported instrument')

    print('Truncated found data shape ', data_f_end.shape)
    print(' and truncated created data shape= ', data_c_start.shape)
    try:
        assert(data_f_end.shape == data_c_start.shape)
    except AssertionError:
        print(' FATAL ERROR: adjusted found data shape ', data_f.shape, \
              ' is not the same as adjusted created data shape= ', data_c.shape)

    neither = (data_c_start == 0.) & (data_f_end == 0.)
    both = (data_c_start != 0.) & (data_f_end != 0.) # created CR was found
    c_only = (data_c_start != 0.) & (data_f_end == 0.) # created CR not found
    f_only = (data_c_start == 0.) & (data_f_end != 0.) # found CR was not created

    try:
        assert(neither.sum() + both.sum() + c_only.sum() + f_only.sum() \
                == data_c_start.size)
    except AssertionError:
        print('FATAL ERROR: sum of components must equal total number of pixels ')

    print(' Within the input dataset cubes:')
    print('   Number of created but not found pixels: ', c_only.sum())
    print('   Number of found but not created pixels: ', f_only.sum())
    print('   Number of pixels that are both found and created: ', both.sum())
    print('   Number of pixels that are neither found nor created: ', neither.sum())
    print('   ')
    print('  The fraction of all pixels that were found only: ', \
        float(f_only.sum()) / float(data_c_start.size))
    print('  The fraction of all pixels that were created only: ', \
        float(c_only.sum()) / float(data_c_start.size))
    print('  The fraction of pixels in the created file having cosmic rays:', \
        float(c_only.sum()) / (data_c_start.shape[-2] * data_c_start.shape[-1]))
    print('   ')

    write_files(neither, both, c_only, f_only, fh_c, data_c_start)

def get_data(found_file, created_file):
    """
    Brief Summary
    ---------
    Open found and created files and get data

    Parameters
    ---------
    found_file: string
       name of file with found crs

    created_file: string
       name of file with created crs

    Returns
    ---------
    fh_f: fits file handle
       handle for found file

    fh_c: fits file handle
       handle for created file

    data_f: int array
       data in found file, 1 indicates a found CR

    data_c: float array
       data in created file, magnitudes of created CRs
    """

    try:
        fh_f = fits.open(found_file)
        print('Found file has: ', fh_f.info())
    except Exception:
        print(' FATAL ERROR: Unable to open found file ', found_file)

    try:
        fh_c = fits.open(created_file)
        print('Created file has: ', fh_c.info())
    except Exception:
        print(' FATAL ERROR: Unable to open created file ', created_file)

    try:
        data_f = fh_f['SCI'].data
    except Exception:
        print(' FATAL ERROR: data for found data was expected in SCI extension')

    try:
        data_c = fh_c['SCI'].data
    except Exception:
        try:
            data_c = fh_c[0].data
        except Exception:
            print(' FATAL ERROR: created data expected in either SCI or 0 extensions')

    return fh_f, fh_c, data_f, data_c

def write_files(neither, both, c_only, f_only, fh_c, data_c_start):
    """
    Brief Summary
    ---------
    Write to files the calculated arrays

    Parameters
    ---------
    neither: int array
       input cube pixels having neither a created or found crs

    both: int array
       input cube pixels having both a created and found cr

    c_only: int array
       input cube pixels having a created but not found cr

    f_only: int array
       input cube pixels having a found but not created cr

    fh_c: fits file handle
       handle for created file

    data_c_start: float array
       data in created file - all but final frame in each integration

    Returns
    ---------
    """

    # output arrays for all pixels in input datasets
    write_to_file(neither.astype(np.int16), 'neither_cube.fits')
    write_to_file(both.astype(np.int16), 'both_cube.fits')
    write_to_file(c_only.astype(np.int16), 'c_only_cube.fits')
    write_to_file(f_only.astype(np.int16), 'f_only_cube.fits')

    # output arrays for pixels in 2d array
    print(' Within the 2d arrays:')
    if (fh_c[0].header['NAXIS'] == 3): # for nirspec with 1 integration
        write_to_file(neither.sum(axis=0), 'neither_2d.fits')
        write_to_file(both.sum(axis=0), 'both_2d.fits')
        write_to_file(c_only.sum(axis=0), 'c_only_2d.fits')
        write_to_file(f_only.sum(axis=0), 'f_only_2d.fits')
        print(' The fraction of pixels in the 2d array having true CRs:',\
              float(len(np.where(both.sum(axis=0) != 0.)[0])) / data_c_start.size)
    elif (fh_c[1].header['NAXIS'] == 4): # for miri or nircam cases
        write_to_file(neither.sum(axis=1).sum(axis=0), 'neither_2d.fits')
        write_to_file(both.sum(axis=1).sum(axis=0), 'both_2d.fits')
        write_to_file(c_only.sum(axis=1).sum(axis=0), 'c_only_2d.fits')
        write_to_file(f_only.sum(axis=1).sum(axis=0), 'f_only_2d.fits')
        print(' The fraction of pixels in the 2d array having true CRs:',\
              float(len(np.where(both.sum(axis=1).sum(axis=0) != 0.)[0])) / \
              data_c_start.size)
    else:
        print('FATAL ERROR - unexpected case in write_file()')

def write_to_file(data, filename):
    """
    Brief Summary
    ---------
    Write the specified data to the specified file name

    Parameters
    ---------
    data: array
       array to write

    filename: string
       file being written

    Returns
    ---------
    """
    fimg = fits.HDUList()
    fimghdu = fits.PrimaryHDU()
    fimghdu.data = data
    fimg.append(fimghdu)
    fimg.writeto(filename, overwrite=True)
    print(' wrote output data to: ', filename)


if __name__ == "__main__":
    """Get found and created files, and call compare_crs.
    """
    usage = "usage: ./compare_crs.py found_file created_file"
    print(' Starting compare_crs at time:', time.asctime())

    found_file = sys.argv[1]
    created_file = sys.argv[2]

    try:
        do_comparison(found_file, created_file)
    except Exception as errmess:
        print('FATAL ERROR: ', errmess)

    print('   ')
    print(' DONE at time:', time.asctime())
