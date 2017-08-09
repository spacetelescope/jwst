import glob
from astropy.io import fits
from astropy import wcs
from astropy.wcs import InvalidTransformError
import six

"""
More Kw to consider 
radesys - I think it should be moved - "OTHER" in HST
velocity aberrration related kw - move them
  DVA_RA< DVA_DEC, VA_SCALE

TARG_RA, TARG_DEC
"""
stsci_wcs_kw = ['PA_V3', 'V2_REF', 'V3_REF', 'PA_APER', 'VPARITY',
                'RA_V1', 'DEC_V1', 'ROLL_REF', 'RA_REF', 'DEC_REF']

wcslib_kw_to_remove = ['LONPOLE', 'LATPOLE', 'MJD-OBS', 'DATE-OBS']


def move_wcs(files):
    for name in files:
        print('Working on file {0}.'.format(name))
        f = fits.open(name, mode='update')
        if f[0].header['telescop'] != 'JWST':
            f.close()
            continue

        new_hdr = _collect_wcs_keywords(f)
        update_sci_extension(f, new_hdr)
        clean_primary_header(f, new_hdr)
        f.close()
    
    
def _collect_wcs_keywords(f):
    # Convert the WCS to a header
    try:
        w = wcs.WCS(f[0].header)
    except InvalidTransformError:
        if f[0].header['cunit3'] and f[0].header['cunit3'] == 'micron':
            f[0].header['cunit3'] = 'um'
            w = wcs.WCS(f[0].header)
    new_hdr = w.to_header()

    # Remove kw added by wcslib from the new header
    for kw in wcslib_kw_to_remove:
        del new_hdr[kw]

    # Add STScI specific kw to new_hdr
    for kw in stsci_wcs_kw:
        try:
            new_hdr[kw] = f[0].header[kw]
        except KeyError:
            pass
    return new_hdr
        

def clean_primary_header(f, new_hdr):
    # f is a fits file opened in 'update' mode.
    for kw in new_hdr:
        try:
            del f[0].header[kw]
        except KeyError:
            pass


def update_sci_extension(f, new_hdr):
    # Update the SCI extension header
    f[('SCI', 1)].header.update(new_hdr)
    for kw in new_hdr:
        try:
            f[('SCI', 1)].header.comments[kw] = f[0].header.comments[kw]
        except KeyError:
            pass




if __name__ == '__main__':
    import argparse
    import os
    import six
    parser = argparse.ArgumentParser(description="Move WCS keywords from Primary to SCI extension.")
    parser.add_argument('files', help='A list of FITS filenames or a path to a directory with FITS files.')
    res = parser.parse_args()
    files = res.files

    if isinstance(files, six.string_types):
        files = os.path.abspath(files)
        if os.path.isdir(files):
            files = glob.glob("files/*.fits")
        else:
            files = glob.glob(files)

    move_wcs(files)
        

