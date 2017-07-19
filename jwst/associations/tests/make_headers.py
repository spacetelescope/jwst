"""
Make header-only images from real images
"""
from glob import glob
import os.path as ospath

from astropy.io import fits

frompath = '$DEVDIR/testdata/jwst/build7/jwstd/info/mary/b7_test/SIC_DIL/sdp.20161228/archive/level2b'
outpath = './data/exposures'


def actual_path(path):
    """Returns fully-qualified path"""
    return ospath.abspath(ospath.expandvars(ospath.expanduser(path)))


def get_headers(files, outpath):
    """Make header-only FITS files from FITS files

    Parameters
    ----------
    files: [str[, ...]]
        The files to read

    outpath: str
        The output path to place the header-only files
    """
    for path in files:
        name = ospath.basename(path)
        hdul = fits.open(path)
        nhdu = fits.PrimaryHDU(header=hdul[0].header)
        nhdul = fits.HDUList([nhdu])
        nhdul.writeto(ospath.join(outpath, name))


if __name__ == '__main__':
    abs_frompath = actual_path(frompath)
    files = glob(ospath.join(abs_frompath, '*.fits'))
    get_headers(files, actual_path(outpath))
