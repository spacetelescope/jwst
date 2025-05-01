"""Make header-only images from real images."""

from pathlib import Path
import os.path as ospath

from astropy.io import fits

frompath = (
    "$DEVDIR/testdata/jwst/build7/jwstd/info/mary/b7_test/SIC_DIL/sdp.20161228/archive/level2b"
)
outpath = "./data/exposures"


def actual_path(path):
    """
    Return fully-qualified path.

    Parameters
    ----------
    path : Path or str
        Path or filename.

    Returns
    -------
    Path
        The fully qualified path.
    """
    return Path(ospath.expandvars(path)).expanduser().resolve()


def get_headers(files, outpath):
    """
    Generate header-only FITS files from FITS files.

    Parameters
    ----------
    files : list[Path or str]
        The list of files to read.
    outpath : Path
        The output path to place the header-only files
    """
    for path in files:
        name = Path(path).name
        hdul = fits.open(path)
        nhdu = fits.PrimaryHDU(header=hdul[0].header)
        nhdul = fits.HDUList([nhdu])
        nhdul.writeto(outpath / name)


if __name__ == "__main__":
    abs_frompath = actual_path(frompath)
    files = Path(abs_frompath).glob("*.fits")
    get_headers(files, actual_path(outpath))
