import numpy as np
from astropy.io import fits

def table_to_hdulist(table):
    '''Create an HDUList from a set of headers and an astropy Table

    Parameters
    ----------
    table: astropy.table.Table
           The table to create the HDUList from. Note that the extra parameter, table.meta,
           if present, is expected to be a dictionary of keyword/value pairs to place into the
           primary HDU as FITS cards.

    Returns
    -------
    astropy.io.fits.HDUlist
        The HDUlist.
    '''

    # Create the Primary HDU
    fits_header = fits.Header()
    fits_header.update(table.meta)
    hdu_primary = fits.PrimaryHDU(header=fits_header)

    # Create the Table HDU
    hdu_table = fits.TableHDU(np.array(table))

    # Put it all together.
    return fits.HDUList([hdu_primary, hdu_table])
