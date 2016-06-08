from .csv_to_table import csv_to_table
from .table_to_hdulist import table_to_hdulist

def csv_to_hdulist(handle, **kwargs):
    '''Convert and CSV file, with keyword/value pairs in leading comments, to an HDUList

    Parameters
    ----------
    handle: File handle
            The file object being read from.

    kwargs: dict
            Named parameters to pass to csv_to_table.

    Returns
    -------
    astropy.io.fits.HDUList
        The HDUList containing the CSV file info.
    '''
    return table_to_hdulist(csv_to_table(handle, **kwargs))
