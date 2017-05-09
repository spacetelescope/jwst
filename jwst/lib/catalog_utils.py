"""
Utilities for naming source catalogs.
"""

import re
from os.path import split, splitext, join, abspath, expanduser
from collections import namedtuple


def replace_suffix_ext(filename, old_suffix_list, new_suffix,
                       output_ext='ecsv', output_dir=None):
    """
    Replace the suffix and extension of a filename.

    If the last suffix in the input filename is in the
    ``old_suffix_list`` then it is replaced by the ``new_suffix``.
    Otherwise, the ``new_suffix`` is appended to the input filename.

    Parameters
    ----------
    filename : str
        The filename to modify.

    old_suffix_list : list of str
        The list of filename suffixes that will be replaced.

    new_suffix : str
        The new filename suffix.

    output_ext : str, optional
        The extension of the output filename.  The default is 'ecsv'.

    output_dir : str or `None`, optional
        The output directory name.  If `None` then the current directory
        will be used.

    Returns
    -------
    result : str
        The modified filename.

    Examples
    --------
    >>> from jwst.lib.catalog_utils import replace_suffix_ext
    >>> replace_suffix_ext('jw12345_nrca_i2d.fits', ['i2d'], 'cat')
    >>> 'jw12345_nrca_cat.ecsv'

    >>> replace_suffix_ext('jw12345_nrca_cal.fits', ['i2d'], 'cat')
    >>> 'jw12345_nrca_cal_cat.ecsv'

    >>> replace_suffix_ext('my_custom_file.fits', ['i2d'], 'cat')
    'my_custom_file_cat.ecsv'

    >>> old_suffixes = ['calints', 'crfints']
    >>> replace_suffix_ext('jw12345_nrca_calints.fits', old_suffixes, 'phot')
    >>> 'jw12345_nrca_phot.ecsv'
    >>> replace_suffix_ext('jw12345_nrca_crfints.fits', old_suffixes, 'phot')
    >>> 'jw12345_nrca_phot.ecsv'

    >>> replace_suffix_ext('jw12345_nrca_i2d.fits', ['i2d'], 'cat',
    ...                    output_dir='/jwst/my_catalogs')
    >>> '/jwst/my_catalogs/jw12345_nrca_cat.ecsv'
    """

    path, filename = split(filename)
    name, ext = splitext(filename)
    remove_suffix = '^(.+?)(_(' + '|'.join(old_suffix_list) + '))?$'
    match = re.match(remove_suffix, name)
    name = match.group(1)

    output_path = '{0}_{1}.{2}'.format(name, new_suffix, output_ext)
    if output_dir is not None:
        output_path = abspath(expanduser(join(output_dir, output_path)))

    return output_path



class SkyObject(namedtuple('SkyObject', ("sid",
                                         "xcentroid",
                                         "ycentroid",
                                         "ra_icrs_centroid",
                                         "dec_icrs_centroid",
                                         "abmag",
                                         "abmag_error",
                                         "ramin",
                                         "decmin",
                                         "ramax",
                                         "decmax"), rename=False)):

    """ Sky Object

    This is a convenience object for storing the catalog information
    as a named tuple. The object has explicit fields to guard for changing column
    locations in the catalog file that's read. Callers should
    validate for the minimum fields they require. This is currently populated for
    the minimum information needed by the WFSS modes in nircam and niriss.
    """

    __slots__ = ()  # prevent instance dictionary creation for low mem

    def __new__(cls, sid=None,
                     xcentroid=None,
                     ycentroid=None,
                     ra_icrs_centroid=None,
                     dec_icrs_centroid=None,
                     abmag=None,
                     abmag_error=None,
                     ramin=None,
                     decmin=None,
                     ramax=None,
                     decmax=None):

        return super(SkyObject, cls).__new__(cls,
                                             sid,
                                             xcentroid,
                                             ycentroid,
                                             ra_icrs_centroid,
                                             dec_icrs_centroid,
                                             abmag,
                                             abmag_error,
                                             ramin,
                                             decmin,
                                             ramax,
                                             decmax)

    def __str__(self):
        """Return a pretty print for the object information."""
        return ("id: {0}\n"
                "xcentroid: {1}\n"
                "ycentroid: {2}\n"
                "ra_icrs_centroid: {3}\n"
                "dec_icrs_centroid: {4}\n"
                "abmag: {5}\n"
                "abmag_error: {6}\n"
                "ramin: {7}\n"
                "decmin: {8}\n"
                "ramax: {9}\n"
                "decmax: {10}\n"
                .format(self.sid,
                        self.xcentroid,
                        self.ycentroid,
                        self.ra_icrs_centroid,
                        self.dec_icrs_centroid,
                        self.abmag,
                        self.abmag_error,
                        self.ramin,
                        self.decmin,
                        self.ramax,
                        self.decmax))
