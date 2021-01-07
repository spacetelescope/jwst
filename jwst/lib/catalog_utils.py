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
        'jw12345_nrca_cat.ecsv'

    >>> replace_suffix_ext('jw12345_nrca_cal.fits', ['i2d'], 'cat')
        'jw12345_nrca_cal_cat.ecsv'

    >>> replace_suffix_ext('my_custom_file.fits', ['i2d'], 'cat')
        'my_custom_file_cat.ecsv'

    >>> old_suffixes = ['calints', 'crfints']
    >>> replace_suffix_ext('jw12345_nrca_calints.fits', old_suffixes, 'phot')
        'jw12345_nrca_phot.ecsv'
    >>> replace_suffix_ext('jw12345_nrca_crfints.fits', old_suffixes, 'phot')
        'jw12345_nrca_phot.ecsv'

    >>> replace_suffix_ext('jw12345_nrca_i2d.fits', ['i2d'], 'cat',
    ...                    output_dir='/jwst/my_catalogs')
        '/jwst/my_catalogs/jw12345_nrca_cat.ecsv'
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


class SkyObject(namedtuple('SkyObject', ("id",
                                         "xcentroid",
                                         "ycentroid",
                                         "sky_centroid",
                                         "isophotal_abmag",
                                         "isophotal_abmag_err",
                                         "sky_bbox_ll",
                                         "sky_bbox_lr",
                                         "sky_bbox_ul",
                                         "sky_bbox_ur",
                                         "is_star",
                                         ), rename=False)):

    """
    Sky Object container for WFSS catalog information.

    This is a convenience object for storing the catalog information
    as a named tuple. The object has explicit fields to guard for changing
    column locations in the catalog file that's read. Callers should
    validate for the minimum fields they require. This is currently populated
    for the minimum information needed by the WFSS modes in nircam and niriss.

    Parameters
    ----------
    id : int
        source identifed
    xcentroid : float
        x center of object in pixels
    ycentroid : float
        y center of object in pixels
    sky_centroid: `~astropy.coordinates.SkyCoord`
        ra and dec of the center of the object
    isophotal_abmag : float
        AB Magnitude of object
    isophotal_abmag_err : float
        Error on the AB magnitude
    sky_bbox_ll : `~astropy.coordinates.SkyCoord`
        Lower left corner of the minimum bounding box
    sky_bbox_lr : `~astropy.coordinates.SkyCoord`
        Lower right corder of the minimum bounding box
    sky_bbox_ul : `~astropy.coordinates.SkyCoord`
        Upper left corner of the minimum bounding box
    sky_bbox_ur : `~astropy.coordinates.SkyCoord`
        Upper right corner of the minimum bounding box
    is_star : bool
        Flag indicating if the object is a point source.
    """

    __slots__ = ()  # prevent instance dictionary creation for lower mem

    def __new__(cls, id=None,
                     xcentroid=None,
                     ycentroid=None,
                     sky_centroid=None,
                     isophotal_abmag=None,
                     isophotal_abmag_err=None,
                     sky_bbox_ll=None,
                     sky_bbox_lr=None,
                     sky_bbox_ul=None,
                     sky_bbox_ur=None,
                     is_star=None,):

        return super(SkyObject, cls).__new__(cls,
                                             id=id,
                                             xcentroid=xcentroid,
                                             ycentroid=ycentroid,
                                             sky_centroid=sky_centroid,
                                             isophotal_abmag=isophotal_abmag,
                                             isophotal_abmag_err=isophotal_abmag_err,
                                             sky_bbox_ll=sky_bbox_ll,
                                             sky_bbox_lr=sky_bbox_lr,
                                             sky_bbox_ul=sky_bbox_ul,
                                             sky_bbox_ur=sky_bbox_ur,
                                             is_star=is_star
                                             )

    def __str__(self):
        """Return a pretty print for the object information."""
        return ("id: {0}\n"
                "xcentroid: {1}\n"
                "ycentroid: {2}\n"
                "sky_centroid: {3}\n"
                "isophotal_abmag: {4}\n"
                "isophotal_abmag_err: {5}\n"
                "sky_bbox_ll: {6}\n"
                "sky_bbox_lr: {7}\n"
                "sky_bbox_ul: {8}\n"
                "sky_bbox_ur: {9}\n"
                "is_star: {10}"
                .format(self.id,
                        self.xcentroid,
                        self.ycentroid,
                        str(self.sky_centroid),
                        self.isophotal_abmag,
                        self.isophotal_abmag_err,
                        str(self.sky_bbox_ll),
                        str(self.sky_bbox_lr),
                        str(self.sky_bbox_ul),
                        str(self.sky_bbox_ur),
                        str(self.is_star)
                        )
                 )
