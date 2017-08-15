"""
Utilities for naming source catalogs.
"""

import re
from os.path import split, splitext, join, abspath, expanduser


def replace_suffix_ext(filename, old_suffix_list, new_suffix,
                       output_ext='escv', output_dir=None):
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
        The extension of the output filename.  The default is 'escv'.

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
    >>> 'jw12345_nrca_cat.escv'

    >>> replace_suffix_ext('jw12345_nrca_cal.fits', ['i2d'], 'cat')
    >>> 'jw12345_nrca_cal_cat.escv'

    >>> replace_suffix_ext('my_custom_file.fits', ['i2d'], 'cat')
    'my_custom_file_cat.escv'

    >>> old_suffixes = ['calints', 'crfints']
    >>> replace_suffix_ext('jw12345_nrca_calints.fits', old_suffixes, 'phot')
    >>> 'jw12345_nrca_phot.escv'
    >>> replace_suffix_ext('jw12345_nrca_crfints.fits', old_suffixes, 'phot')
    >>> 'jw12345_nrca_phot.escv'

    >>> replace_suffix_ext('jw12345_nrca_i2d.fits', ['i2d'], 'cat',
    ...                    output_dir='/jwst/my_catalogs')
    >>> '/jwst/my_catalogs/jw12345_nrca_cat.escv'
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
