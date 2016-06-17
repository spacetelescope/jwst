""" General utilities used throughout the resample package.

"""
import numpy as np

from astropy import wcs as fitswcs
from gwcs import wcs
from gwcs import wcstools

from jwst import assign_wcs

DEFAULT_DOMAIN = {'lower': None, 'upper': None, 'includes_lower': True, 'includes_upper': False}

def make_output_wcs(input_models):
    """ Generate output WCS here based on footprints of all input WCS objects
    Parameters
    ----------
    input_models : ModelContainer or list
        ModelContainer or list of ImageModels for all input images to be combined

    Returns
    -------
    output_wcs : object
        WCS object, with defined domain, covering entire set of input frames

    """
    wcslist = [i.meta.wcs for i in input_models]
    for w, i in zip(wcslist, input_models):
        if w.domain is None:
            w.domain = create_domain(w, i.data.shape)

    output_wcs = assign_wcs.util.wcs_from_footprints(wcslist)
    data_size = build_size_from_domain(output_wcs.domain)
    output_wcs.data_size = (data_size[1], data_size[0])
    return output_wcs

def define_wcslin(model):
    """ Create an undistorted version of the input WCS.
    """
    from jwst.assign_wcs import jwcs
    hdr = model.storage.get_fits_header('PRIMARY')
    return jwcs.JWSTWCS(model.meta, hdr)

def create_domain(wcs, shape):
    """ Create domain for WCS based on shape of model data.
    """
    wcs_domain = []
    for s in shape:
        domain = DEFAULT_DOMAIN.copy()
        domain['lower'] = 0
        domain['upper'] = s
        wcs_domain.append(domain)
    return wcs_domain

def build_size_from_domain(domain):
    """ Return the size of the frame based on the provided domain
    """
    size = []
    for axs in domain:
        delta = axs['upper'] - axs['lower']
        #for i in [axs['includes_lower'], axs['includes_upper']]: delta += 1
        size.append(int(delta + 0.5))
    return tuple(size)

def calc_gwcs_pixmap(in_wcs, out_wcs):

    g = wcstools.grid_from_domain(in_wcs.domain)

    pixmap_tuple = reproject(in_wcs, out_wcs)(g[1], g[0])
    pixmap = np.dstack(pixmap_tuple)
    return pixmap

def reproject(wcs1, wcs2, origin=0):
    """
    Given two WCSs return a function which takes pixel coordinates in
    the first WCS and computes them in the second one.
    It performs the forward transformation of ``wcs1`` followed by the
    inverse of ``wcs2``.

    NOTE
    ----
    This function should be removed when it can be found in another
    accessible package for JWST pipeline use.

    Parameters
    ----------
    wcs1, wcs2 : `~astropy.wcs.WCS` or `~gwcs.wcs.WCS`
        WCS objects.
    origin : {0, 1}
        Whether to use 0- or 1-based pixel coordinates.
    Returns
    -------
    _reproject : func
        Function to compute the transformations.  It takes x, y
        positions in ``wcs1`` and returns x, y positions in ``wcs2``.
    """

    args = []
    if isinstance(wcs1, fitswcs.WCS):
        forward = wcs1.all_pix2world
        args = [origin]
    elif isinstance(wcs2, wcs.WCS):
        forward = wcs1.forward_transform
    else:
        raise ValueError("Expected astropy.wcs.WCS or gwcs.WCS object.")

    if isinstance(wcs2, fitswcs.WCS):
        args = [origin]
        inverse = wcs2.all_world2pix
    elif isinstance(wcs2, wcs.WCS):
        #inverse = wcs2.forward_transform.inverse
        inverse = wcs2.backward_transform
    else:
        raise ValueError("Expected astropy.wcs.WCS or gwcs.WCS object.")

    def _reproject(x, y):
        forward_args = [x, y] + args
        sky = forward(*forward_args)
        inverse_args = list(sky) + args
        return inverse(*inverse_args)
    return _reproject
