import logging
import warnings

import numpy as np
from astropy import wcs as fitswcs
from astropy.coordinates import SkyCoord
from astropy.modeling.models import Scale, AffineTransformation2D
from astropy.modeling import Model
from gwcs import WCS, wcstools
from astropy.nddata.bitmask import interpret_bit_flags

from ..assign_wcs.util import wcs_from_footprints, wcs_bbox_from_shape

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def make_output_wcs(input_models):
    """ Generate output WCS here based on footprints of all input WCS objects
    Parameters
    ----------
    wcslist : list of gwcs.WCS objects

    Returns
    -------
    output_wcs : object
        WCS object, with defined domain, covering entire set of input frames

    """

    # The API needing input_models instead of just wcslist is because
    # currently the domain is not defined in any of imaging modes for NIRCam
    # NIRISS or MIRI
    #
    # TODO: change the API to take wcslist instead of input_models and
    #       remove the following block
    wcslist = [i.meta.wcs for i in input_models]
    for w, i in zip(wcslist, input_models):
        if w.bounding_box is None:
            w.bounding_box = wcs_bbox_from_shape(i.data.shape)
    naxes = wcslist[0].output_frame.naxes

    if naxes == 3:
        # THIS BLOCK CURRENTLY ISN"T USED BY resample_spec
        pass
    elif naxes == 2:
        output_wcs = wcs_from_footprints(input_models)
        output_wcs.data_size = shape_from_bounding_box(output_wcs.bounding_box)

    # Check that the output data shape has no zero length dimensions
    if not np.product(output_wcs.data_size):
        raise ValueError("Invalid output frame shape: "
                         "{}".format(output_wcs.data_size))

    return output_wcs


def compute_output_transform(refwcs, filename, fiducial):
    """Compute a simple FITS-type WCS transform
    """
    x0, y0 = refwcs.backward_transform(*fiducial)
    x1 = x0 + 1
    y1 = y0 + 1
    ra0, dec0 = refwcs(x0, y0)
    ra_xdir, dec_xdir = refwcs(x1, y0)
    ra_ydir, dec_ydir = refwcs(x0, y1)

    position0 = SkyCoord(ra=ra0, dec=dec0, unit='deg')
    position_xdir = SkyCoord(ra=ra_xdir, dec=dec_xdir, unit='deg')
    position_ydir = SkyCoord(ra=ra_ydir, dec=dec_ydir, unit='deg')
    offset_xdir = position0.spherical_offsets_to(position_xdir)
    offset_ydir = position0.spherical_offsets_to(position_ydir)

    xscale = np.abs(position0.separation(position_xdir).value)
    yscale = np.abs(position0.separation(position_ydir).value)
    scale = np.sqrt(xscale * yscale)

    c00 = offset_xdir[0].value / scale
    c01 = offset_xdir[1].value / scale
    c10 = offset_ydir[0].value / scale
    c11 = offset_ydir[1].value / scale
    pc_matrix = AffineTransformation2D(matrix=[[c00, c01], [c10, c11]])
    cdelt = Scale(scale) & Scale(scale)

    return pc_matrix | cdelt


def shape_from_bounding_box(bounding_box):
    """ Return a numpy shape based on the provided bounding_box
    """
    size = []
    for axs in bounding_box:
        delta = axs[1] - axs[0]
        size.append(int(delta + 0.5))
    return tuple(reversed(size))


def calc_gwcs_pixmap(in_wcs, out_wcs, shape=None):
    """ Return a pixel grid map from input frame to output frame.
    """
    if shape:
        bb = wcs_bbox_from_shape(shape)
        log.debug("Bounding box from data shape: {}".format(bb))
    else:
        bb = in_wcs.bounding_box
        log.debug("Bounding box from WCS: {}".format(in_wcs.bounding_box))

    grid = wcstools.grid_from_bounding_box(bb)
    pixmap = np.dstack(reproject(in_wcs, out_wcs)(grid[0], grid[1]))
    pixmap[np.isnan(pixmap)] = -1

    return pixmap


def reproject(wcs1, wcs2):
    """
    Given two WCSs or transforms return a function which takes pixel
    coordinates in the first WCS or transform and computes them in the second
    one. It performs the forward transformation of ``wcs1`` followed by the
    inverse of ``wcs2``.

    Parameters
    ----------
    wcs1, wcs2 : `~astropy.wcs.WCS` or `~gwcs.wcs.WCS` or `~astropy.modeling.Model`
        WCS objects.

    Returns
    -------
    _reproject : func
        Function to compute the transformations.  It takes x, y
        positions in ``wcs1`` and returns x, y positions in ``wcs2``.
    """

    if isinstance(wcs1, fitswcs.WCS):
        forward_transform = wcs1.all_pix2world
    elif isinstance(wcs1, WCS):
        forward_transform = wcs1.forward_transform
    elif issubclass(wcs1, Model):
        forward_transform = wcs1
    else:
        raise TypeError("Expected input to be astropy.wcs.WCS or gwcs.WCS "
            "object or astropy.modeling.Model subclass")

    if isinstance(wcs2, fitswcs.WCS):
        backward_transform = wcs2.all_world2pix
    elif isinstance(wcs2, WCS):
        backward_transform = wcs2.backward_transform
    elif issubclass(wcs2, Model):
        backward_transform = wcs2.inverse
    else:
        raise TypeError("Expected input to be astropy.wcs.WCS or gwcs.WCS "
            "object or astropy.modeling.Model subclass")

    def _reproject(x, y):
        sky = forward_transform(x, y)
        flat_sky = []
        for axis in sky:
            flat_sky.append(axis.flatten())
        # Filter out RuntimeWarnings due to computed NaNs in the WCS
        warnings.simplefilter("ignore")
        det = backward_transform(*tuple(flat_sky))
        warnings.resetwarnings()
        det_reshaped = []
        for axis in det:
            det_reshaped.append(axis.reshape(x.shape))
        return tuple(det_reshaped)
    return _reproject


def build_driz_weight(model, weight_type=None, good_bits=None):
    """ Create input weighting image
    """
    dqmask = build_mask(model.dq, good_bits)
    exptime = model.meta.exposure.exposure_time

    if weight_type == 'error':
        err_model = np.nan_to_num(model.err)
        inwht = (exptime / err_model)**2 * dqmask
        log.debug("DEBUG weight mask: {} {}".format(type(inwht), np.sum(inwht)))
    # elif weight_type == 'ivm':
    #     _inwht = img.buildIVMmask(chip._chip,dqarr,pix_ratio)
    elif weight_type == 'exptime':
        inwht = exptime * dqmask
    else:
        inwht = np.ones(model.data.shape, dtype=model.data.dtype)
    return inwht


def build_mask(dqarr, bitvalue):
    """ Builds a bit-mask from an input DQ array and a bitvalue flag
    """
    bitvalue = interpret_bit_flags(bitvalue)

    if bitvalue is None:
        return (np.ones(dqarr.shape, dtype=np.uint8))
    return np.logical_not(np.bitwise_and(dqarr, ~bitvalue)).astype(np.uint8)
