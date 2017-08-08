from __future__ import (division, print_function, unicode_literals,
    absolute_import)

import numpy as np
import numpy.ma as ma
from scipy import interpolate

from astropy.utils.misc import isiterable
from astropy import wcs as fitswcs
from astropy.coordinates import SkyCoord
from astropy.modeling.models import (Shift, Scale, Mapping, Rotation2D,
    Pix2Sky_TAN, RotateNative2Celestial, Tabular1D, AffineTransformation2D)
from gwcs import WCS, utils, wcstools


from .. import assign_wcs

import logging
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
            w.bounding_box = bounding_box_from_shape(i.data.shape)

    output_frame = wcslist[0].output_frame
    naxes = wcslist[0].output_frame.naxes

    if naxes == 3:
        # THIS BLOCK CURRENTLY ISN"T USED BY resample_spec
        output_wcs = wcs_from_spec_footprints(wcslist)
        data_size = shape_from_bounding_box(output_wcs.bounding_box)
    elif naxes == 2:
        output_wcs = assign_wcs.util.wcs_from_footprints(input_models)
        data_size = shape_from_bounding_box(output_wcs.bounding_box)

    output_wcs.data_size = (data_size[1], data_size[0])
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
    #coeff_x = np.mean([np.abs(offset_x[1].value), np.abs(offset_y[0].value)])
    #coeff_y = np.mean([np.abs(offset_x[0].value), np.abs(offset_y[1].value)])

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


def bounding_box_from_shape(shape):
    """ Create bounding_box for WCS based on shape of model data.
    """
    bb = []
    for s in reversed(shape):
        bb.append((-0.5, s - 0.5))
    return tuple(bb)


def shape_from_bounding_box(bounding_box):
    """ Return the size of the frame based on the provided domain
    """
    size = []
    for axs in bounding_box:
        delta = axs[1] - axs[0]
        #for i in [axs['includes_lower'], axs['includes_upper']]: delta += 1
        size.append(int(delta + 0.5))
    return tuple(reversed(size))


def calc_gwcs_pixmap(in_wcs, out_wcs, shape=None):
    """ Return a pixel grid map from input frame to output frame.
    """
    if shape:
        bb = bounding_box_from_shape(shape)
        log.debug("Bounding box from data shape: {}".format(bb))
    else:
        bb = in_wcs.bounding_box
        log.debug("Bounding box from WCS: {}".format(in_wcs.bounding_box))

    grid = wcstools.grid_from_bounding_box(bb, step=(1, 1), center=True)
    pixmap = np.dstack(reproject(in_wcs, out_wcs)(grid[0], grid[1]))
    return pixmap


def reproject(wcs1, wcs2):
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
        forward_transform = wcs1.all_pix2world
    elif isinstance(wcs2, WCS):
        forward_transform = wcs1.forward_transform
    else:
        raise ValueError("Expected astropy.wcs.WCS or gwcs.WCS object.")

    if isinstance(wcs2, fitswcs.WCS):
        backward_transform = wcs2.all_world2pix
    elif isinstance(wcs2, WCS):
        #inverse = wcs2.forward_transform.inverse
        backward_transform = wcs2.backward_transform
    else:
        raise ValueError("Expected astropy.wcs.WCS or gwcs.WCS object.")

    def _reproject(x, y):
        sky = forward_transform(x, y)
        flat_sky = []
        for axis in sky:
            flat_sky.append(axis.flatten())
        det = backward_transform(*tuple(flat_sky))
        det_reshaped = []
        for axis in det:
            det_reshaped.append(axis.reshape(x.shape))
        return tuple(det_reshaped)
    return _reproject

# Following function can be deprecated
def spec_bounding_box(inwcs, bounding_box=None):
    """
    Returns a boolean mask of the pixels where the wcs is defined.

    Build-7 workaround.

    Parameters
    ----------
    inwcs : gwcs.WCS object

    bounding_box : optional bounding_box

    Returns
    -------
    Boolean mask where pixels in the input that produced NaNs in the output
    are set to True.  Valid wcs transform pixels are set to False.
    """
    if not bounding_box:
        bounding_box = inwcs.bounding_box
    x, y = wcstools.grid_from_bounding_box(in_wcs.bounding_box, step=(1, 1),
        center=True)
    ra, dec, lam = inwcs(x, y)
    return np.isnan(lam)

# Following function can be deprecated
def grid_from_spec_bounding_box(inwcs, bounding_box=None, mask=None):
    """
    Returns grid of x, y coordinates using a bounding_box or mask.

    Build-7 workaround.

    Parameters
    ----------
    inwcs : gwcs.WCS object

    bounding_box : optional bounding_box

    mask : option Boolean bounding_box mask

    Returns
    -------
    ndarray
    Grid array if y, x inputs for the included or specified domain mask.
    """
    if not mask:
        try:
            mask = inwcs.bounding_box_mask
        except AttributeError:
            mask = spec_bounding_box(inwcs)
    if not bounding_box:
        bounding_box = inwcs.bounding_box
    x, y = wcstools.grid_from_bounding_box(in_wcs.bounding_box, step=(1, 1),
        center=True)
    return x, y


def spec_footprint(in_wcs, bounding_box=None, mask=None):
    """
    Returns wcs footprint grid coordinates where NaNs are masked.

    Build-7 workaround.
    """
    x, y = wcstools.grid_from_bounding_box(in_wcs.bounding_box, step=(1, 1),
        center=True)
    ra, dec, lam = in_wcs(x, y)
    m = np.isnan(lam)
    return ma.array([ra, dec, lam], mask=[m, m, m])


def wcs_from_spec_footprints(wcslist, refwcs=None, transform=None,
    bounding_box=None):
    """
    Create a WCS from a list of spatial/spectral WCS.

    Build-7 workaround.
    """
    if not isiterable(wcslist):
        raise ValueError("Expected 'wcslist' to be an iterable of gwcs.WCS")
    if not all([isinstance(w, WCS) for w in wcslist]):
        raise TypeError("All items in 'wcslist' must have instance of gwcs.WCS")
    if refwcs is None:
        refwcs = wcslist[0]
    else:
        if not isinstance(refwcs, WCS):
            raise TypeError("Expected refwcs to be an instance of gwcs.WCS.")

    # TODO: generalize an approach to do this for more than one wcs.  For
    # now, we just do it for one, using the api for a list of wcs.
    # Compute a fiducial point for the output frame at center of input data
    fiducial = compute_spec_fiducial(wcslist, bounding_box=bounding_box)
    # Create transform for output frame
    transform = compute_spec_transform(fiducial, refwcs)
    output_frame = refwcs.output_frame
    wnew = WCS(output_frame=output_frame, forward_transform=transform)

    # Build the bounding_box in the output frame wcs object by running the
    # input wcs footprints through the backward transform of the output wcs
    sky = [spec_footprint(w) for w in wcslist]
    bounding_box_grid = [wnew.backward_transform(*f) for f in sky]

    sky0 = sky[0]
    det = bounding_box_grid[0]
    offsets = []
    input_frame = refwcs.input_frame
    for axis in input_frame.axes_order:
        axis_min = np.nanmin(det[axis])
        offsets.append(axis_min)
    transform = Shift(offsets[0]) & Shift(offsets[1]) | transform
    wnew = WCS(output_frame=output_frame, input_frame=input_frame,
        forward_transform=transform)

    # Build the bounding_box in the output frame wcs object
    bounding_box = []
    for axis in input_frame.axes_order:
        axis_min = np.nanmin(bounding_box_grid[axis])
        axis_max = np.nanmax(bounding_box_grid[axis])
        bounding_box.append((axis_min, axis_max))
    wnew.bounding_box = tuple(bounding_box)
    return wnew


def compute_spec_transform(fiducial, refwcs):
    """
    Compute a simple transform given a fidicial point in a spatial-spectral wcs.
    """
    cdelt1 = refwcs.wcsinfo.cdelt1 / 3600.
    cdelt2 = refwcs.wcsinfo.cdelt2 / 3600.
    cdelt3 = refwcs.wcsinfo.cdelt3
    roll_ref = refwcs.wcsinfo.roll_ref

    y, x = grid_from_spec_domain(refwcs)
    ra, dec, lam = refwcs(x, y)

    min_lam = np.nanmin(lam)

    offset = Shift(0.) & Shift(0.)
    rot = Rotation2D(roll_ref)
    scale = Scale(cdelt1) & Scale(cdelt2)
    tan = Pix2Sky_TAN()
    skyrot = RotateNative2Celestial(fiducial[0][0], fiducial[0][1], 180.0)
    spatial = offset | rot | scale | tan | skyrot
    spectral = Scale(cdelt3) | Shift(min_lam)
    mapping = Mapping((1, 1, 0),)
    mapping.inverse = Mapping((2, 1))
    transform = mapping | spatial & spectral
    transform.outputs = ('ra', 'dec', 'lamda')
    return transform


def compute_spec_fiducial(wcslist):
    """
    For a celestial footprint this is the center.
    For a spectral footprint, it is the beginning of the range.

    This function assumes all WCSs have the same output coordinate frame.

    Build-7 workaround.
    """
    output_frame = wcslist[0].output_frame
    axes_types = wcslist[0].output_frame.axes_type
    spatial_axes = np.array(axes_types) == 'SPATIAL'
    spectral_axes = np.array(axes_types) == 'SPECTRAL'
    footprints = ma.hstack([spec_footprint(w,
        bounding_box=w.bounding_box) for w in wcslist])
    spatial_footprint = footprints[spatial_axes]
    spectral_footprint = footprints[spectral_axes]
    # Compute center of footprint
    fiducial = np.empty(len(axes_types))
    if (spatial_footprint).any():
        lon, lat = spatial_footprint
        lon, lat = np.deg2rad(lon), np.deg2rad(lat)
        x_mean = np.mean(np.cos(lat) * np.cos(lon))
        y_mean = np.mean(np.cos(lat) * np.sin(lon))
        z_mean = np.mean(np.sin(lat))
        lon_fiducial = np.rad2deg(np.arctan2(y_mean, x_mean)) % 360.0
        lat_fiducial = np.rad2deg(np.arctan2(z_mean, np.sqrt(x_mean ** 2 +
            y_mean ** 2)))
        fiducial[spatial_axes] = lon_fiducial, lat_fiducial
    if (spectral_footprint).any():
        fiducial[spectral_axes] = spectral_footprint.min()
    return ((fiducial[spatial_axes]), fiducial[spectral_axes])
