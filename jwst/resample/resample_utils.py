from __future__ import (division, print_function, unicode_literals, 
    absolute_import)

import logging
import numpy as np
import numpy.ma as ma
from scipy import interpolate

from astropy.utils.misc import isiterable
from astropy import wcs as fitswcs
from astropy.modeling.models import (Shift, Scale, Mapping, Rotation2D,
    Pix2Sky_TAN, RotateNative2Celestial, Tabular1D)
from gwcs import WCS, utils, wcstools

from .. import assign_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


DEFAULT_DOMAIN = {'lower': None,
                  'upper': None,
                  'includes_lower': True,
                  'includes_upper': False}


def make_output_wcs(wcslist):
    """ Generate output WCS here based on footprints of all input WCS objects
    Parameters
    ----------
    wcslist : list of gwcs.WCS objects

    Returns
    -------
    output_wcs : object
        WCS object, with defined domain, covering entire set of input frames

    """
    # wcslist = [i.meta.wcs for i in input_models]
    # for w, i in zip(wcslist, input_models):
    #     if w.domain is None:
    #         w.domain = create_domain(w, i.data.shape)
    output_frame = getattr(wcslist[0], 'output_frame')
    naxes = getattr(wcslist[0], output_frame).naxes
    if naxes == 3:
        output_wcs = wcs_from_spec_footprints(wcslist)
        data_size = build_size_from_spec_domain(output_wcs.domain)
    else:
        output_wcs = assign_wcs.util.wcs_from_footprints(wcslist)
        data_size = build_size_from_domain(output_wcs.domain)
    output_wcs.data_size = (data_size[1], data_size[0])
    return output_wcs

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

def build_size_from_spec_domain(domain):
    """ Return the size of the frame based on the provided domain
    """
    size = []
    for axs in domain:
        delta = axs['upper'] - axs['lower']
        #for i in [axs['includes_lower'], axs['includes_upper']]: delta += 1
        size.append(int(delta + 0.5))
    print('DEBUG size_from_domain', size)
    return tuple(size)

def calc_gwcs_pixmap(in_wcs, out_wcs):
    """ Return a pixel grid map from input frame to output frame.
    """
    if in_wcs.forward_transform.n_outputs == 3:
        # Spectral wcs with ra, dec, lamda
        grid = grid_from_spec_domain(in_wcs)
    else:
        grid = wcstools.grid_from_domain(in_wcs.domain)

    pixmap_tuple = reproject(in_wcs, out_wcs)(grid[1], grid[0])
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
    elif isinstance(wcs2, WCS):
        forward = wcs1.forward_transform
    else:
        raise ValueError("Expected astropy.wcs.WCS or gwcs.WCS object.")

    if isinstance(wcs2, fitswcs.WCS):
        args = [origin]
        inverse = wcs2.all_world2pix
    elif isinstance(wcs2, WCS):
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


def spec_domain(inwcs, domain=None):
    """
    Returns a boolean mask of the pixels where the wcs is defined.

    Build-7 workaround.

    Parameters
    ----------
    inwcs : gwcs.WCS object

    domain : optional domain

    Returns
    -------
    Boolean mask where pixels in the input that produced NaNs in the output 
    are set to True.  Valid wcs transform pixels are set to False.
    """
    if not domain:
        domain = inwcs.domain
    xstart, xend = domain[0]['lower'], domain[0]['upper']
    ystart, yend = domain[1]['lower'], domain[1]['upper']
    y, x = np.mgrid[ystart: yend, xstart: xend]
    ra, dec, lam = inwcs(x, y)
    return np.isnan(lam)


def grid_from_spec_domain(inwcs, domain=None, domain_mask=None):
    """
    Returns grid of x, y coordinates using a domain or domain mask.

    Build-7 workaround.

    Parameters
    ----------
    inwcs : gwcs.WCS object

    domain : optional domain

    domain_mask : option Boolean domain mask

    Returns
    -------
    ndarray
    Grid array if y, x inputs for the included or specified domain mask.
    """
    if not domain_mask:
        try:
            domain_mask = inwcs.domain_mask
        except AttributeError:
            domain_mask = spec_domain(inwcs)
    if not domain:
        domain = inwcs.domain
    xstart, xend = domain[0]['lower'], domain[0]['upper']
    ystart, yend = domain[1]['lower'], domain[1]['upper']
    return np.mgrid[ystart: yend, xstart: xend]
    # y, x = np.mgrid[ystart: yend, xstart: xend]
    # m = domain_mask
    # return ma.array([y, x], mask=[m, m])


def spec_footprint(wcs, domain_mask=None, domain=None):
    """
    Returns wcs footprint grid coordinates where NaNs are masked.

    Build-7 workaround.
    """
    y, x = grid_from_spec_domain(wcs, domain_mask=domain_mask)
    ra, dec, lam = wcs(x, y)
    m = np.isnan(lam)
    return ma.array([ra, dec, lam], mask=[m, m, m])


def wcs_from_spec_footprints(wcslist, refwcs=None, transform=None, domain=None):
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
    fiducial = compute_spec_fiducial(wcslist, domain=domain)
    # Create transform for output frame
    transform = compute_spec_transform(fiducial, refwcs)
    output_frame = getattr(refwcs, refwcs.output_frame)
    wnew = WCS(output_frame=output_frame, forward_transform=transform)

    # Build the domain in the output frame wcs object by running the input wcs
    # footprints through the backward transform of the output wcs
    sky = [spec_footprint(w) for w in wcslist]
    domain_grid = [wnew.backward_transform(*f) for f in sky]
    domain = []
    in_frame = getattr(refwcs, getattr(wnew, 'input_frame'))
    for axis in in_frame.axes_order:
        axis_min = np.nanmin(domain_grid[0][axis])
        axis_max = np.nanmax(domain_grid[0][axis])
        domain.append({'lower': axis_min,
            'upper': axis_max,
            'includes_lower': True, 'includes_upper': True})
    wnew.domain = domain
    return wnew


# def compute_spec_transform(fiducial, refwcs):
#     """
#     Compute a simple transform given a fidicial point in a spatial-spectral wcs.
#     """
#     cdelt1 = refwcs.wcsinfo.cdelt1 / 3600.
#     cdelt2 = refwcs.wcsinfo.cdelt2 / 3600.
#     cdelt3 = refwcs.wcsinfo.cdelt3
#     roll_ref = refwcs.wcsinfo.roll_ref

#     y, x = grid_from_spec_domain(w)       
#     ra, dec, lam = refwcs(x, y)

#     min_lam = np.nanmin(lam)

#     offset = Shift(0.) & Shift(0.)
#     rot = Rotation2D(roll_ref)
#     scale = Scale(cdelt1) & Scale(cdelt2)
#     tan = Pix2Sky_TAN()
#     skyrot = RotateNative2Celestial(fiducial[0][0], fiducial[0][1], 180.0)
#     spatial = offset | rot | scale | tan | skyrot
#     spectral = Scale(cdelt3) | Shift(min_lam)
#     mapping = Mapping((1, 1, 0),)
#     mapping.inverse = Mapping((2, 1))
#     transform = mapping | spatial & spectral
#     transform.outputs = ('ra', 'dec', 'lamda')
#     return transform


def compute_spec_transform(fiducial, refwcs):
    """
    Compute a simple transform given a fidicial point in a spatial-spectral wcs.
    """
    cdelt3 = refwcs.wcsinfo.cdelt3
    
    # Build tabular lookup table model for the slit spatial transform
    w = refwcs
    y, x = grid_from_spec_domain(w)       
    ra, dec, lam = refwcs(x, y)
    # Find dispersion axis
    disp_axis = int(ra.shape[0] > ra.shape[1])
    xdisp_axis = int(ra.shape[0] < ra.shape[1])
    # Map dispersion axis to the x or y grids
    if disp_axis == 0:
        grid_disp = x
        grid_xdisp = y
    else:
        grid_disp = y
        grid_xdisp = x
    # Find middle of dispersion by looking at the domain
    mid_disp = int((w.domain[disp_axis]['upper'] - w.domain[disp_axis]['lower'])/2)
    mid_xdisp = int((w.domain[xdisp_axis]['upper'] - w.domain[xdisp_axis]['lower'])/2)
    # Make a cutout 4 pix in the dispersion range, of the full slit, trimming nans
    lam_subarr = lam[:,mid_disp-2:mid_disp+2]
    # lam_subarr = lam_subarr[~np.isnan(lam_subarr)]
    # lam_subarr = lam_subarr.reshape(-1, 4)
    # And use this midpoint to walk up the slit from bottom to top
    # interpolating to find the input x,y that gives constant output lambda
    fiducial_lambda = lam[mid_xdisp,mid_disp]
    xpos = []
    for row in lam_subarr:
        if any(np.isnan(row)):
            xpos.append(np.nan)
        else:
            f = interpolate.interp1d(row, x[0,mid_disp-2:mid_disp+2])
            xpos.append(f(fiducial_lambda))
    xpos = np.array(xpos)

    # Now compute RA/DEC values for these new X positions; remove NaNs
    ra_slit, dec_slit, lam_slit = w(xpos,y[:,mid_disp])
    ra_slit, dec_slit = ra_slit[~np.isnan(ra_slit)], dec_slit[~np.isnan(dec_slit)]

    # Build the tabular models now
    points = (np.arange(len(ra_slit)).tolist())
    ra_tabular = Tabular1D(points, lookup_table=ra_slit, bounds_error=False, fill_value=None)
    dec_tabular = Tabular1D(points, lookup_table=dec_slit, bounds_error=False, fill_value=None)
    ra_tabular.inverse = Tabular1D(ra_slit, lookup_table=points, bounds_error=False, fill_value=None)
    dec_tabular.inverse = Tabular1D(np.flipud(dec_slit), lookup_table=np.flipud(points), bounds_error=False, fill_value=None)

    # radec_tabular = Tabular2D(lookup_table='', bounds_error=False, fill_value=None)
    # radec_tabular.inverse = Tabular2D(np.array(lookup_table='', bounds_error=False, fill_value=None)

    min_lam = np.nanmin(lam)

    spatial = ra_tabular & dec_tabular
    spectral = Scale(cdelt3) | Shift(min_lam)
    mapping = Mapping((1, 1, 0),)
    mapping.inverse = Mapping((2, 1))
    transform = mapping | spatial & spectral
    transform.outputs = ('ra', 'dec', 'lamda')
    return transform


def compute_spec_fiducial(wcslist, domain=None):
    """
    For a celestial footprint this is the center.
    For a spectral footprint, it is the beginning of the range.

    This function assumes all WCSs have the same output coordinate frame.

    Build-7 workaround.
    """
    output_frame = getattr(wcslist[0], 'output_frame')
    axes_types = getattr(wcslist[0], output_frame).axes_type
    spatial_axes = np.array(axes_types) == 'SPATIAL'
    spectral_axes = np.array(axes_types) == 'SPECTRAL'
    footprints = ma.hstack([spec_footprint(w,
        domain=domain) for w in wcslist])
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
    #    c = coord.SkyCoord(lon_fiducial, lat_fiducial, unit='deg')
    if (spectral_footprint).any():
        fiducial[spectral_axes] = spectral_footprint.min()
    return ((fiducial[spatial_axes]), fiducial[spectral_axes])
    #return (c, spectral_footprint.min())


