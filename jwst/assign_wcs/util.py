"""
Utility function for WCS

"""

import logging
import functools
import numpy as np
import numpy.ma as ma

from astropy.utils.misc import isiterable
from astropy.io import fits
from astropy.modeling import projections
from astropy.modeling import models as astmodels
from astropy.modeling.core import Model
from astropy import coordinates as coord
import astropy.units as u

from gwcs import WCS, utils, coordinate_frames
from gwcs.wcstools import frame2transform

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def reproject(wcs1, wcs2, origin=0):
    """
    Given two WCSs return a function which takes pixel coordinates in
    the first WCS and computes their location in the second one.

    It performs the forward transformation of ``wcs1`` followed by the
    inverse of ``wcs2``.

    Parameters
    ----------
    wcs1, wcs2 : `~gwcs.wcs.WCS`
        WCS objects.

    Returns
    -------
    _reproject : func
        Function to compute the transformations.  It takes x, y
        positions in ``wcs1`` and returns x, y positions in ``wcs2``.
    """

    def _reproject(x, y):
        sky = wcs1.forward_transform(x, y)
        return wcs2.backward_transform(*sky)
    return _reproject


def wcs_from_footprints(wcslist, refwcs=None, transform=None, domain=None):
    """
    Create a WCS from a list of WCS objects.

    A fiducial point in the output coordinate frame is created from  the
    footprints of all WCS objects. For a spatial frame this is the center
    of the union of the footprints. For a spectral frame the fiducial is in
    the beginning of the footprint range.
    If ``refwcs`` is not specified, the first WCS object in the list is considered
    a reference. The output coordinate frame and projection (for celestial frames)
    is taken from ``refwcs``.
    If ``transform`` is not suplied, a compound transform comprised of
    scaling and rotation is copied from ``refwcs``.
    If ``domain`` is not supplied, the domain of the new WCS is computed
    from the domains of all input WCSs

    Parameters
    ----------
    wcslist : list of `~gwcs.wcs.WCS`
        A list of WCS objects.
    refwcs : `~gwcs.wcs.WCS`, optional
        Reference WCS. The output coordinate frame, the projection and a
        scaling and rotation transform is created from it. If not supplied
        the first WCS in the list is used as ``refwcs``.
    transform : `~astropy.modeling.core.Model`, optional
        A transform, passed to :class_method:`~gwcs.WCS.wcs_from_fiducial`
        If not supplied Scaling | Rotation is computed from ``refwcs``.
    domain : list of dicts, optional
        Domain of the new WCS.
        If not supplied it is computed from the domain of all inputs.
    """
    if not isiterable(wcslist):
        raise ValueError("Expected 'wcslist' to be an iterable of WCS objects.")
    if not all([isinstance(w, WCS) for w in wcslist]):
        raise TypeError("All items in wcslist are to be instances of gwcs.WCS.")
    if refwcs is None:
        refwcs = wcslist[0]
    else:
        if not isinstance(refwcs, WCS):
            raise TypeError("Expected refwcs to be an instance of gwcs.WCS.")

    fiducial = compute_fiducial(wcslist, domain)
    prj = np.array([isinstance(m, projections.Projection) for m \
                    in refwcs.forward_transform]).nonzero()[0]
    if prj:
        prj = refwcs.forward_transform[prj[0]]
    else:
        prj = astmodels.Pix2Sky_TAN()
    trans = []
    scales = [m for m in refwcs.forward_transform if isinstance(m, astmodels.Scale)]
    if scales:
        trans.append(functools.reduce(lambda x, y: x & y, scales))
    rotation = [m for m in refwcs.forward_transform if \
                isinstance(m, astmodels.AffineTransformation2D)]
    if rotation:
        trans.append(rotation[0])
    if trans:
        tr = functools.reduce(lambda x, y: x | y, trans)
    else:
        tr = None
    out_frame = refwcs.output_frame
    wnew = wcs_from_fiducial(fiducial, coordinate_frame=out_frame,
                             projection=prj, transform=tr)

    #domain_bounds = np.hstack([gwutils._domain_to_bounds(d) for d in \
    #[w.domain for w in wcslist]])
    domain_footprints = [w.footprint() for w in wcslist]
    domain_bounds = np.hstack([wnew.backward_transform(*f) for f in domain_footprints])
    for axs in domain_bounds:
        axs -= axs.min()
    domain = []
    for axis in out_frame.axes_order:
        axis_min, axis_max = domain_bounds[axis].min(), domain_bounds[axis].max()
        domain.append({'lower': axis_min, 'upper': axis_max,
                       'includes_lower': True, 'includes_upper': True})

    wnew.domain = domain
    return wnew


def spec_domain(slit, domain=None):
    """
    Returns a boolean mask of the pixels where the wcs is definied (not NaN).

    Build-7 workaround.

    Parameters
    ----------
    slit : datamodel

    Returns
    -------
    Boolean mask of valid wcs transform pixels are set to False.  Pixels
    in the input that produced NaNs in the output are set to True.
    """
    xstart, xend = slit.xstart - 1, slit.xstart - 1 + slit.xsize
    ystart, yend = slit.ystart - 1, slit.ystart - 1 + slit.ysize
    y, x = np.mgrid[ystart: yend, xstart: xend]
    ra, dec, lam = slit.meta.wcs(x, y)
    ra_masked = ma.array(ra, mask=np.isnan(ra))
    return ra_masked.mask


def grid_from_spec_domain(slit, domain=None):
    xstart, xend = slit.xstart - 1, slit.xstart - 1 + slit.xsize
    ystart, yend = slit.ystart - 1, slit.ystart - 1 + slit.ysize
    y, x = np.mgrid[ystart: yend, xstart: xend]
    mask = spec_domain(slit, domain=domain)
    return ma.array([x, y], mask=[mask, mask])


def spec_footprint(slit, domain=None):
    """
    Returns wcs footprint grid coordinates where NaNs are masked.

    Build-7 workaround.
    """
    x, y = grid_from_spec_domain(slit, domain=domain)
    ra, dec, lam = slit.meta.wcs(x, y)
    mask = y.mask
    return ma.array([ra, dec, lam], mask=[mask, mask, mask])


def wcs_from_spec_footprints(slitlist, refwcs=None, transform=None, domain=None):
    """
    Create a WCS from a list of datamodels with spectral WCS.

    Build-7 workaround.
    """
    if not isiterable(slitlist):
        raise ValueError("Expected 'slitlist' to be an iterable of datamodels.")
    if not all([isinstance(s.meta.wcs, WCS) for s in slitlist]):
        raise TypeError("All items in 'slitlist' must have instance of gwcs.WCS.")
    if refwcs is None:
        refwcs = slitlist[0].meta.wcs
        refslit = slitlist[0]
    else:
        if not isinstance(refwcs, WCS):
            raise TypeError("Expected refwcs to be an instance of gwcs.WCS.")
    fiducial = compute_spec_fiducial(slitlist, domain=domain)
     # trans = []
    # scales = [m for m in refwcs.forward_transform if isinstance(m, astmodels.Scale)]
    # if scales:
    #     trans.append(functools.reduce(lambda x, y: x & y, scales))
    # rotation = [m for m in refwcs.forward_transform if \
    #             isinstance(m, astmodels.AffineTransformation2D)]
    # if rotation:
    #     trans.append(rotation[0])
    # if trans:
    #     tr = functools.reduce(lambda x, y: x | y, trans)
    # else:
    #     tr = None
    proj = astmodels.Pix2Sky_TAN()
    transform = compute_spec_transform(fiducial, refslit)
    output_frame = getattr(refwcs, refwcs.output_frame)
    wnew = WCS(output_frame=output_frame, forward_transform=transform)
    footprints = [spec_footprint(s) for s in slitlist]
    # domain_bounds = np.hstack([wnew.backward_transform(f[0],f[1],f[2]) for f in footprints])
    # for axs in domain_bounds:
    #     axs -= axs.min()
    # domain = []
    # for axis in output_frame.axes_order:
    #     axis_min, axis_max = domain_bounds[axis].min(), domain_bounds[axis].max()
    #     domain.append({'lower': axis_min, 'upper': axis_max,
    #                    'includes_lower': True, 'includes_upper': True})

    # wnew.domain = domain
    return wnew, footprints


def compute_spec_fiducial(slitlist, domain=None):
    """
    For a celestial footprint this is the center.
    For a spectral footprint, it is the beginning of the range.

    This function assumes all WCSs have the same output coordinate frame.

    Build-7 workaround.
    """
    output_frame = getattr(slitlist[0].meta.wcs, 'output_frame')
    axes_types = getattr(slitlist[0].meta.wcs, output_frame).axes_type
    spatial_axes = np.array(axes_types) == 'SPATIAL'
    spectral_axes = np.array(axes_types) == 'SPECTRAL'
    footprints = ma.hstack([spec_footprint(s,
        domain=domain) for s in slitlist])
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


def compute_spec_transform(fiducial, slitmodel):
    """
    Compute a simple transform given a fidicial point in a spatial-spectral frame.
    """
    offset = 0.
    shift = astmodels.Shift(offset)
    rot = astmodels.Rotation2D(slitmodel.meta.wcsinfo.roll_ref)
    pixel_scale = 1.0
    scale = astmodels.Scale(pixel_scale)
    tan = astmodels.Pix2Sky_TAN()
    skyrot = astmodels.RotateNative2Celestial(fiducial[0], fiducial[1], 180.0)
    spatial = rot | tan | skyrot
    dispersion = 1.0
    spectral = astmodels.Scale(dispersion)
    mapping = astmodels.Mapping((0, 1, 0))
    transform = mapping | spatial & spectral
    return transform


def compute_fiducial(wcslist, domain=None):
    """
    For a celestial footprint this is the center.
    For a spectral footprint, it is the beginning of the range.

    This function assumes all WCSs have the same output coordinate frame.
    """
    output_frame = wcslist[0].output_frame
    axes_types = wcslist[0].output_frame.axes_type
    spatial_axes = np.array(axes_types) == 'SPATIAL'
    spectral_axes = np.array(axes_types) == 'SPECTRAL'
    footprints = np.hstack([w.footprint(domain=domain) for w in wcslist])
    spatial_footprint = footprints[spatial_axes]
    spectral_footprint = footprints[spectral_axes]

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
    return fiducial


def is_fits(input):
    """
    Returns
    --------
    isFits: tuple
        An ``(isfits, fitstype)`` tuple.  The values of ``isfits`` and
        ``fitstype`` are specified as:

         - ``isfits``: True|False
         - ``fitstype``: if True, one of 'waiver', 'mef', 'simple'; if False, None

    Notes
    -----
    Input images which do not have a valid FITS filename will automatically
    result in a return of (False, None).

    In the case that the input has a valid FITS filename but runs into some
    error upon opening, this routine will raise that exception for the calling
    routine/user to handle.
    """

    isfits = False
    fitstype = None
    names = ['fits', 'fit', 'FITS', 'FIT']
    #determine if input is a fits file based on extension
    # Only check type of FITS file if filename ends in valid FITS string
    f = None
    fileclose = False
    if isinstance(input, fits.HDUList):
        isfits = True
        f = input
    else:
        isfits = True in [input.endswith(l) for l in names]

    # if input is a fits file determine what kind of fits it is
    #waiver fits len(shape) == 3
    if isfits:
        if not f:
            try:
                f = fits.open(input, mode='readonly')
                fileclose = True
            except Exception:
                if f is not None:
                    f.close()
                raise
        data0 = f[0].data
        if data0 is not None:
            try:
                if isinstance(f[1], fits.TableHDU):
                    fitstype = 'waiver'
            except IndexError:
                fitstype = 'simple'

        else:
            fitstype = 'mef'
        if fileclose:
            f.close()

    return isfits, fitstype


def not_implemented_mode(input_model, ref):
    exp_type = input_model.meta.exposure.type
    message = "WCS for EXP_TYPE of {0} is not implemented.".format(exp_type)
    log.critical(message)
    #raise AttributeError(message)
    return None
