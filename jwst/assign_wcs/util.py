"""
Utility function for WCS

"""

import logging
import functools
import numpy as np

from astropy.utils.misc import isiterable
from astropy.io import fits
from astropy.modeling import projections
from astropy.modeling import models as astmodels

from gwcs import WCS
from gwcs.wcstools import wcs_from_fiducial

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


def wcs_from_footprints(dmodels, refmodel=None, transform=None, domain=None):
    """
    Create a WCS from a list of input models.

    A fiducial point in the output coordinate frame is created from  the
    footprints of all WCS objects. For a spatial frame this is the center
    of the union of the footprints. For a spectral frame the fiducial is in
    the beginning of the footprint range.
    If ``refmodel`` is None, the first WCS object in the list is considered
    a reference. The output coordinate frame and projection (for celestial frames)
    is taken from ``refmodel``.
    If ``transform`` is not suplied, a compound transform is made using
    CDELTs and rotation by ROLL_REF.
    If ``domain`` is not supplied, the domain of the new WCS is computed
    from the domains of all input WCSs.

    Parameters
    ----------
    dmodels : list of `~jwst.datamodels.DataModel`
        A list of data models.
    refmodel : `~jwst.datamodels.DataModel`, optional
        This model's WCS is used as a reference.
        WCS. The output coordinate frame, the projection and a
        scaling and rotation transform is created from it. If not supplied
        the first model in the list is used as ``refmodel``.
    transform : `~astropy.modeling.core.Model`, optional
        A transform, passed to :class_method:`~gwcs.WCS.wcs_from_fiducial`
        If not supplied Scaling | Rotation is computed from ``refmodel``.
    domain : list of dicts, optional
        Domain of the new WCS.
        If not supplied it is computed from the domain of all inputs.
    """
    wcslist = [im.meta.wcs for im in dmodels]
    if not isiterable(wcslist):
        raise ValueError("Expected 'wcslist' to be an iterable of WCS objects.")
    if not all([isinstance(w, WCS) for w in wcslist]):
        raise TypeError("All items in wcslist are to be instances of gwcs.WCS.")
    if refmodel is None:
        refmodel = dmodels[0]
    else:
        if not isinstance(refmodel, DataModel):
            raise TypeError("Expected refmodel to be an instance of DataModel.")

    fiducial = compute_fiducial(wcslist, domain)

    prj = astmodels.Pix2Sky_TAN()
    trans = []

    scales = astmodels.Scale(refmodel.meta.wcsinfo.cdelt1) & \
        astmodels.Scale(refmodel.meta.wcsinfo.cdelt2)
    trans.append(scales)

    if trans:
        tr = functools.reduce(lambda x, y: x | y, trans)
    else:
        tr = None
    rotation = astmodels.Rotation2D._compute_matrix(np.deg2rad(refmodel.meta.wcsinfo.roll_ref))
    trans.append(rotation)
    out_frame = refmodel.meta.wcs.output_frame
    wnew = wcs_from_fiducial(fiducial, coordinate_frame=out_frame,
                             projection=prj, transform=tr)

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
