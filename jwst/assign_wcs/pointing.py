from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import six
import numpy as np
from astropy import units as u
from astropy import coordinates as coords
from astropy.modeling import models as astmodels
from ..datamodels import fits_support, DataModel
from gwcs import utils as gwutils
from gwcs import coordinate_frames as cf
from gwcs import wcs
from ..transforms.models import V23ToSky


def v23tosky(input_model):
    v2_ref = input_model.meta.wcsinfo.v2_ref / 3600
    v3_ref = input_model.meta.wcsinfo.v3_ref / 3600
    roll_ref = input_model.meta.wcsinfo.roll_ref
    ra_ref = input_model.meta.wcsinfo.ra_ref
    dec_ref = input_model.meta.wcsinfo.dec_ref
    
    angles = [-v2_ref, v3_ref, -roll_ref, -dec_ref, ra_ref]
    axes = "zyxyz"
    sky_rotation = V23ToSky(angles, axes_order=axes, name="v23tosky")
    return sky_rotation


def compute_roll_ref(v2_ref, v3_ref, roll_ref, ra_ref, dec_ref, new_v2_ref, new_v3_ref):
    """
    Computes the position of V3 (measured N to E) at the center af an aperture.

    Parameters
    ----------
    v2_ref, v3_ref : float
        Reference point in the V2, V3 frame [in arcsec] (FITS keywords V2_REF and V3_REF)
    roll_ref : float
        Position angle of V3 at V2_REF, V3_REF, [in deg]
        When ROLL_REF == PA_V2, V2_REF, V3_REF = (0, 0)
    ra_ref, dec_ref : float
        RA and DEC corresponding to V2_REF and V3_REF, [in deg]
    new_v2_ref, new_v3_ref : float
        The new position in V2, V3 where the position of V3 is computed, [in arcsec]
        The center of the aperture in V2,V3

    Returns
    -------
    new_roll : float
        The value of ROLL_REF (in deg)

    """
    v2 = np.deg2rad(new_v2_ref / 3600)
    v3 = np.deg2rad(new_v3_ref / 3600)

    if np.isclose(v2_ref, 0, atol=1e-13) and np.isclose(v3_ref, 0, atol=1e-13):
        angles = [-roll_ref, -dec_ref, - ra_ref]
        axes = "xyz"
    else:
        angles = [-v2_ref, v3_ref, -roll_ref, -dec_ref, ra_ref]
        axes = "zyxyz"
    M = V23ToSky._compute_matrix(np.deg2rad(angles), axes)
    return _roll_angle_from_matrix(M, v2, v3)


def _roll_angle_from_matrix(matrix, v2, v3):
    X = -(matrix[2, 0] * np.cos(v2) + matrix[2, 1] * np.sin(v2)) * \
        np.sin(v3) + matrix[2, 2] * np.cos(v3)
    Y = (matrix[0, 0] * matrix[1, 2] - matrix[1, 0] * matrix[0, 2]) * np.cos(v2) + \
        (matrix[0, 1] * matrix[1, 2] - matrix[1, 1] * matrix[0, 2]) * np.sin(v2)
    new_roll = np.rad2deg(np.arctan2(Y, X))
    if new_roll < 0:
        new_roll += 360
    return new_roll


def create_fitswcs_transform(input_model):
    ff = fits_support.to_fits(input_model._instance, input_model._schema)
    hdu = fits_support.get_hdu(ff._hdulist, "PRIMARY")
    header = hdu.header
    transform = gwutils.make_fitswcs_transform(header)
    return transform


def wcsinfo_from_model(input_model):
    defaults = {'crpix': 0, 'crval': 0, 'cdelt': 1}
    wcsinfo = {}
    wcsaxes = input_model.meta.wcsinfo.wcsaxes
    wcsinfo['WCSAXES'] = wcsaxes
    for key in ['CRPIX', 'CRVAL', 'CDELT']:
        val = []
        for ax in range(1, wcsaxes + 1):
            k = (key + "{0}".format(ax)).lower()
            try:
                v = getattr(input_model.meta.wcsinfo, k)
                if v is None:
                    v = defaults[key]
                val.append(v)
            except KeyError:
                val.append(defaults[key])
        wcsinfo[key.upper()] = np.array(val)

    ctypes = []
    cunits = []
    for ax in range(1, wcsaxes + 1):
        try:
            ctypes.append(getattr(input_model.meta.wcsinfo, ("CTYPE{0}".format(ax)).lower()))
        except KeyError:
            ctypes.append("")
        try:
            v = getattr(input_model.meta.wcsinfo, ("CUNIT{0}".format(ax)).lower())
            cunits.append(u.Unit(v))
        except KeyError:
            cunits.append(u.Unit(""))

    wcsinfo["CTYPE"] = np.array(ctypes)
    wcsinfo["CUNIT"] = np.array(cunits)

    pc = np.zeros((wcsaxes, 2))

    for i in range(1, wcsaxes + 1):
        for j in range(1, 3):
            pc[i - 1, j - 1] = getattr(input_model.meta.wcsinfo, 'pc{0}_{1}'.format(i, j))
    wcsinfo['PC'] = pc
    wcsinfo['RADESYS'] = input_model.meta.coordinates.reference_frame
    wcsinfo['has_cd'] = False
    return wcsinfo


def fitswcs_transform_from_model(wcsinfo):
    """
    Create a fits wcs from a datamodel.meta.wcsinfo.
    """
    spatial_axes, spectral_axes = gwutils.get_axes(wcsinfo)

    transform = gwutils.make_fitswcs_transform(wcsinfo)
    if wcsinfo['WCSAXES'] == 3:
        spectral_transform = astmodels.Shift(-wcsinfo['CRPIX'][spectral_axes]) | \
                           astmodels.Scale(wcsinfo['CDELT'][spectral_axes]) | \
                           astmodels.Shift(wcsinfo['CRVAL'][spectral_axes])
        transform = transform & spectral_transform

    return transform


def frame_from_fits(ff):
    raise NotImplementedError()


def frame_from_model(wcsinfo):
    wcsaxes = wcsinfo['WCSAXES']
    spatial_axes, spectral_axes = gwutils.get_axes(wcsinfo)
    cunit = wcsinfo['CUNIT']
    ref_frame = coords.ICRS()
    sky = cf.CelestialFrame(name='sky', axes_order=tuple(spatial_axes), reference_frame=ref_frame,
                            unit=cunit[spatial_axes], axes_names=('RA', 'DEC'))
    if wcsaxes == 2:
        return sky
    elif wcsaxes == 3:
        spec = cf.SpectralFrame(name='spectral', axes_order=tuple(spectral_axes), unit=cunit[spectral_axes[0]],
                                axes_names=('wavelength',))
        world = cf.CompositeFrame([sky, spec], name='world')
        return world
    else:
        raise ValueError("WCSAXES can be 2 or 3, git {0}".format(wcsaxes))


def create_fitswcs(inp):
    if isinstance(inp, DataModel):
        wcsinfo = wcsinfo_from_model(inp)
        transform = fitswcs_transform_from_model(wcsinfo)
        output_frame = frame_from_model(wcsinfo)
    elif isinstance(inp, six.string_types):
        transform = create_fitswcs_transform(inp)
        output_frame = frame_from_fits(inp)
    else:
        raise TypeError("Input is expected to be a DataModel instance or a FITS file.")

    input_frame = cf.Frame2D(name='detector', axes_order=(0, 1))
    pipeline = [(input_frame, transform),
               (output_frame, None)]

    wcsobj = wcs.WCS(pipeline)
    return wcsobj
