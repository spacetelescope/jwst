import numpy as np
from astropy import units as u
from astropy import coordinates as coords
from astropy.modeling import models as astmodels
from ..datamodels import DataModel
from gwcs import utils as gwutils
from gwcs import coordinate_frames as cf
from gwcs import wcs
from ..transforms.models import V23ToSky


__all__ = ["compute_roll_ref", "frame_from_model", "fitswcs_transform_from_model"]


def v23tosky(input_model):
    v2_ref = input_model.meta.wcsinfo.v2_ref / 3600
    v3_ref = input_model.meta.wcsinfo.v3_ref / 3600
    roll_ref = input_model.meta.wcsinfo.roll_ref
    ra_ref = input_model.meta.wcsinfo.ra_ref
    dec_ref = input_model.meta.wcsinfo.dec_ref

    angles = [-v2_ref, v3_ref, -roll_ref, -dec_ref, ra_ref]
    axes = "zyxyz"
    sky_rotation = V23ToSky(angles, axes_order=axes, name="v23tosky")
    # The sky rotation expects values in deg.
    # This should be removed when models work with quantities.
    return astmodels.Scale(1/3600) & astmodels.Scale(1/3600) | sky_rotation


def compute_roll_ref(v2_ref, v3_ref, roll_ref, ra_ref, dec_ref, new_v2_ref, new_v3_ref):
    """
    Computes the position of V3 (measured N to E) at the center af an aperture.

    Parameters
    ----------
    v2_ref, v3_ref : float
        Reference point in the V2, V3 frame [in arcsec] (FITS keywords V2_REF and V3_REF)
    roll_ref : float
        Position angle of V3 at V2_REF, V3_REF, [in deg]
        When ROLL_REF == PA_V3, then (V2_REF, V3_REF) = (0, 0)
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

    v2_ref = v2_ref / 3600
    v3_ref = v3_ref / 3600

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


#def create_fitswcs_transform(input_model):
    #"""

    #"""
    #ff = fits_support.to_fits(input_model._instance, input_model._schema)
    #hdu = fits_support.get_hdu(ff._hdulist, "PRIMARY")
    #header = hdu.header
    #transform = gwutils.make_fitswcs_transform(header)
    #return transform


def wcsinfo_from_model(input_model):
    """
    Create a dict {wcs_keyword: array_of_values} pairs from a data model.

    Parameters
    ----------
    input_model : `~jwst.datamodels.model_base.DataModel`
        The input data model

    """
    defaults = {'CRPIX': 0, 'CRVAL': 0, 'CDELT': 1., 'CTYPE': "", 'CUNIT': u.Unit("")}
    wcsinfo = {}
    wcsaxes = input_model.meta.wcsinfo.wcsaxes
    wcsinfo['WCSAXES'] = wcsaxes
    for key in ['CRPIX', 'CRVAL', 'CDELT', 'CTYPE', 'CUNIT']:
        val = []
        for ax in range(1, wcsaxes + 1):
            k = (key + "{0}".format(ax)).lower()
            v = getattr(input_model.meta.wcsinfo, k, defaults[key])
            val.append(v)
        wcsinfo[key] = np.array(val)

    pc = np.zeros((wcsaxes, wcsaxes))
    for i in range(1, wcsaxes + 1):
        for j in range(1, wcsaxes + 1):
            pc[i - 1, j - 1] = getattr(input_model.meta.wcsinfo, 'pc{0}_{1}'.format(i, j), 1)
    wcsinfo['PC'] = pc
    wcsinfo['RADESYS'] = input_model.meta.coordinates.reference_frame
    wcsinfo['has_cd'] = False
    return wcsinfo


def fitswcs_transform_from_model(wcsinfo):
    """
    Create a WCS object using from datamodel.meta.wcsinfo.
    Transforms assume 0-based coordinates.

    Parameters
    ----------
    wcsinfo : dict-like
        ``~jwst.meta.wcsinfo`` structure.

    Return
    ------
    transform : `~astropy.modeling.core.Model`
        WCS forward transform - from pixel to world coordinates.

    """
    spatial_axes, spectral_axes, unknown = gwutils.get_axes(wcsinfo)
    #sp_axis = spectral_axes[0]

    transform = gwutils.make_fitswcs_transform(wcsinfo)
    #if wcsinfo['WCSAXES'] == 3:
    if spectral_axes:
        sp_axis = spectral_axes[0]
        # Subtract one from CRPIX which is 1-based.
        spectral_transform = astmodels.Shift(-(wcsinfo['CRPIX'][sp_axis] - 1)) | \
                           astmodels.Scale(wcsinfo['CDELT'][sp_axis]) | \
                           astmodels.Shift(wcsinfo['CRVAL'][sp_axis])
        transform = transform & spectral_transform

    return transform


def frame_from_fits(ff):
    raise NotImplementedError


def frame_from_model(wcsinfo):
    """
    Initialize a coordinate frame based on values in model.meta.wcsinfo.

    Parameters
    ----------
    wcsinfo : `~jwst.datamodels.model_base.DataModel` or dict
        Either one of the JWST data moels or a dict with model.meta.wcsinfo.

    Returns
    -------
    frame : `~coordinate_frames.CoordinateFrame`

    """
    if isinstance(wcsinfo, DataModel):
        wcsinfo = wcsinfo_from_model(wcsinfo)

    wcsaxes = wcsinfo['WCSAXES']
    celestial_axes, spectral_axes, other = gwutils.get_axes(wcsinfo)
    cunit = wcsinfo['CUNIT']
    frames = []
    if celestial_axes:
        ref_frame = coords.ICRS()
        celestial = cf.CelestialFrame(name='sky', axes_order=tuple(celestial_axes),
                                    reference_frame=ref_frame, unit=cunit[celestial_axes],
                                    axes_names=('RA', 'DEC'))
        frames.append(celestial)
    if spectral_axes:
        spec = cf.SpectralFrame(name='spectral', axes_order=tuple(spectral_axes),
                                unit=cunit[spectral_axes],
                                axes_names=('wavelength',))
        frames.append(spec)
    if other:
        # Make sure these are strings and not np.str_ objects.
        axes_names = tuple([str(name) for name in wcsinfo['CTYPE'][other]])
        name = "_".join(wcsinfo['CTYPE'][other])
        spatial = cf.Frame2D(name=name, axes_order=tuple(other), unit=cunit[other],
                             axes_names=axes_names)
        frames.append(spatial)
    if wcsaxes == 2:
        return frames[0]
    elif wcsaxes == 3:
        world = cf.CompositeFrame(frames, name='world')
        return world
    else:
        raise ValueError("WCSAXES can be 2 or 3, got {0}".format(wcsaxes))


def create_fitswcs(inp, input_frame=None):
    if isinstance(inp, DataModel):
        wcsinfo = wcsinfo_from_model(inp)
        transform = fitswcs_transform_from_model(wcsinfo)
        output_frame = frame_from_model(wcsinfo)
    #elif isinstance(inp, str):
        #transform = create_fitswcs_transform(inp)
        #output_frame = frame_from_fits(inp)
    else:
        raise TypeError("Input is expected to be a DataModel instance or a FITS file.")

    if input_frame is None:
        input_frame = "detector"
    pipeline = [(input_frame, transform),
               (output_frame, None)]

    wcsobj = wcs.WCS(pipeline)
    return wcsobj
