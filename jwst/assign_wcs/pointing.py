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
        transform = fitswcs_transform_from_model(wcsinfo) #inp)
        output_frame = frame_from_model(wcsinfo) #inp)
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
