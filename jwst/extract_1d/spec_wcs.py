from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from astropy.modeling.models import Mapping, Const1D, Tabular1D
from astropy import units as u
from astropy import coordinates as coord

from gwcs.wcs import WCS
from gwcs import coordinate_frames as cf


def create_spectral_wcs(ra, dec, wavelength):
    """Assign a WCS for sky coordinates and a table of wavelengths

    Parameters
    ----------
    ra: float
        The right ascension (in degrees) at the nominal location of the
        entrance aperture (slit).

    dec: float
        The declination (in degrees) at the nominal location of the
        entrance aperture.

    wavelength: ndarray
        The wavelength in microns at each pixel of the extracted spectrum.

    Returns:
    --------
    wcs: a gwcs.wcs.WCS object
        This takes a float or sequence of float and returns a tuple of
        the right ascension, declination, and wavelength (or sequence of
        wavelengths) at the pixel(s) specified by the input argument.
    """

    # Only the first coordinate is used.
    input_frame = cf.Frame2D(axes_order=(0, 1), unit=(u.pix, u.pix),
                             name="pixel_frame")

    sky = cf.CelestialFrame(name='sky', axes_order=(0, 1),
                            reference_frame=coord.ICRS())
    spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.micron,),
                            axes_names=('wavelength',))

    pixel = np.arange(len(wavelength), dtype=np.float)
    tab = Mapping((0, 0, 0)) | \
          Const1D(ra) & Const1D(dec) & Tabular1D(pixel, wavelength)

    world = cf.CompositeFrame([sky, spec], name='world')

    pipeline = [(input_frame, tab),
                (world, None)]

    return WCS(pipeline)
