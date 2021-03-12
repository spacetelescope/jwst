import logging
import numpy as np

from astropy.modeling.models import Mapping, Const1D, Tabular1D
from astropy import units as u
from astropy import coordinates as coord

from gwcs.wcs import WCS
from gwcs import coordinate_frames as cf

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


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

    input_frame = cf.CoordinateFrame(naxes=1, axes_type=("SPATIAL",),
                                     axes_order=(0,), unit=(u.pix,),
                                     name="pixel_frame")

    sky = cf.CelestialFrame(name='sky', axes_order=(0, 1),
                            reference_frame=coord.ICRS())
    spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.micron,),
                            axes_names=('wavelength',))
    world = cf.CompositeFrame([sky, spec], name='world')

    pixel = np.arange(len(wavelength), dtype=float)
    tab = Mapping((0, 0, 0)) | \
        Const1D(ra) & Const1D(dec) & Tabular1D(points=pixel,
                                               lookup_table=wavelength,
                                               bounds_error=False,
                                               fill_value=None)
    tab.name = "pixel_to_world"

    if all(np.diff(wavelength) > 0):
        tab.inverse = Mapping((2,)) | Tabular1D(points=wavelength,
                                                lookup_table=pixel,
                                                bounds_error=False,
                                                )
    elif all(np.diff(wavelength) < 0):
        tab.inverse = Mapping((2,)) | Tabular1D(points=wavelength[::-1],
                                                lookup_table=pixel[::-1],
                                                bounds_error=False,
                                                )
    else:
        log.warning("Wavelengths are not strictly monotonic, inverse transform is not set")

    pipeline = [(input_frame, tab),
                (world, None)]

    return WCS(pipeline)
