from os import path

import numpy as np
import asdf
from astropy.coordinates import SkyCoord
from astropy import units as u

from jwst.tweakreg.tests import data
from jwst.tweakreg.utils import adjust_wcs


data_path = path.split(path.abspath(data.__file__))[0]


def test_adjust_wcs():
    wcs_file = path.join(data_path, 'nrcb1-wcs.asdf')
    w0 = asdf.open(wcs_file)['wcs']

    crval20, crval10 = w0.pipeline[-2].transform.angles_3.value.tolist()[-2:]
    crval10 = -crval10

    crpix10, crpix20 = w0.numerical_inverse(crval10, crval20)

    wa = adjust_wcs(
        w0,
        delta_ra=0.0135,
        delta_dec=0.0208,
        delta_roll=25.7,
        scale_factor=1.003
    )

    crval1a, crval2a = wa(crpix10, crpix20)

    assert np.allclose(
        [crval1a - crval10, crval2a - crval20],
        [0.0135, 0.0208],
        rtol=0,
        atol=1e-13
    )

    offset_ra10, offset_dec10 = w0(crpix10, crpix20 + 10)
    offset_ra1a, offset_dec1a = wa(crpix10, crpix20 + 10)

    ca0 = SkyCoord(crval1a * u.deg, crval2a * u.deg, frame='icrs')
    ca1 = SkyCoord(offset_ra1a * u.deg, offset_dec1a * u.deg, frame='icrs')
    c0 = SkyCoord(crval10 * u.deg, crval20 * u.deg, frame='icrs')
    c01 = SkyCoord(offset_ra10 * u.deg, offset_dec10 * u.deg, frame='icrs')

    sep0 = c0.separation(c01).degree
    sepa = ca0.separation(ca1).degree

    # test scale:
    assert np.allclose(sepa / sep0, 1.003, rtol=0, atol=1e-8)

    # test roll:
    assert np.allclose(
        ca1.position_angle(ca0).degree - c01.position_angle(c0).degree,
        25.7,
        rtol=1.0e-5,
        atol=0.0
    )
