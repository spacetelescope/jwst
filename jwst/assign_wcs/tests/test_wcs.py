from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from numpy.testing import utils
from .. import pointing
from ...transforms import models


def test_roll_angle():
    """
    A sanity test - when V2_REF = 0 and V3_REF = 0,
    for V2, V3 close to he origin
    ROLL_REF should be approximately PA_V3 .

    (Test taken from SIAF report.)
    """
    ra_ref = 165 # in deg
    dec_ref = 54 # in deg
    v2_ref = 0
    v3_ref = 0
    r0 = 37 # in deg

    v2 = .01 # in arcsec
    v3 = .01 # in arcsec
    roll_angle = pointing.compute_roll_ref(v2_ref, v3_ref, r0, ra_ref, dec_ref, v2, v3)
    utils.assert_allclose(roll_angle, r0, atol=1e-3)


def test_v23_to_sky():
    """
    Test taken from INS report.
    """
    ra_ref = 165 # in deg
    dec_ref = 54 # in deg
    v2_ref = -503.654472 / 3600 # in deg
    v3_ref = -318.742464 / 3600 # in deg
    r0 = 37 # in deg

    v2 = 210 # in deg
    v3 = -75 # in deg
    expected_ra_dec = (107.12810484789563, -35.97940247128502) # in deg
    angles = [-v2_ref, v3_ref, -r0, -dec_ref, ra_ref]
    axes = "zyxyz"
    v2s = models.V23ToSky(angles, axes_order=axes)
    radec = v2s(v2, v3)
    utils.assert_allclose(radec, expected_ra_dec, atol=1e-10)
