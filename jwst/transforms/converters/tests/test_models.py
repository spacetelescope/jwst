# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.modeling.models import Shift, Rotation2D
from asdf_astropy.converters.transform.tests.test_transform import (
    assert_model_roundtrip)

from jwst.transforms.models import (
    Rotation3DToGWA, Gwa2Slit, Logical, Slit,
    DirCos2Unitless, Unitless2DirCos, Snell, AngleFromGratingEquation,
    WavelengthFromGratingEquation)
import pytest


m1 = Shift(1) & Shift(2) | Rotation2D(3.1)
m2 = Shift(2) & Shift(2) | Rotation2D(23.1)


test_models = [DirCos2Unitless(), Unitless2DirCos(),
               Rotation3DToGWA(angles=[12.1, 1.3, 0.5, 3.4], axes_order='xyzx'),
               AngleFromGratingEquation(20000, -1), WavelengthFromGratingEquation(25000, 2),
               Logical('GT', 5, 10), Logical('LT', np.ones((10,)) * 5, np.arange(10)),
               Snell(angle=-16.5, kcoef=[0.583, 0.462, 3.891], lcoef=[0.002526, 0.01, 1200.556],
                     tcoef=[-2.66e-05, 0.0, 0.0, 0.0, 0.0, 0.0], tref=35, pref=0,
                     temperature=35, pressure=0),
               ]


@pytest.mark.parametrize(('model'), test_models)
def test_model(tmpdir, model, version=None):
    assert_model_roundtrip(model, tmpdir)


def test_gwa_to_slit(tmpdir):
    transforms = [m1, m2]
    s0 = Slit("s0", 1, 2, 3, 4, 5, 6, 7, 8)
    s1 = Slit("s1", 10, 20, 30, 40, 50, 60, 70, 80)
    slits = [s0, s1]
    m = Gwa2Slit(slits, transforms)
    assert_model_roundtrip(m, tmpdir)

    slits = [1, 2]
    m = Gwa2Slit(slits, transforms)
    assert_model_roundtrip(m, tmpdir)
