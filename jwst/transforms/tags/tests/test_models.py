# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function
import numpy as np
from astropy.modeling.models import Shift, Rotation2D
from asdf.tests import helpers
from ...import jwextension
from ...models import *
from astropy.tests.helper import pytest


m1 = Shift(1) & Shift(2) | Rotation2D(3.1)
m2 = Shift(2) & Shift(2) | Rotation2D(23.1)


test_models = [DirCos2Unitless(), Unitless2DirCos(), NRSZCoord(),
               Rotation3DToGWA(angles=[12.1, 1.3, 0.5, 3.4], axes_order='xyzx'),
               AngleFromGratingEquation(20000, -1), WavelengthFromGratingEquation(25000, 2),
               Logical('GT', 5, 10), Logical('LT', np.ones((10,))* 5, np.arange(10))
               ]


@pytest.mark.parametrize(('model'), test_models)
def test_model(tmpdir, model):
    tree = {'model': model}
    helpers.assert_roundtrip_tree(tree, tmpdir, extensions=jwextension.JWSTExtension())


def test_gwa_to_slit(tmpdir):
    transforms = [m1, m2]
    s0 = Slit("s0", 1, 2, 3, 4, 5, 6, 7, 8)
    s1 = Slit("s1", 10, 20, 30, 40, 50, 60, 70, 80)
    slits = [s0, s1]
    m = Gwa2Slit(slits, transforms)
    tree = {'model': m}
    helpers.assert_roundtrip_tree(tree, tmpdir, extensions=jwextension.JWSTExtension())

    slits = [1, 2]
    m = Gwa2Slit(slits, transforms)
    tree = {'model': m}
    helpers.assert_roundtrip_tree(tree, tmpdir, extensions=jwextension.JWSTExtension())
