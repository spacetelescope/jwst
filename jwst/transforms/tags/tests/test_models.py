# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function

from astropy.modeling.models import Shift, Rotation2D
from asdf.tests import helpers
from ...import jwextension
from ...models import *
from astropy.tests.helper import pytest


m1 = Shift(1) & Shift(2) | Rotation2D(3.1)
m2 = Shift(2) & Shift(2) | Rotation2D(23.1)
gwa2slit_models = {1: m1, 2: m2}

test_models = [DirCos2Unitless(), Unitless2DirCos(), NRSZCoord(),
               Rotation3DToGWA(angles=[12.1, 1.3, 0.5, 3.4], axes_order='xyzx'),
               AngleFromGratingEquation(20000, -1), WavelengthFromGratingEquation(25000, 2),
               Gwa2Slit(gwa2slit_models)]


@pytest.mark.parametrize(('model'), test_models)
def test_model(tmpdir, model):
    tree = {'model': model}
    helpers.assert_roundtrip_tree(tree, tmpdir, extensions=jwextension.JWSTExtension())
