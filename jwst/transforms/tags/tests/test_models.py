# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function

from numpy.testing import assert_allclose
from asdf import AsdfFile
from asdf.tests import helpers
from ...import jwextension
from ...models import *
from astropy.tests.helper import pytest


test_models = [DirCos2Unitless(), Unitless2DirCos(), NRSZCoord(),
          Rotation3DToGWA(angles=[12.1, 1.3, 0.5, 3.4], axes_order='xyzx'),
          AngleFromGratingEquation(20000, -1), WavelengthFromGratingEquation(25000, 2)
          ]


@pytest.mark.parametrize(('model'), test_models)
def test_model(tmpdir, model):
    tree = {'model': model}
    helpers.assert_roundtrip_tree(tree, tmpdir, extensions=jwextension.JWSTExtension())


