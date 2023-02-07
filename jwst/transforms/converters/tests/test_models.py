# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
import os
import warnings
from astropy.modeling.models import Shift, Rotation2D, Const1D
from asdf_astropy.converters.transform.tests.test_transform import (
    assert_model_roundtrip)

from jwst.transforms.models import (
    NirissSOSSModel, Rotation3DToGWA, Gwa2Slit, Logical, Slit,
    DirCos2Unitless, Unitless2DirCos, Snell, AngleFromGratingEquation,
    WavelengthFromGratingEquation)
import asdf
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


def test_niriss_soss(tmpdir):
    """Regression test for bugs discussed in issue #7401"""
    spectral_orders = [1, 2, 3]
    models = [
        Const1D(1.0) & Const1D(2.0) & Const1D(3.0),
        Const1D(4.0) & Const1D(5.0) & Const1D(6.0),
        Const1D(7.0) & Const1D(8.0) & Const1D(9.0),
    ]

    soss_model = NirissSOSSModel(spectral_orders, models)

    # Check that no warning is issued when serializing
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert_model_roundtrip(soss_model, tmpdir)

    # Check tag is the latest version
    path = str(tmpdir / "test.asdf")
    with asdf.AsdfFile({"model": soss_model}) as af:
        af.write_to(path)

    with asdf.open(path, _force_raw_types=True) as af:
        assert af.tree["model"]._tag == "tag:stsci.edu:jwst_pipeline/niriss_soss-1.0.0"


def test_niriss_soss_legacy():
    data_path = os.path.join(os.path.dirname(__file__), 'data')
    data = os.path.join(data_path, 'niriss_soss.asdf')

    # confirm that the file contains the legacy tag
    with asdf.open(data, _force_raw_types=True) as af:
        assert af.tree["model"]._tag == "tag:stsci.edu:jwst_pipeline/niriss-soss-0.7.0"

    # test that it opens with the legacy tag
    with asdf.open(data) as af:
        model = af['model']
        assert model.spectral_orders == [1, 2, 3]
        assert (model.models[1].parameters == (1.0, 2.0, 3.0)).all()
        assert (model.models[2].parameters == (4.0, 5.0, 6.0)).all()
        assert (model.models[3].parameters == (7.0, 8.0, 9.0)).all()
