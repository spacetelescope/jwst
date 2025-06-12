"""Unit tests for AMI mask_definition_ami module."""

import pytest
import numpy as np

from jwst.ami.mask_definition_ami import NRMDefinition


@pytest.fixture
def ctrs_asdesigned(nrm_model):
    return np.array(
        [
            [nrm_model.x_a1, nrm_model.y_a1],
            [nrm_model.x_a2, nrm_model.y_a2],
            [nrm_model.x_a3, nrm_model.y_a3],
            [nrm_model.x_a4, nrm_model.y_a4],
            [nrm_model.x_a5, nrm_model.y_a5],
            [nrm_model.x_a6, nrm_model.y_a6],
            [nrm_model.x_a7, nrm_model.y_a7],
        ]
    )


def test_NRMDefinition(nrm_model, ctrs_asdesigned):
    nrm = NRMDefinition(nrm_model)

    # test copy over of attributes
    assert nrm.maskname == "jwst_ami"
    assert nrm.hdia == nrm_model.flat_to_flat
    assert nrm.active_D == nrm_model.diameter
    assert nrm.OD == nrm_model.pupil_circumscribed

    # test hole centers
    # this tests NRMDefinition.read_nrm_model implicitly because
    # it calculates self.ctrs during initialization
    assert isinstance(nrm.ctrs, np.ndarray)
    assert nrm.ctrs.shape == (7, 2)
    assert np.allclose(nrm.ctrs, -np.fliplr(ctrs_asdesigned))


def test_NRMDefinition_chooseholes(nrm_model, ctrs_asdesigned):
    # test the options for chooseholes
    chooseholes = ["B2", "B4", "B5", "B6"]  # indices [3, 0, 2, 5]
    nrm = NRMDefinition(nrm_model, chooseholes=chooseholes)
    assert nrm.ctrs.shape == (4, 2)
    assert np.allclose(nrm.ctrs, -np.fliplr(ctrs_asdesigned)[[3, 0, 2, 5]])

    with pytest.raises(ValueError):
        NRMDefinition(nrm_model, maskname="jwst_unsupported")
