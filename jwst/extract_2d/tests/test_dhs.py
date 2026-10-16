"""Tests for extract_2d on NIRCam DHS (multi-stripe) mode."""

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.helpers import (
    NRCA1_DHS_STRIPE_IDS,
    make_mock_dhs_nrca1_rate,
    make_mock_dhs_nrcalong_rate,
)
from jwst.extract_2d.extract_2d_step import Extract2dStep
from jwst.extract_2d.grisms import extract_tso_object


def get_reference_files(datamodel):
    """Return the extract_2d reference files for the input exposure."""
    refs = {}
    step = Extract2dStep()
    for reftype in Extract2dStep.reference_file_types:
        refs[reftype] = step.get_reference_file(datamodel, reftype)
    return refs


@pytest.fixture
def mock_dhs_nrca1_rate():
    """Mock DHS NRCA1 rate CubeModel."""
    return make_mock_dhs_nrca1_rate()


@pytest.fixture
def mock_dhs_nrcalong_rate():
    """Mock DHS NRCALONG rate CubeModel."""
    return make_mock_dhs_nrcalong_rate()


@pytest.fixture
def dhs_nrca1_wcs_model(mock_dhs_nrca1_rate):
    """NRCA1 DHS CubeModel with WCS assigned."""
    return AssignWcsStep.call(mock_dhs_nrca1_rate)


@pytest.fixture
def dhs_nrcalong_wcs_model(mock_dhs_nrcalong_rate):
    """NRCALONG DHS CubeModel with WCS assigned."""
    return AssignWcsStep.call(mock_dhs_nrcalong_rate)


def test_extract_tso_dhs_nrca1(dhs_nrca1_wcs_model):
    """Extract NRCA1 DHS data."""
    refs = get_reference_files(dhs_nrca1_wcs_model)
    result = extract_tso_object(dhs_nrca1_wcs_model, reference_files=refs)

    assert isinstance(result, datamodels.MultiSlitModel)
    assert len(result.slits) == len(NRCA1_DHS_STRIPE_IDS)
    assert {slit.name for slit in result.slits} == {str(sid) for sid in NRCA1_DHS_STRIPE_IDS}
    full_width = dhs_nrca1_wcs_model.meta.subarray.xsize
    for slit in result.slits:
        assert slit.xsize == full_width
        assert slit.ysize > 0
        # Verify the correct detector region was extracted: each stripe
        # row range was filled with the stripe ID value in the fixture.
        assert np.all(slit.data == int(slit.name))


def test_extract_tso_dhs_nrcalong(dhs_nrcalong_wcs_model):
    """
    Extract NRCALONG DHS data.

    Output is a MultiSlitModel with four slits (one per readout)
    representing the same physical location on the detector.
    """
    refs = get_reference_files(dhs_nrcalong_wcs_model)
    result = extract_tso_object(dhs_nrcalong_wcs_model, reference_files=refs)

    assert isinstance(result, datamodels.MultiSlitModel)
    assert len(result.slits) == 4
    full_width = dhs_nrcalong_wcs_model.meta.subarray.xsize
    ref_data = result.slits[0].data
    for slit in result.slits:
        assert slit.xsize == full_width
        assert slit.ysize > 0
        # All NRCALONG stripes read from the same physical detector band
        # ensure they are all identical
        assert np.all(slit.data == ref_data)


def test_extract_2d_step_dhs_nrca1(dhs_nrca1_wcs_model):
    """Extract2dStep completes on NRCA1 DHS data and does not modify the input."""
    result = Extract2dStep.call(dhs_nrca1_wcs_model)

    assert result.meta.cal_step.extract_2d == "COMPLETE"
    assert isinstance(result, datamodels.MultiSlitModel)
    assert len(result.slits) == len(NRCA1_DHS_STRIPE_IDS)
    assert result is not dhs_nrca1_wcs_model
    assert dhs_nrca1_wcs_model.meta.cal_step.extract_2d is None


def test_extract_2d_step_dhs_nrcalong(dhs_nrcalong_wcs_model):
    """Extract2dStep completes on NRCALONG DHS data and does not modify the input."""
    result = Extract2dStep.call(dhs_nrcalong_wcs_model)

    assert result.meta.cal_step.extract_2d == "COMPLETE"
    assert isinstance(result, datamodels.MultiSlitModel)
    assert len(result.slits) == len(NRCA1_DHS_STRIPE_IDS)
    assert result is not dhs_nrcalong_wcs_model
    assert dhs_nrcalong_wcs_model.meta.cal_step.extract_2d is None
