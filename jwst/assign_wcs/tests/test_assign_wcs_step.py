"""Test the AssignWCSStep."""

import asdf
import numpy as np
import pytest
from gwcs import coordinate_frames as cf
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_miri import create_hdul as create_miri
from jwst.assign_wcs.tests.test_nircam import create_hdul as create_nircam
from jwst.assign_wcs.tests.test_niriss import create_hdul as create_niriss
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_ifu_file
from jwst.assign_wcs.util import NoDataOnDetectorError


def test_assign_wcs_step_miri_ifu():
    hdul = create_miri(detector="MIRIFULONG", channel="34", band="MEDIUM")
    hdul[1].data = np.zeros((3, 40, 50))
    model = datamodels.CubeModel(hdul)
    result = AssignWcsStep.call(model)
    assert result is not model
    assert result.meta.cal_step.assign_wcs == "COMPLETE"
    assert model.meta.cal_step.assign_wcs is None


def test_assign_wcs_step_nis_wfss():
    hdul = create_niriss(filtername="GR150R", pupil="F200W", exptype="NIS_WFSS")
    model = datamodels.ImageModel(hdul)
    result = AssignWcsStep.call(model)
    assert result is not model
    assert result.meta.cal_step.assign_wcs == "COMPLETE"
    assert model.meta.cal_step.assign_wcs is None


def test_assign_wcs_step_nrc_wfss():
    hdul = create_nircam(exptype="NRC_WFSS", filtername="F444W", pupil="GRISMR")
    model = datamodels.ImageModel(hdul)
    model.data = np.zeros((10, 10))
    result = AssignWcsStep.call(model)
    assert result is not model
    assert result.meta.cal_step.assign_wcs == "COMPLETE"
    assert model.meta.cal_step.assign_wcs is None


def test_unsupported_input(caplog):
    model = datamodels.SlitModel()
    result = AssignWcsStep.call(model)
    assert result is not model
    assert "type is not supported" in caplog.text
    assert result.meta.cal_step.assign_wcs == "SKIPPED"
    assert model.meta.cal_step.assign_wcs is None


def test_assign_wcs_step_nrs_ifu_coord_wcs():
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    model = datamodels.IFUImageModel(hdul)
    hdul.close()

    result = AssignWcsStep.call(model, nrs_ifu_slice_wcs=False)
    assert result is not model
    assert result.meta.cal_step.assign_wcs == "COMPLETE"
    assert model.meta.cal_step.assign_wcs is None

    # The first frame is an identity transform from coordinates to detector
    assert result.meta.wcs.available_frames[0] == "coordinates"
    assert isinstance(result.meta.wcs.pipeline[0].frame, cf.Frame2D)
    assert result.meta.wcs.transform("coordinates", "detector", 1, 1) == (1, 1)


@pytest.mark.parametrize("valid_data", [True, False])
def test_assign_wcs_step_nrs_ifu_m_grating(tmp_path, valid_data):
    hdul = create_nirspec_ifu_file(grating="G235M", filter="F170LP")
    model = datamodels.IFUImageModel(hdul)
    model.meta.instrument.detector = "NRS2"
    hdul.close()

    # Mock a wavelength range file to either cover the NRS2 wavelengths or not
    if valid_data:
        wrange = [1.66e-06, 5.5e-06]
    else:
        wrange = [1.66e-06, 3.17e-06]
    wavelengthrange = {
        "order": [-1],
        "wavelengthrange": [wrange],
        "waverange_selector": ["F170LP_G235M"],
    }
    mock_reffile = str(tmp_path / "wrange.asdf")
    af = asdf.AsdfFile(tree=wavelengthrange)
    af.write_to(mock_reffile)
    af.close()

    if valid_data:
        result = AssignWcsStep.call(model, override_wavelengthrange=mock_reffile)
        assert result.meta.cal_step.assign_wcs == "COMPLETE"
    else:
        with pytest.raises(NoDataOnDetectorError, match="NRS2"):
            AssignWcsStep.call(model, override_wavelengthrange=mock_reffile)
