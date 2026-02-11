"""Test the AssignWCSStep."""

import numpy as np
from gwcs import coordinate_frames as cf
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_miri import create_hdul as create_miri
from jwst.assign_wcs.tests.test_nircam import create_hdul as create_nircam
from jwst.assign_wcs.tests.test_niriss import create_hdul as create_niriss
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_ifu_file


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
