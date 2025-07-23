"""Test the AssignWCSStep."""

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_miri import create_hdul as create_miri
from jwst.assign_wcs.tests.test_nircam import create_hdul as create_nircam
from jwst.assign_wcs.tests.test_niriss import create_hdul as create_niriss


def test_assign_wcs_step_miri_ifu():
    hdul = create_miri(detector="MIRIFULONG", channel="34", band="MEDIUM")
    model = datamodels.CubeModel(hdul)
    model.data = np.zeros((3, 40, 50))
    result = AssignWcsStep.call(model)
    assert result is not model
    assert result.meta.cal_step.assign_wcs == "COMPLETE"


def test_assign_wcs_step_nis_wfss():
    hdul = create_niriss(filtername="GR150R", pupil="F200W", exptype="NIS_WFSS")
    model = datamodels.ImageModel(hdul)
    result = AssignWcsStep.call(model)
    assert result is not model
    assert result.meta.cal_step.assign_wcs == "COMPLETE"


def test_assign_wcs_step_nrc_wfss():
    hdul = create_nircam(exptype="NRC_WFSS", filtername="F444W", pupil="GRISMR")
    model = datamodels.ImageModel(hdul)
    result = AssignWcsStep.call(model)
    assert result is not model
    assert result.meta.cal_step.assign_wcs == "COMPLETE"


def test_unsupported_input(caplog):
    model = datamodels.SlitModel()
    result = AssignWcsStep.call(model)
    assert result is not model
    assert "type is not supported" in caplog.text
    assert result.meta.cal_step.assign_wcs == "SKIPPED"
