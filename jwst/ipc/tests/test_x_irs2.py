"""Test the x_irs2 module in the ipc package."""

from collections import namedtuple
import re

import numpy as np
import pytest

from jwst import datamodels
from jwst.ipc import x_irs2

from jwst.ipc.x_irs2 import ReadoutParam


def test_get_irs2_parameters_ndarray():
    # Numpy.ndarray as first parameter, always returns the same answer
    input_model = np.zeros((100, 100), dtype=np.float32)
    param = x_irs2._get_irs2_parameters(input_model)
    assert param == ReadoutParam(refout=640, n=16, r=4)


def test_get_irs2_parameters_blank_rampmodel():
    # Blank RampModel returns nrs_normal and nrs_reference = None
    # As does non-IRS2 data
    # So raises exception trying to combine NoneType with int
    input_model = datamodels.RampModel()
    with pytest.raises(
        TypeError, match=re.escape("unsupported operand type(s) for //: 'int' and 'NoneType'")
    ):
        param = x_irs2._get_irs2_parameters(input_model)


def test_get_irs2_parameters_non_default():
    # With some non-default values inserted into the irs2 attributes
    input_model = datamodels.RampModel()
    input_model.meta.exposure.nrs_normal = 8
    input_model.meta.exposure.nrs_reference = 2
    param = x_irs2._get_irs2_parameters(input_model)
    assert param == ReadoutParam(640, 8, 2)


def test_get_irs2_parameters_override_both():
    # Override datamodel attributes
    input_model = datamodels.RampModel()
    input_model.meta.exposure.nrs_normal = 8
    input_model.meta.exposure.nrs_reference = 2
    param = x_irs2._get_irs2_parameters(input_model, n=12, r=6)
    assert param == ReadoutParam(refout=764, n=12, r=6)


def test_get_irs2_parameters_remove_one():
    # Remove one of the IRS2 attributes
    input_model = datamodels.RampModel()
    input_model.meta.exposure.nrs_normal = 8
    input_model.meta.exposure.nrs_reference = 2
    del input_model.meta.exposure.nrs_normal
    # Even if you do this, accessing nrs_normal won't throw an AttributeError, it will return None
    # So the TypeError exception is thrown like above
    with pytest.raises(
        TypeError, match=re.escape("unsupported operand type(s) for //: 'int' and 'NoneType'")
    ):
        param = x_irs2._get_irs2_parameters(input_model)


def test_get_irs2_parameters_override_one():
    # Keep nrs_reference=2 in model attribute, use n=16 in function call
    input_model = datamodels.RampModel()
    input_model.meta.exposure.nrs_normal = 8
    input_model.meta.exposure.nrs_reference = 2
    del input_model.meta.exposure.nrs_normal
    param = x_irs2._get_irs2_parameters(input_model, n=16)
    assert param == ReadoutParam(refout=576, n=16, r=2)


def test_get_irs2_parameters_override_empty_values():
    input_model = datamodels.RampModel()
    input_model.meta.exposure.nrs_normal = 8
    input_model.meta.exposure.nrs_reference = 2
    del input_model.meta.exposure.nrs_normal
    del input_model.meta.exposure.nrs_reference
    # Now there's no parameters in the model, so both values are provided in function call
    param = x_irs2._get_irs2_parameters(input_model, n=12, r=6)
    assert param == ReadoutParam(refout=764, n=12, r=6)


def test_normal_shape_ndarray():
    # Input_model as numpy.ndarray, not IRS2 data
    input_model = np.zeros((512, 512), dtype=np.float32)
    assert (512, 512) == x_irs2.normal_shape(input_model)


def test_normal_shape_irs2_ndarray():
    # IRS2 sized numpy.ndarray
    input_model = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    assert (1, 1, 2048, 2048) == x_irs2.normal_shape(input_model, detector="NRS1")


def test_normal_shape_ndarray_irs2_nodetector():
    # If detector is not specified, will use the last dimension as the one with the reference pixels
    # For data in detector coordinates
    input_model = np.zeros((1, 1, 2048, 3200), dtype=np.float32)
    assert (1, 1, 2048, 2048) == x_irs2.normal_shape(input_model)


def test_normal_shape_model_irs2():
    # Simplest IRS2 datamodel with required attributes
    input_model = datamodels.RampModel()
    input_model.meta.instrument.detector = "NRS2"
    input_model.data = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    input_model.meta.exposure.nrs_normal = 16
    input_model.meta.exposure.nrs_reference = 4
    assert (1, 1, 2048, 2048) == x_irs2.normal_shape(input_model)


def test_normal_shape_model_not_nrs():
    # If detector is not NRS1 or NRS2, will raise an exception
    input_model = datamodels.RampModel()
    input_model.data = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    input_model.meta.exposure.nrs_normal = 16
    input_model.meta.exposure.nrs_reference = 4
    input_model.meta.instrument.detector = "NRCA1"
    with pytest.raises(RuntimeError, match="Detector NRCA1 is not supported for IRS2 data"):
        returned_shape = x_irs2.normal_shape(input_model)


def test_make_mask_default():
    input_model = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    mask = x_irs2.make_mask(input_model)
    # First 640 elements are reference pixels
    assert ~np.all(mask[:640])
    # 641-648 are science pixels
    assert np.all(mask[640:648])


def test_make_mask_branch():
    # n=13 traverses a different branch
    input_model = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    mask = x_irs2.make_mask(input_model, n=13)
    # The first 640 are still reference pixels
    assert ~np.all(mask[:640])


def test_from_irs2_detector():
    # Detector format data (3200 coluns)
    irs2_array = np.zeros((1, 1, 2048, 3200), dtype=np.float32)
    mask = x_irs2.make_mask(irs2_array)
    irs2_array += 1000.0
    irs2_array[:, :, :, ~mask] = 100.0
    result = x_irs2.from_irs2(irs2_array, mask, detector=None)
    assert np.allclose(result, 1000.0, rtol=1.0e-7)


def test_from_irs2_nrs1():
    # NRS1
    irs2_array = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    mask = x_irs2.make_mask(irs2_array)
    irs2_array += 1000.0
    irs2_array[:, :, ~mask, :] = 100.0
    result = x_irs2.from_irs2(irs2_array, mask, detector="NRS1")
    assert np.allclose(result, 1000.0, rtol=1.0e-7)


def test_from_irs2_nrs2():
    # NRS2
    irs2_array = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    mask = x_irs2.make_mask(irs2_array)
    irs2_array += 1000.0
    irs2_array[:, :, ~mask[::-1], :] = 100.0
    result = x_irs2.from_irs2(irs2_array, mask, detector="NRS2")
    assert np.allclose(result, 1000.0, rtol=1.0e-7)


def test_from_irs2_not_nrs():
    # Not NRS
    irs2_array = np.zeros((1, 1, 2048, 2048), dtype=np.float32)
    mask = x_irs2.make_mask(irs2_array)
    with pytest.raises(RuntimeError, match="Detector NRCA3 is not supported for IRS2 data"):
        result = x_irs2.from_irs2(irs2_array, mask, detector="NRCA3")


def test_to_irs2_detector():
    # Detector frame data, use detector=None
    # Last 8 columns should be data from normal_array
    irs2_array = np.zeros((1, 1, 2048, 3200), dtype=np.float32)
    normal_array = np.zeros((1, 1, 2048, 2048), dtype=np.float32) + 1000.0
    mask = x_irs2.make_mask(irs2_array)
    x_irs2.to_irs2(irs2_array, normal_array, mask, detector=None)
    assert np.allclose(irs2_array[:, :, :, 3192:], 1000.0, rtol=1.0e-7)


def test_to_irs2_nrs1():
    # NRS1 - last 8 rows should be data from normal array
    irs2_array = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    normal_array = np.zeros((1, 1, 2048, 2048), dtype=np.float32) + 1000.0
    mask = x_irs2.make_mask(irs2_array)
    x_irs2.to_irs2(irs2_array, normal_array, mask, detector="NRS1")
    assert np.allclose(irs2_array[:, :, 3192:], 1000.0, rtol=1.0e-7)


def test_to_irs2_nrs2():
    # NRS2 - first 8 rows should be data from normal array
    irs2_array = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    normal_array = np.zeros((1, 1, 2048, 2048), dtype=np.float32) + 1000.0
    mask = x_irs2.make_mask(irs2_array)
    x_irs2.to_irs2(irs2_array, normal_array, mask, detector="NRS2")
    assert np.allclose(irs2_array[:, :, :8], 1000.0, rtol=1.0e-7)


def test_to_irs2_not_nrs():
    # Non NRS2: raise exception
    irs2_array = np.zeros((1, 1, 3200, 2048), dtype=np.float32)
    normal_array = np.zeros((1, 1, 2048, 2048), dtype=np.float32) + 1000.0
    mask = x_irs2.make_mask(irs2_array)
    with pytest.raises(RuntimeError, match="Detector NRCA3 is not supported for IRS2 data"):
        result = x_irs2.to_irs2(irs2_array, normal_array, mask, detector="NRCA3")
