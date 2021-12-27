import pytest
import numpy as np

from jwst import datamodels
from jwst.flatfield import FlatFieldStep
from jwst.flatfield.flat_field_step import NRS_IMAGING_MODES, NRS_SPEC_MODES


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("NIRCAM", "NRC_IMAGE"),
        ("NIRCAM", "NRC_WFSS"),
        ("MIRI", "MIR_IMAGE"),
        ("MIRI", "MIR_LRS-FIXEDSLIT"),
        ("MIRI", "MIR_LRS-SLITLESS"),
        ("MIRI", "MIR_MRS"),
        ("NIRISS", "NIS_IMAGE"),
        ("NIRISS", "NIS_WFSS"),
        ("NIRISS", "NIS_SOSS"),
        ("NIRISS", "NIS_AMI"),
        ("FGS", "FGS_IMAGE"),
    ] + [("NIRSPEC", exptype) for exptype in NRS_IMAGING_MODES]
)
@pytest.mark.skip(reason="modifying reference_file_types caused other tests to fail")
def test_flatfield_step_interface(instrument, exptype):
    """Test that the basic interface works for data requiring a FLAT reffile"""

    shape = (20, 20)

    data = datamodels.ImageModel(shape)
    data.meta.instrument.name = instrument
    data.meta.exposure.type = exptype
    data.meta.subarray.xstart = 1
    data.meta.subarray.ystart = 1
    data.meta.subarray.xsize = shape[1]
    data.meta.subarray.ysize = shape[0]

    flat = datamodels.FlatModel(shape)
    flat.meta.instrument.name = instrument
    flat.meta.subarray.xstart = 1
    flat.meta.subarray.ystart = 1
    flat.meta.subarray.xsize = shape[1]
    flat.meta.subarray.ysize = shape[0]
    flat.data += 1
    flat.data[0, 0] = np.nan
    flat.err = np.random.random(shape) * 0.05

    # override class attribute so only the `flat` type needs to be overridden
    # in the step call.  Otherwise CRDS calls will be made for the other 3
    # types of flat reference file not used in this test.
    FlatFieldStep.reference_file_types = ["flat"]
    result = FlatFieldStep.call(data, override_flat=flat)

    assert (result.data == data.data).all()
    assert result.var_flat.shape == shape
    assert result.meta.cal_step.flat_field == 'COMPLETE'


def exptypes():
    """Generate NRS EXPTYPES from the schema enum, removing spec types"""
    model = datamodels.ImageModel()
    alltypes = set(model.meta.exposure._schema['properties']['type']['enum'])
    spectypes = set(NRS_SPEC_MODES)
    return sorted([i for i in (alltypes - spectypes)])


@pytest.mark.parametrize(
    "exptype",
    exptypes()
)
def test_nirspec_flatfield_step_interface(exptype):
    """Test that the interface works all NIRSpec types"""

    shape = (20, 20)

    data = datamodels.ImageModel(shape)
    data.meta.observation.date = "2019-01-01"
    data.meta.observation.time = "00:00:00"
    data.meta.instrument.name = "NIRSPEC"
    data.meta.instrument.detector = "NRS1"
    data.meta.instrument.filter = "CLEAR"
    data.meta.instrument.grating = "MIRROR"
    data.meta.exposure.type = exptype
    data.meta.subarray.xstart = 1
    data.meta.subarray.ystart = 1
    data.meta.subarray.xsize = shape[1]
    data.meta.subarray.ysize = shape[0]

    FlatFieldStep.call(data)
