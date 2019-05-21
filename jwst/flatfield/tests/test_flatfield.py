import pytest
import numpy as np

from jwst import datamodels
from jwst.flatfield import FlatFieldStep


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("NIRCAM", "NRC_IMAGE"),
        ("NIRCAM", "NRC_WFSS"),
        ("NIRSPEC", "NRS_FOCUS"),
        ("MIRI", "MIR_IMAGE"),
        ("MIRI", "MIR_LRS-FIXEDSLIT"),
        ("MIRI", "MIR_MRS"),
        ("NIRISS", "NIS_IMAGE"),
        ("FGS", "FGS_IMAGE"),
    ]
)
def test_flatfield_step_interface(tmpdir, instrument, exptype):
    """Test that the basic inferface works for data requiring a FLAT reffile"""
    flat_file = str(tmpdir.join('flat.fits'))

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
    flat.data[0,0] = np.nan
    flat.err = np.random.random(shape) * 0.05
    flat.save(flat_file)

    # override class attribute so only the `flat` type needs to be overriden
    # in the step call.  Otherwise CRDS calls will be made for the other 3
    # types of flat reference file.
    FlatFieldStep.reference_file_types = ["flat"]
    result = FlatFieldStep.call(data, override_flat=flat_file)

    assert (result.data == data.data).all()
    assert result.var_flat.shape == shape
