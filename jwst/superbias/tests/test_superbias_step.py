import numpy as np

from jwst.superbias.superbias_step import SuperBiasStep
from jwst.superbias.tests.helpers import nrc_full_ramp


def test_full_step():
    """Test full run of the SuperBiasStep."""

    # Create inputs, data, and superbiases
    ngroups = 5
    nrows = 2048
    ncols = 2048

    data = nrc_full_ramp(ngroups, nrows, ncols)

    # Add signal values and bias values
    # Use signal = 0 ADU so value will be negative after superbias step
    data.data[0, :, 500, 500] = 0

    # Run the pipeline
    output = SuperBiasStep.call(data)

    # Check that pixel value is negative after bias is subtracted
    assert np.sign(output.data[0, 0, 500, 500]) == -1

    # Step is marked complete
    assert output.meta.cal_step.superbias == "COMPLETE"

    # Input is not modified
    assert output is not data
    assert data.meta.cal_step.superbias is None


def test_missing_reffile():
    """Test that SuperBiasStep is skipped for missing reference file."""

    # Create inputs, data, and superbiases
    ngroups = 5
    nrows = 20
    ncols = 20
    data = nrc_full_ramp(ngroups, nrows, ncols)

    # Run the pipeline
    output = SuperBiasStep.call(data, override_superbias="N/A")

    # Step is marked skipped
    assert output.meta.cal_step.superbias == "SKIPPED"

    # Input is not modified
    assert output is not data
    assert data.meta.cal_step.superbias is None
