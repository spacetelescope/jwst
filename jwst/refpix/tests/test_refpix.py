import numpy as np
import pytest

from stdatamodels.jwst.datamodels import RampModel, dqflags

from jwst.refpix import RefPixStep
from jwst.refpix.reference_pixels import (
    Dataset,
    NIRDataset,
    correct_model,
    create_dataset,
    NRS_edgeless_subarrays,
)


conv_kernel_params = {
    "refpix_algorithm": "median",
    "sirs_kernel_model": None,
    "sigreject": 4,
    "gaussmooth": 1,
    "halfwidth": 30,
}


def test_refpix_subarray_miri():
    """
    Check that the correction is skipped for MIR subarray data.

    For MIRI, no reference pixel correction is performed on subarray data
    No changes should be seen in the data arrays before and after correction

    """
    ngroups = 3
    ysize = 22
    xsize = 28

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize)
    im.meta.subarray.name = "MASK1550"
    im.meta.subarray.ystart = 467

    # set reference pixel values left side
    im.data[:, 1:, :, 0] = 1.0
    im.data[:, 1:, :, 1] = 2.0
    im.data[:, 1:, :, 2] = 3.0
    im.data[:, 1:, :, 3] = 4.0

    outim = RefPixStep.call(im)

    # test that the science data are not changed
    np.testing.assert_array_equal(im.data, outim.data)


@pytest.mark.parametrize("subarray,ysize,xsize", [("SUB512", 32, 512), ("SUBS200A1", 64, 2048)])
def test_refpix_subarray_nirspec(subarray, ysize, xsize):
    """
    Check that the correction is performed for NRS subarray data.

    For NIRSpec, reference pixel correction is performed for all
    subarrays, with and without edges.

    Parameters
    ----------
    subarray : str
        Subarray name
    ysize : int
        Number of rows in input data
    xsize : int
        Number of columns in input data
    """
    # create input data
    ngroups = 3
    im = make_rampmodel(ngroups, ysize, xsize, instrument="NIRSPEC", fill_value=0.0)
    im.meta.subarray.name = subarray
    im.meta.subarray.xstart = 1
    im.meta.subarray.ystart = 1

    # set reference pixel values top and bottom, left and right
    im.data[:, :, :4, :] = 1.0
    im.data[:, :, -4:] = 2.0
    im.data[:, :, :, :4] = 3.0
    im.data[:, :, :, -4:] = 4.0

    # set reference pixels to "REFERENCE_PIXEL"
    if subarray not in NRS_edgeless_subarrays:
        im.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
        im.pixeldq[:, -4:] = dqflags.pixel["REFERENCE_PIXEL"]

    # run the step
    out = RefPixStep.call(im)

    # value subtracted should be the average of the left and right
    # reference pixels
    assert np.allclose(out.data[:, :, 4:-5, 4:-5], -3.5)


def test_each_amp():
    """
    Test that each amp is calculated separately.

    Using the average of left and right pixels.
    """
    ngroups = 7
    ysize = 1024
    xsize = 1032

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize)

    # set reference pixel values left and right side
    im.data[:, 1:, :, 0] = 1.0
    im.data[:, 1:, :, 1] = 2.0
    im.data[:, 1:, :, 2] = 3.0
    im.data[:, 1:, :, 3] = 4.0
    im.data[:, 1:, :, 1028] = 1.0
    im.data[:, 1:, :, 1029] = 2.0
    im.data[:, 1:, :, 1030] = 3.0
    im.data[:, 1:, :, 1031] = 4.0

    # set reference pixels to "REFERENCE_PIXEL"
    im.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
    im.pixeldq[:, 1028:] = dqflags.pixel["REFERENCE_PIXEL"]

    # run the step
    out = RefPixStep.call(im)

    # for amp 1, value subtracted should be 1, for amp 2, value should be 2, etc.
    assert out.data[0, 5, 100, 4] == 4.0  # pick a random pixel in the 4th column
    assert out.data[0, 5, 100, 5] == 3.0
    assert out.data[0, 5, 100, 6] == 2.0
    assert out.data[0, 5, 100, 7] == 1.0


def test_firstframe_sub():
    """
    Check that MIRI data has first group subtracted first.

    For MIR data, check that the first group is subtracted from each group in an integration
    and added back in after the correction.

    This was found in testing the amp step. Make sure that the first frame is
    subtracted from each group and added back in afterwards. If the reference pixels
    in the first group match the reference pixels in all other groups, then the
    subtraction will result in zeros, leaving zeros to be calculated as the reference
    pixel values, and the output data will match the input data after the frame is
    added back in. So there should be no change to the data.
    """
    ngroups = 5
    ysize = 1024
    xsize = 1032

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize)

    # set reference pixel values left and right side
    im.data[:, :, :, 0] = 1.0
    im.data[:, :, :, 1] = 2.0
    im.data[:, :, :, 2] = 3.0
    im.data[:, :, :, 3] = 4.0
    im.data[:, :, :, 1028] = 1.0
    im.data[:, :, :, 1029] = 2.0
    im.data[:, :, :, 1030] = 3.0
    im.data[:, :, :, 1031] = 4.0

    # set reference pixels to "REFERENCE_PIXEL"
    im.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
    im.pixeldq[:, 1028:] = dqflags.pixel["REFERENCE_PIXEL"]

    # run the step
    outim = RefPixStep.call(im)

    # test that the science data are not changed
    np.testing.assert_array_equal(im.data, outim.data)


def test_odd_even():
    """
    Check that odd/even rows are applied when flag is set.

    Test that odd and even rows are calculated separately.
    """
    ngroups = 7
    ysize = 1024
    xsize = 1032

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize)

    im.data = im.data * 10

    # set reference pixel values left and right side, odd and even rows.
    im.data[:, 1:, 1:ysize:2, 0] = 1.0
    im.data[:, 1:, 1:ysize:2, 1] = 2.0
    im.data[:, 1:, 1:ysize:2, 2] = 3.0
    im.data[:, 1:, 1:ysize:2, 3] = 4.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 0] = 5.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1] = 6.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 2] = 7.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 3] = 8.0
    im.data[:, 1:, 1:ysize:2, 1028] = 1.0
    im.data[:, 1:, 1:ysize:2, 1029] = 2.0
    im.data[:, 1:, 1:ysize:2, 1030] = 3.0
    im.data[:, 1:, 1:ysize:2, 1031] = 4.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1028] = 5.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1029] = 6.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1030] = 7.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1031] = 8.0

    # set reference pixels to "REFERENCE_PIXEL"
    im.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
    im.pixeldq[:, 1028:] = dqflags.pixel["REFERENCE_PIXEL"]

    # run the step
    out = RefPixStep.call(im)

    # values should be different by amp and by odd/even row
    # value of data in 5th frame is 50, ref values are subtracted from that
    assert out.data[0, 5, 100, 4] == 45.0  # pick a random pixel in the 4th column
    assert out.data[0, 5, 100, 5] == 44.0
    assert out.data[0, 5, 100, 6] == 43.0
    assert out.data[0, 5, 100, 7] == 42.0
    assert out.data[0, 5, 101, 4] == 49.0
    assert out.data[0, 5, 101, 5] == 48.0
    assert out.data[0, 5, 101, 6] == 47.0
    assert out.data[0, 5, 101, 7] == 46.0


@pytest.mark.parametrize(
    "detector,ysize,odd_even",
    [
        ("NRS1", 2048, True),
        ("NRS1", 3200, True),
        ("NRS1", 2048, False),
        ("NRS1", 3200, False),
        ("NRS2", 2048, True),
        ("NRS2", 3200, True),
        ("NRS2", 2048, False),
        ("NRS2", 3200, False),
    ],
)
def test_odd_even_amp_nirspec(detector, ysize, odd_even):
    """
    Check that odd/even columns are applied when flag is set.

    Test that odd and even rows are calculated separately

    Parameters
    ----------
    detector : str
        Name of detector
    ysize : int
        Number of rows in input data
    odd_even : bool
        Calculate odd and even column reference values separately?
    """
    ngroups = 1
    xsize = 2048

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize, instrument="NIRSPEC", fill_value=0.0)
    im.meta.instrument.detector = detector

    # check for irs2 data
    if ysize == 3200:
        n_amp = 5
        is_irs = True
    else:
        n_amp = 4
        is_irs = False
    amp_size = ysize // n_amp

    # set reference pixel values left and right side, odd and even rows.
    rval = []
    for i in range(n_amp):
        start_y = i * amp_size
        end_y = start_y + amp_size
        amp_val = 10 * i + 1.0
        im.data[:, :, start_y:end_y:2, :4] = amp_val
        im.data[:, :, start_y + 1 : end_y : 2, :4] = amp_val + 1
        im.data[:, :, start_y:end_y:2, -4:] = amp_val
        im.data[:, :, start_y + 1 : end_y : 2, -4:] = amp_val + 1

        rval.append(amp_val)

    # set reference pixels to "REFERENCE_PIXEL"
    im.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
    im.pixeldq[:, -4:] = dqflags.pixel["REFERENCE_PIXEL"]

    # run the step
    out = RefPixStep.call(
        im, use_side_ref_pixels=False, odd_even_columns=odd_even, irs2_mean_subtraction=True
    )

    # values should be different by amp and by odd/even row if specified
    # pick a random pixel to test
    # Note: only 4 amps in output, regardless of how many are in input
    test_row = 256
    test_col = xsize // 2
    for i in range(n_amp):
        if n_amp == 5:
            # skip ref section: first for nrs1, last for nrs2
            if i == 0 and detector == "NRS1":
                continue
            elif i == n_amp - 1 and detector == "NRS2":
                continue
        if odd_even:
            if is_irs:
                # odd/even is a little different for irs2: values get mixed
                assert np.isclose(out.data[0, 0, test_row, test_col], -rval[i] - 0.1)
                assert np.isclose(out.data[0, 0, test_row + 1, test_col], -rval[i] - 0.9)
            else:
                assert out.data[0, 0, test_row, test_col] == -rval[i]
                assert out.data[0, 0, test_row + 1, test_col] == -rval[i] - 1
        else:
            assert out.data[0, 0, test_row, test_col] == -rval[i] - 0.5
            assert out.data[0, 0, test_row + 1, test_col] == -rval[i] - 0.5

        test_row += 2048 // 4


def test_no_odd_even():
    """
    Check that odd/even rows are not applied if flag is set to False.

    Test that odd and even rows are calculated together
    """
    ngroups = 7
    ysize = 1024
    xsize = 1032

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize)

    im.data = im.data * 10

    # set reference pixel values left and right side, odd and even rows.
    im.data[:, 1:, 1:ysize:2, 0] = 1.0
    im.data[:, 1:, 1:ysize:2, 1] = 2.0
    im.data[:, 1:, 1:ysize:2, 2] = 3.0
    im.data[:, 1:, 1:ysize:2, 3] = 4.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 0] = 5.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1] = 6.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 2] = 7.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 3] = 8.0
    im.data[:, 1:, 1:ysize:2, 1028] = 1.0
    im.data[:, 1:, 1:ysize:2, 1029] = 2.0
    im.data[:, 1:, 1:ysize:2, 1030] = 3.0
    im.data[:, 1:, 1:ysize:2, 1031] = 4.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1028] = 5.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1029] = 6.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1030] = 7.0
    im.data[:, 1:, 0 : ysize - 1 : 2, 1031] = 8.0

    # set reference pixels to "REFERENCE_PIXEL"
    im.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
    im.pixeldq[:, 1028:] = dqflags.pixel["REFERENCE_PIXEL"]

    # run the step
    out = RefPixStep.call(im, odd_even_rows=False)

    # values should be different by amp and not by odd/even row
    # value of data in 5th frame is 50, ref values are subtracted from that
    # odd+even/2 -> (1+5)/2=3, (2+6)/2=4, (3+7)/2=5, (4+8)/2=6
    assert out.data[0, 5, 100, 4] == 47.0  # pick a random pixel in the 4th column
    assert out.data[0, 5, 100, 5] == 46.0
    assert out.data[0, 5, 100, 6] == 45.0
    assert out.data[0, 5, 100, 7] == 44.0
    assert out.data[0, 5, 101, 4] == 47.0
    assert out.data[0, 5, 101, 5] == 46.0
    assert out.data[0, 5, 101, 6] == 45.0
    assert out.data[0, 5, 101, 7] == 44.0


def test_side_averaging():
    """
    Test MIRI side averaging operation.

    For MIRI data, check that the mean value in the reference pixels is calculated for each
    amplifier using the average of the left and right side reference pixels.

    Test that the left and right side pixels are averaged.
    """
    ngroups = 7
    ysize = 1024
    xsize = 1032

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize)

    im.data = im.data * 10

    # set reference pixel values left and right side, odd and even rows.
    im.data[:, 1:, :, :4] = 1.0
    im.data[:, 1:, :, 1028:] = 2.0

    # set reference pixels to "REFERENCE_PIXEL"
    im.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
    im.pixeldq[:, 1028:] = dqflags.pixel["REFERENCE_PIXEL"]

    # run the step
    out = RefPixStep.call(im)

    # average reference pixel value should be 1.5 (all 1's on left, all 2's on right)
    assert out.data[0, 5, 100, 50] == 48.5


def test_above_sigma():
    """
    Test sigma rejection operation.

    Test that a value greater than 3 sigma above mean of reference pixels is rejected
    in the averaging of the reference pixels to be subtracted.
    """
    ngroups = 5
    ysize = 1024
    xsize = 1032

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize)

    im.data = im.data * 10

    # set reference pixel values left and right side, odd and even rows.
    im.data[:, 1:, :, :4] = 1.0
    im.data[:, 1:, :, 1028:] = 2.0
    im.data[0, 3, 50, 3] = 35.0

    # set reference pixels to "REFERENCE_PIXEL"
    im.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
    im.pixeldq[:, 1028:] = dqflags.pixel["REFERENCE_PIXEL"]

    # run the step
    out = RefPixStep.call(im)

    # average reference pixel value should be 1.5 (all 1's on left, all 2's on right)
    assert out.data[0, 3, 50, 7] == 28.5


def test_nan_refpix():
    """
    Verify that the reference pixels flagged DO_NOT_USE are not used in the calculation.

    Test that flagging a reference pixel with DO_NOT_USE does not use the pixel in the
    average. Set the pixel to NaN, which results in a NaN average value if used. If the test
    passes, then the NaN was correctly flagged and rejected from the average.
    """
    ngroups = 5
    ysize = 1024
    xsize = 1032

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize)

    im.data = im.data * 10

    # set reference pixel values left and right side, odd and even rows.
    im.data[:, 1:, :, :4] = 1.0
    im.data[:, 1:, :, 1028:] = 2.0
    im.data[0, 3, 50, 3] = np.nan

    # set reference pixels to "REFERENCE_PIXEL"
    im.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
    im.pixeldq[:, 1028:] = dqflags.pixel["REFERENCE_PIXEL"]
    im.pixeldq[50, 3] = dqflags.pixel["DO_NOT_USE"]

    # run the step
    out = RefPixStep.call(im)

    # average reference pixel value should be 1.5 (all 1's on left, all 2's on right)
    assert out.data[0, 3, 50, 7] == 28.5


def test_do_corrections_subarray_no_oddEven(setup_subarray_cube):
    """Test all corrections for subarray data with no even/odd."""

    # Create inputs and subarray SUB320A335R data, and set correction parameters
    ngroups = 3
    nrows = 160
    ncols = 160
    xstart = 1
    ystart = 1

    odd_even_columns = False
    use_side_ref_pixels = False
    side_smoothing_length = 11
    side_gain = 1.0
    odd_even_rows = False

    left_rpix = 5
    bottom_rpix = 7
    dataval = 150
    rmean = np.mean([left_rpix, bottom_rpix])

    input_model = setup_subarray_cube("SUB160", "NRCB1", xstart, ystart, ngroups, nrows, ncols)
    input_model.data[0, 0, :, :] = dataval
    input_model.data[0, 0, :4, :] = bottom_rpix
    input_model.data[0, 0, :, :4] = left_rpix
    input_model.pixeldq[:4, :] = dqflags.pixel["REFERENCE_PIXEL"]
    input_model.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]

    init_dataset = create_dataset(
        input_model,
        odd_even_columns,
        use_side_ref_pixels,
        side_smoothing_length,
        side_gain,
        odd_even_rows,
        conv_kernel_params,
    )

    init_dataset.do_corrections()

    np.testing.assert_almost_equal(
        np.mean(input_model.data[0, 0, :4, 4:-4]), bottom_rpix - rmean, decimal=0
    )
    np.testing.assert_almost_equal(
        np.mean(input_model.data[0, 0, 4:-4, :4]), left_rpix - rmean, decimal=0
    )
    np.testing.assert_almost_equal(
        np.mean(input_model.data[0, 0, 4:-4, 4:-4]), dataval - rmean, decimal=0
    )


def test_do_corrections_subarray(setup_subarray_cube):
    """Test all corrections for subarray data."""
    # Create inputs and subarray SUB320A335R data, and set correction parameters
    ngroups = 3
    nrows = 160
    ncols = 160
    xstart = 1
    ystart = 1

    odd_even_columns = True
    use_side_ref_pixels = False
    side_smoothing_length = 11
    side_gain = 1.0
    odd_even_rows = False

    left_rpix = 5
    bottom_rpix = 7
    dataval = 150
    rmean = np.mean([left_rpix, bottom_rpix])

    input_model = setup_subarray_cube("SUB160", "NRCB1", xstart, ystart, ngroups, nrows, ncols)
    input_model.data[0, 0, :, :] = dataval
    input_model.data[0, 0, :4, :] = bottom_rpix
    input_model.data[0, 0, :, :4] = left_rpix
    input_model.pixeldq[:4, :] = dqflags.pixel["REFERENCE_PIXEL"]
    input_model.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]

    init_dataset = create_dataset(
        input_model,
        odd_even_columns,
        use_side_ref_pixels,
        side_smoothing_length,
        side_gain,
        odd_even_rows,
        conv_kernel_params,
    )

    init_dataset.do_corrections()

    np.testing.assert_almost_equal(
        np.mean(input_model.data[0, 0, :4, 4:-4]), bottom_rpix - rmean, decimal=0
    )
    np.testing.assert_almost_equal(
        np.mean(input_model.data[0, 0, 4:-4, :4]), left_rpix - rmean, decimal=0
    )
    np.testing.assert_almost_equal(
        np.mean(input_model.data[0, 0, 4:-4, 4:-4]), dataval - rmean, decimal=0
    )


def test_do_corrections_subarray_4amp(setup_subarray_cube):
    """Test all corrections for subarray data."""
    # Create inputs and subarray SUBGRISM64 data, and set correction parameters
    ngroups = 3
    nrows = 64
    ncols = 2048
    xstart = 1
    ystart = 1

    odd_even_columns = True
    use_side_ref_pixels = True
    side_smoothing_length = 11
    side_gain = 1.0
    odd_even_rows = False

    left_rpix = 0
    right_rpix = 1
    side_rpix_mean = 0.5 * (left_rpix + right_rpix)
    bottom_rpix_a_odd = 7
    bottom_rpix_a_even = 7.4
    bottom_rpix_b_odd = 9
    bottom_rpix_b_even = 9.4
    bottom_rpix_c_odd = 6
    bottom_rpix_c_even = 6.4
    bottom_rpix_d_odd = 8
    bottom_rpix_d_even = 8.4
    dataval = 150

    input_model = setup_subarray_cube("SUBGRISM64", "NRCA1", xstart, ystart, ngroups, nrows, ncols)
    input_model.meta.exposure.noutputs = 4
    input_model.data[0, 0, 4:-4, 4:512:2] = dataval + bottom_rpix_a_odd + side_rpix_mean
    input_model.data[0, 0, 4:-4, 5:512:2] = dataval + bottom_rpix_a_even + side_rpix_mean
    input_model.data[0, 0, 4:-4, 512:1024:2] = dataval + bottom_rpix_b_odd + side_rpix_mean
    input_model.data[0, 0, 4:-4, 513:1024:2] = dataval + bottom_rpix_b_even + side_rpix_mean
    input_model.data[0, 0, 4:-4, 1024:1536:2] = dataval + bottom_rpix_c_odd + side_rpix_mean
    input_model.data[0, 0, 4:-4, 1025:1536:2] = dataval + bottom_rpix_c_even + side_rpix_mean
    input_model.data[0, 0, 4:-4, 1536:2044:2] = dataval + bottom_rpix_d_odd + side_rpix_mean
    input_model.data[0, 0, 4:-4, 1537:2044:2] = dataval + bottom_rpix_d_even + side_rpix_mean

    input_model.data[0, 0, :4, 0:512:2] = bottom_rpix_a_odd
    input_model.data[0, 0, :4, 1:512:2] = bottom_rpix_a_even
    input_model.data[0, 0, :4, 512:1024:2] = bottom_rpix_b_odd
    input_model.data[0, 0, :4, 513:1024:2] = bottom_rpix_b_even
    input_model.data[0, 0, :4, 1024:1536:2] = bottom_rpix_c_odd
    input_model.data[0, 0, :4, 1025:1536:2] = bottom_rpix_c_even
    input_model.data[0, 0, :4, 1536:2048:2] = bottom_rpix_d_odd
    input_model.data[0, 0, :4, 1537:2048:2] = bottom_rpix_d_even

    input_model.data[0, 0, :, :4:2] = left_rpix + bottom_rpix_a_odd
    input_model.data[0, 0, :, 1:4:2] = left_rpix + bottom_rpix_a_even
    input_model.data[0, 0, :, -4::2] = right_rpix + bottom_rpix_d_odd
    input_model.data[0, 0, :, -3::2] = right_rpix + bottom_rpix_d_even
    input_model.pixeldq[:4, :] = dqflags.pixel["REFERENCE_PIXEL"]
    input_model.pixeldq[:, :4] = dqflags.pixel["REFERENCE_PIXEL"]
    input_model.pixeldq[:, -4:] = dqflags.pixel["REFERENCE_PIXEL"]

    init_dataset = create_dataset(
        input_model,
        odd_even_columns,
        use_side_ref_pixels,
        side_smoothing_length,
        side_gain,
        odd_even_rows,
        conv_kernel_params,
    )

    init_dataset.do_corrections()

    np.testing.assert_almost_equal(np.mean(input_model.data[0, 0, 4:-4, 4:-4]), dataval, decimal=4)


def test_get_restore_group_subarray(setup_subarray_cube):
    """Test subarray input model data is replaced with group data."""

    # Create inputs and subarray SUB320A335R data, and set correction parameters
    ngroups = 3
    nrows = 320
    ncols = 320
    xstart = 486
    ystart = 1508

    odd_even_columns = True
    use_side_ref_pixels = True
    side_smoothing_length = 11
    side_gain = 1.0
    odd_even_rows = False

    input_model = setup_subarray_cube(
        "SUB320A335R", "NRCALONG", xstart, ystart, ngroups, nrows, ncols
    )
    input_model.data[0, 0, :, :] = 150

    init_dataset = Dataset(
        input_model,
        odd_even_columns,
        use_side_ref_pixels,
        side_smoothing_length,
        side_gain,
        conv_kernel_params,
        odd_even_rows,
    )

    # Make sure get_group properly copied the subarray
    assert np.all(init_dataset.input_model.data[0, 0, :, :] == 150)
    init_dataset.get_group(0, 0)
    assert np.shape(init_dataset.group.data) == (2048, 2048)
    assert np.all(init_dataset.group[: ystart - 1, : xstart - 1] == 0)
    assert np.all(init_dataset.group[ystart - 1 + nrows :, xstart - 1 + ncols :] == 0)
    assert np.all(
        init_dataset.group[ystart - 1 : ystart - 1 + nrows, xstart - 1 : xstart - 1 + ncols] == 150
    )

    init_dataset.group[:, :] = 20
    init_dataset.restore_group(0, 0)
    assert np.all(init_dataset.input_model.data[0, 0, :, :] == 20)


def test_do_top_bottom_correction(setup_cube):
    """Test top/bottom correction for NIRCam data."""

    ngroups = 3
    nrows = 2048
    ncols = 2048

    odd_even_columns = True
    use_side_ref_pixels = True
    side_smoothing_length = 11
    side_gain = 1.0

    input_model = setup_cube("NIRCAM", "NRCALONG", ngroups, nrows, ncols)
    input_model.meta.subarray.name = "FULL"
    init_dataset = NIRDataset(
        input_model,
        odd_even_columns,
        use_side_ref_pixels,
        side_smoothing_length,
        side_gain,
        conv_kernel_params,
    )

    abounds = [0, 512, 1024, 1536, 2048]
    top_even_amps = [12, 13, 14, 15]
    top_odd_amps = [16, 17, 18, 19]
    bottom_even_amps = [20, 21, 22, 23]
    bottom_odd_amps = [24, 25, 26, 27]
    dataval = [50, 51, 52, 53]

    for i in np.arange(0, len(abounds) - 1):
        # bottom, odd
        input_model.data[0, 0, :4, abounds[i] : abounds[i + 1] : 2] = bottom_even_amps[i]

        # bottom, even
        input_model.data[0, 0, :4, abounds[i] + 1 : abounds[i + 1] - 1 : 2] = bottom_odd_amps[i]

        # top, odd
        input_model.data[0, 0, -4:, abounds[i] : abounds[i + 1] : 2] = top_even_amps[i]

        # top, even
        input_model.data[0, 0, -4:, abounds[i] + 1 : abounds[i + 1] - 1 : 2] = top_odd_amps[i]

        # data
        input_model.data[0, 0, 4:-4, abounds[i] : abounds[i + 1]] = dataval[i]

    refpix = init_dataset.get_refvalues(input_model.data[0, 0, :, :])
    init_dataset.do_top_bottom_correction(input_model.data[0, 0, :, :], refpix)

    for i in np.arange(0, len(abounds) - 1):
        even_rmean = np.mean([bottom_even_amps[i], top_even_amps[i]])
        odd_rmean = np.mean([bottom_odd_amps[i], top_odd_amps[i]])
        rmean = np.mean([even_rmean, odd_rmean])
        np.testing.assert_almost_equal(
            np.mean(input_model.data[0, 0, :4, abounds[i] : abounds[i + 1] : 2]),
            bottom_even_amps[i] - even_rmean,
            decimal=1,
        )
        np.testing.assert_almost_equal(
            np.mean(input_model.data[0, 0, :4, abounds[i] + 1 : abounds[i + 1] : 2]),
            bottom_odd_amps[i] - odd_rmean,
            decimal=1,
        )
        np.testing.assert_almost_equal(
            np.mean(input_model.data[0, 0, -4:, abounds[i] : abounds[i + 1] : 2]),
            top_even_amps[i] - even_rmean,
            decimal=1,
        )
        np.testing.assert_almost_equal(
            np.mean(input_model.data[0, 0, -4:, abounds[i] + 1 : abounds[i + 1] : 2]),
            top_odd_amps[i] - odd_rmean,
            decimal=1,
        )
        np.testing.assert_almost_equal(
            np.mean(input_model.data[0, 0, 4:-4, abounds[i] : abounds[i + 1]]),
            dataval[i] - rmean,
            decimal=1,
        )


def test_do_top_bottom_correction_no_even_odd(setup_cube):
    """Test top/bottom correction with no even/odd."""

    ngroups = 3
    nrows = 2048
    ncols = 2048

    odd_even_columns = False
    use_side_ref_pixels = True
    side_smoothing_length = 11
    side_gain = 1.0

    input_model = setup_cube("NIRCAM", "NRCALONG", ngroups, nrows, ncols)
    input_model.meta.subarray.name = "FULL"
    init_dataset = NIRDataset(
        input_model,
        odd_even_columns,
        use_side_ref_pixels,
        side_smoothing_length,
        side_gain,
        conv_kernel_params,
    )

    abounds = [0, 512, 1024, 1536, 2048]
    top_amps = [12, 13, 14, 15]
    bottom_amps = [16, 17, 18, 19]
    dataval = [50, 51, 52, 53]

    for i in np.arange(0, len(abounds) - 1):
        # bottom
        input_model.data[0, 0, :4, abounds[i] : abounds[i + 1]] = bottom_amps[i]

        # top
        input_model.data[0, 0, -4:, abounds[i] : abounds[i + 1]] = top_amps[i]

        # data
        input_model.data[0, 0, 4:-4, abounds[i] : abounds[i + 1]] = dataval[i]

    refpix = init_dataset.get_refvalues(input_model.data[0, 0, :, :])
    init_dataset.do_top_bottom_correction(input_model.data[0, 0, :, :], refpix)

    for i in np.arange(0, len(abounds) - 1):
        rmean = np.mean([top_amps[i], bottom_amps[i]])
        np.testing.assert_almost_equal(
            np.mean(input_model.data[0, 0, :4, abounds[i] : abounds[i + 1]]),
            bottom_amps[i] - rmean,
            decimal=1,
        )
        np.testing.assert_almost_equal(
            np.mean(input_model.data[0, 0, -4:, abounds[i] : abounds[i + 1]]),
            top_amps[i] - rmean,
            decimal=1,
        )
        np.testing.assert_almost_equal(
            np.mean(input_model.data[0, 0, 4:-4, abounds[i] : abounds[i + 1]]),
            dataval[i] - rmean,
            decimal=1,
        )


def make_rampmodel(ngroups, ysize, xsize, instrument="MIRI", fill_value=None):
    """
    Make MIRI or NIRSpec ramp model for testing.

    Parameters
    ----------
    ngroups : int
        Number of groups in created data
    ysize : int
        Number of rows in created data
    xsize : int
        Number of columns in created data
    instrument : str
        Instrument name
    full_value : float or None
        Fill value for data
    Returns
    -------
    dm_ramp : jwst RampModel
        The jwst RampModel created
    """

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)

    # create JWST datamodel
    dm_ramp = RampModel(csize)

    # set each frame equal to frame number if fill value is not provided
    if fill_value is None:
        for i in range(0, ngroups - 1):
            dm_ramp.data[0, i, :, :] = i
    else:
        dm_ramp.data[:] = fill_value

    # populate header of data model
    if instrument == "NIRSPEC":
        dm_ramp.meta.instrument.name = "NIRSPEC"
        dm_ramp.meta.instrument.detector = "NRS1"
        dm_ramp.meta.exposure.type = "NRS_FIXEDSLIT"
        if ysize > 2048:
            dm_ramp.meta.exposure.readpatt = "NRSIRS2"
            dm_ramp.meta.exposure.nrs_normal = 16
            dm_ramp.meta.exposure.nrs_reference = 4
        else:
            dm_ramp.meta.exposure.readpatt = "NRS"
    else:
        dm_ramp.meta.instrument.name = "MIRI"
        dm_ramp.meta.instrument.detector = "MIRIMAGE"
        dm_ramp.meta.instrument.filter = "F560W"
        dm_ramp.meta.exposure.type = "MIR_IMAGE"

    dm_ramp.meta.instrument.band = "N/A"
    dm_ramp.meta.observation.date = "2016-06-01"
    dm_ramp.meta.observation.time = "00:00:00"
    dm_ramp.meta.subarray.name = "FULL"
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    return dm_ramp


@pytest.fixture(scope="function")
def setup_cube():
    """
    Set up fake data to test.

    Returns
    -------
    _cube : function
        The function that creates the datamodel
    """

    def _cube(instr, detector, ngroups, nrows, ncols):
        nints = 1

        # create a JWST datamodel for any instrument's FULL data
        data_model = RampModel((nints, ngroups, nrows, ncols))
        data_model.meta.subarray.name = "FULL"
        data_model.meta.subarray.xstart = 1
        data_model.meta.subarray.ystart = 1
        data_model.meta.subarray.xsize = ncols
        data_model.meta.subarray.ysize = nrows
        data_model.meta.exposure.ngroups = ngroups
        data_model.meta.instrument.name = instr
        data_model.meta.instrument.detector = detector
        data_model.meta.observation.date = "2019-10-14"
        data_model.meta.observation.time = "16:44:12.000"

        return data_model

    return _cube


@pytest.fixture(scope="function")
def setup_subarray_cube():
    """
    Set up fake NIRCam subarray data to test.

    Returns
    -------
    _cube : function
        The cube creating function
    """

    def _cube(name, detector, xstart, ystart, ngroups, nrows, ncols):
        nints = 1

        # create a JWST datamodel for NIRCam subarray data
        data_model = RampModel((nints, ngroups, nrows, ncols))
        data_model.meta.subarray.name = name
        data_model.meta.subarray.xstart = xstart
        data_model.meta.subarray.ystart = ystart
        data_model.meta.subarray.xsize = ncols
        data_model.meta.subarray.ysize = nrows
        data_model.meta.exposure.ngroups = ngroups
        data_model.meta.instrument.name = "NIRCAM"
        data_model.meta.instrument.detector = detector
        data_model.meta.observation.date = "2019-10-14"
        data_model.meta.observation.time = "16:44:12.000"

        return data_model

    return _cube


@pytest.mark.parametrize(
    "instr, det",
    [
        ("NIRCAM", "NRCA1"),
        ("NIRCAM", "NRCA2"),
        ("NIRCAM", "NRCA3"),
        ("NIRCAM", "NRCA4"),
        ("NIRCAM", "NRCALONG"),
        ("NIRCAM", "NRCB1"),
        ("NIRCAM", "NRCB2"),
        ("NIRCAM", "NRCB3"),
        ("NIRCAM", "NRCB4"),
        ("NIRCAM", "NRCBLONG"),
        ("FGS", "GUIDER1"),
        ("FGS", "GUIDER2"),
    ],
)
def test_correct_model(setup_cube, instr, det):
    """Test all corrections for full frame data for all detectors."""

    ngroups = 2
    nrows = 2048
    ncols = 2048
    # Set correction parameters
    odd_even_columns = True
    use_side_ref_pixels = True
    side_smoothing_length = 11
    side_gain = 1.0
    odd_even_rows = False

    rpix = 7
    dataval = 150

    input_model = setup_cube(instr, det, ngroups, nrows, ncols)
    input_model.data[0, 0, :, :] = rpix
    input_model.data[0, 0, 4:-4, 4:-4] = dataval

    correct_model(
        input_model,
        odd_even_columns,
        use_side_ref_pixels,
        side_smoothing_length,
        side_gain,
        odd_even_rows,
        conv_kernel_params,
    )

    np.testing.assert_almost_equal(np.mean(input_model.data[0, 0, :4, 4:-4]), 0, decimal=0)
    np.testing.assert_almost_equal(
        np.mean(input_model.data[0, 0, 4:-4, 4:-4]), dataval - rpix, decimal=0
    )


def test_zero_frame(setup_cube):
    """
    Test ZEROFRAME refpix processing.
    """

    ngroups = 2
    nrows = 2048
    ncols = 2048
    # Set correction parameters
    odd_even_columns = True
    use_side_ref_pixels = True
    side_smoothing_length = 11
    side_gain = 1.0
    odd_even_rows = False

    rpix = 7
    dataval = 150
    instr = "NIRCAM"
    det = "NRCA1"

    # Setup RampModel
    input_model = setup_cube(instr, det, ngroups, nrows, ncols)
    nints = input_model.shape[0]
    input_model.data[0, 0, :, :] = rpix
    input_model.data[0, 0, 4:-4, 4:-4] = dataval

    # Setup ZEROFRAME
    input_model.zeroframe = np.zeros((1, nrows, ncols), dtype=float)
    input_model.zeroframe[0, :, :] = rpix / 2.0
    input_model.zeroframe[0, 4:-4, 4:-4] = dataval / 2.0
    input_model.zeroframe[0, 5, 5] = 0.0  # Test a bad pixel.
    input_model.meta.exposure.zero_frame = True

    correct_model(
        input_model,
        odd_even_columns,
        use_side_ref_pixels,
        side_smoothing_length,
        side_gain,
        odd_even_rows,
        conv_kernel_params,
    )

    # Make sure the SCI data is as expected.
    data = np.zeros(input_model.data.shape, dtype=input_model.data.dtype)
    data[0, 0, 4:-4, 4:-4] = dataval - rpix
    np.testing.assert_almost_equal(input_model.data, data, decimal=5)

    # Check the ZEROFRAME
    zeroframe = np.zeros((nints, nrows, ncols), dtype=float)
    zeroframe[0, 4:-4, 4:-4] = dataval / 2.0 - rpix / 2.0
    zeroframe[0, 5, 5] = 0.0  # Make sure this pixel is zero.
    np.testing.assert_almost_equal(input_model.zeroframe, zeroframe, decimal=5)


@pytest.mark.parametrize(
    "detector,irs2,preserve",
    [
        ("NRS1", False, False),
        ("NRS1", False, True),
        ("NRS1", True, True),
        ("NRS1", True, False),
        ("NRS2", False, False),
        ("NRS2", False, True),
        ("NRS2", True, True),
        ("NRS2", True, False),
    ],
)
def test_preserve_refpix(detector, irs2, preserve):
    # make some nirspec data
    ngroups = 1
    xsize = 2048
    if irs2:
        ysize = 3200
    else:
        ysize = 2048

    # make ramp model
    im = make_rampmodel(ngroups, ysize, xsize, instrument="NIRSPEC", fill_value=0.0)
    im.meta.instrument.detector = detector

    # run the step
    out = RefPixStep.call(im, preserve_irs2_refpix=preserve)
    if not irs2:
        # parameter ignored for non-irs2 data
        assert out.data.shape == (1, ngroups, ysize, xsize)
        assert out.pixeldq.shape == (ysize, xsize)
    elif preserve:
        # output data shape is the same
        assert out.data.shape == (1, ngroups, ysize, xsize)
        assert out.pixeldq.shape == (ysize, xsize)
    else:
        # output data is trimmed to remove interleaved refpix
        assert out.data.shape == (1, ngroups, xsize, xsize)
        assert out.pixeldq.shape == (xsize, xsize)
