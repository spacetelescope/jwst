import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.assign_wcs.tests.helpers import make_mock_dhs_nrca1_rate
from jwst.dq_init.tests.helpers import make_superstripe_mask_model, make_superstripe_model
from jwst.lib import stripe_utils
from jwst.lib.reffile_utils import science_detector_frame_transform


@pytest.fixture(scope="module")
def substripe_model():
    return make_mock_dhs_nrca1_rate()


@pytest.fixture(scope="module")
def superstripe_model():
    return make_superstripe_model()


@pytest.mark.parametrize("science_frame", [True, False])
def test_generate_substripe_ranges(substripe_model, science_frame):
    all_ranges = stripe_utils.generate_substripe_ranges(
        substripe_model, science_frame=science_frame
    )
    expected = {
        "full": {0: [1527, 1567], 1: [1652, 1692], 2: [1777, 1817], 3: [1902, 1942]},
        "subarray": {0: [1, 41], 1: [42, 82], 2: [83, 123], 3: [124, 164]},
        "reference_full": {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]},
        "reference_subarray": {0: [0, 1], 1: [41, 42], 2: [82, 83], 3: [123, 124]},
    }

    # Since slowaxis is positive, science and detector frames are the same
    assert all_ranges == expected


def test_generate_substripe_ranges_swap_slow(substripe_model):
    model = substripe_model.copy()
    model.meta.subarray.slowaxis *= -1

    # Detector orientation: same as for positive slowaxis
    all_ranges = stripe_utils.generate_substripe_ranges(model)
    expected = {
        "full": {0: [1527, 1567], 1: [1652, 1692], 2: [1777, 1817], 3: [1902, 1942]},
        "subarray": {0: [1, 41], 1: [42, 82], 2: [83, 123], 3: [124, 164]},
        "reference_full": {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]},
        "reference_subarray": {0: [0, 1], 1: [41, 42], 2: [82, 83], 3: [123, 124]},
    }
    assert all_ranges == expected

    # Since slowaxis is negative, science ranges will be swapped
    all_ranges = stripe_utils.generate_substripe_ranges(model, science_frame=True)
    expected = {
        "full": {0: [481, 521], 1: [356, 396], 2: [231, 271], 3: [106, 146]},
        "subarray": {0: [123, 163], 1: [82, 122], 2: [41, 81], 3: [0, 40]},
        "reference_full": {0: [2047, 2048], 1: [2047, 2048], 2: [2047, 2048], 3: [2047, 2048]},
        "reference_subarray": {0: [163, 164], 1: [122, 123], 2: [81, 82], 3: [40, 41]},
    }
    assert all_ranges == expected


def test_generate_substripe_ranges_invalid_input(substripe_model):
    model = substripe_model.copy()
    model.meta.subarray.multistripe_reads2 = 0

    with pytest.raises(ValueError, match="Invalid value for multistripe_reads2"):
        stripe_utils.generate_substripe_ranges(model)


def test_generate_substripe_ranges_wrong_size(substripe_model):
    model = substripe_model.copy()
    model.meta.subarray.multistripe_reads2 = 38

    with pytest.raises(ValueError, match="readout does not match science array shape"):
        stripe_utils.generate_substripe_ranges(model)


def test_generate_superstripe_ranges(superstripe_model):
    all_ranges = stripe_utils.generate_superstripe_ranges(superstripe_model)

    expected_full = {
        0: [(4, 208)],
        1: [(208, 412)],
        2: [(412, 616)],
        3: [(616, 820)],
        4: [(820, 1024)],
        5: [(1024, 1228)],
        6: [(1228, 1432)],
        7: [(1432, 1636)],
        8: [(1636, 1840)],
        9: [(1840, 2044)],
    }
    expected_sub = dict.fromkeys(range(10), [(4, 208)])
    expected_ref_full = dict.fromkeys(range(10), [(0, 4)])
    expected_ref_sub = dict.fromkeys(range(10), [(0, 4)])
    assert all_ranges["full"] == expected_full
    assert all_ranges["subarray"] == expected_sub
    assert all_ranges["reference_full"] == expected_ref_full
    assert all_ranges["reference_subarray"] == expected_ref_sub


def test_generate_superstripe_ranges_science_frame(superstripe_model):
    all_ranges = stripe_utils.generate_superstripe_ranges(superstripe_model, science_frame=True)

    # Swapped from the detector frame values
    expected_full = {
        0: [(1840, 2044)],
        1: [(1636, 1840)],
        2: [(1432, 1636)],
        3: [(1228, 1432)],
        4: [(1024, 1228)],
        5: [(820, 1024)],
        6: [(616, 820)],
        7: [(412, 616)],
        8: [(208, 412)],
        9: [(4, 208)],
    }
    expected_sub = dict.fromkeys(range(10), [(0, 204)])
    expected_ref_full = dict.fromkeys(range(10), [(2044, 2048)])
    expected_ref_sub = dict.fromkeys(range(10), [(204, 208)])
    assert all_ranges["full"] == expected_full
    assert all_ranges["subarray"] == expected_sub
    assert all_ranges["reference_full"] == expected_ref_full
    assert all_ranges["reference_subarray"] == expected_ref_sub


def test_generate_superstripe_ranges_invalid_input(superstripe_model):
    model = superstripe_model.copy()
    model.meta.subarray.multistripe_reads2 = 0

    with pytest.raises(ValueError, match="Invalid value for multistripe_reads2"):
        stripe_utils.generate_superstripe_ranges(model)


def test_generate_superstripe_ranges_wrong_size(superstripe_model):
    model = superstripe_model.copy()
    model.meta.subarray.multistripe_reads2 = 38

    with pytest.raises(ValueError, match="readout does not match science array shape"):
        stripe_utils.generate_superstripe_ranges(model)


def _assign_metadata(metanode, keys, vals):
    """
    Assign a list of values to parameters in a datamodel
    for ease of testing.
    Stripe params: xsize_sci, ysize_sci, nreads1, nreads2, nskips1,
                   nskips2, repeat_stripe, interleave_reads1, fastaxis, slowaxis
    """
    for key, val in zip(keys, vals):
        setattr(metanode.subarray, key, val)


def _multistripe_subarray_keys():
    """
    Returns the metadata keys that we'll assign values to for testing purposes.
    """
    return (
        "xsize",
        "ysize",
        "multistripe_reads1",
        "multistripe_reads2",
        "multistripe_skips1",
        "multistripe_skips2",
        "repeat_stripe",
        "interleave_reads1",
        "superstripe_step",
        "num_superstripe",
        "fastaxis",
        "slowaxis",
    )


def test_generate_stripe_reference_substripe():
    # Mock model for metadata
    model = datamodels.RampModel()
    sci_meta = model.meta
    subarray_keys = _multistripe_subarray_keys()

    # Generate test array with pixel values
    # equal to row number in detector frame.
    test_array = (np.ones((2048, 2048), dtype=int) * np.arange(2048)).T

    # SUB41STRIPE1_DHS nrca1 case
    stripe_params = (2048, 41, 1, 40, 1901, 0, 1, 1, 0, 0, -1, 2)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)

    # Function presumes input in science frame, so move test array to science frame
    # before supplying to function.
    stripe1_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe1_array.shape == (41, 2048)
    assert stripe1_array[0, 1024] == 0
    assert stripe1_array[1, 1024] == 1902  # nreads1 + nskips1

    # Test swapped axes
    # Note: If axes are swapped, substriped array will also be rotated, e.g.
    # xsize corresponds to slowaxis, such that ysize will be full detector size.
    stripe_params = (41, 2048, 1, 40, 1901, 0, 1, 1, 0, 0, 2, 1)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)
    stripe1swap_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe1swap_array.shape == (2048, 41)
    assert stripe1swap_array[1024, 1] == 1902  # nreads1 + nskips1

    # SUB82STRIPE2_DHS nrca2 case
    stripe_params = (2048, 82, 1, 40, 1662, 82, 1, 1, 0, 0, 1, -2)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)

    stripe2_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe2_array.shape == (82, 2048)
    # nrca2 has flipped row direction, so in science frame the row indices are flipped.
    assert stripe2_array[-1, 1024] == 0
    assert stripe2_array[-2, 1024] == 1663  # nreads1 + nskips1
    assert stripe2_array[-42, 1024] == 0
    assert stripe2_array[-43, 1024] == 1785  # nreads1 + nskips1 + nreads2 + nskips2

    # SUB164STRIPE4_DHS nrcalong case
    stripe_params = (2048, 164, 1, 40, 971, 0, 1, 0, 0, 0, -1, 2)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)

    stripe4_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe4_array.shape == (164, 2048)
    assert stripe4_array[0, 1024] == 0
    assert stripe4_array[1, 1024] == 972  # nreads1 + nskips1
    assert stripe4_array[42, 1024] == 972  # nreads1 + nskips1, stripe 2
    assert stripe4_array[83, 1024] == 972  # nreads1 + nskips1, stripe 3


def test_generate_stripe_reference_substripe_repeat_stripe_zero():
    """
    Test generate_stripe_reference in DHS mode when repeat_stripe=0.

    When repeat_stripe=0 the read head does not reset between stripes;
    successive nreads2 blocks are separated by nskips2 rows advancing
    monotonically through the detector.
    """
    # Array where every pixel value equals its row index.
    test_array = (np.ones((2048, 2048), dtype=int) * np.arange(2048)).T

    xsize_sci = 2048
    nreads1, nreads2 = 1, 20
    nskips1, nskips2 = 10, 5
    repeat_stripe = 0
    interleave_reads1 = 0
    superstripe_step, num_superstripe = 0, 0
    fastaxis, slowaxis, ngroups = -1, 2, 1
    # Detector read sequence: row 0 (nreads1), skip 10 (nskips1),
    # rows 11-30 (nreads2), skip 5 (nskips2), rows 36-55 (nreads2),
    # skip 5 (nskips2), rows 61-80 (nreads2) → 61 output rows total.
    ysize_sci = nreads1 + nreads2 * 3  # 1 + 20*3 = 61

    stripe_params = (
        xsize_sci,
        ysize_sci,
        nreads1,
        nreads2,
        nskips1,
        nskips2,
        repeat_stripe,
        interleave_reads1,
        superstripe_step,
        num_superstripe,
        fastaxis,
        slowaxis,
    )
    model = datamodels.RampModel()
    sci_meta = model.meta
    _assign_metadata(sci_meta, _multistripe_subarray_keys(), stripe_params)
    stripe_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe_array.shape == (ysize_sci, xsize_sci)
    # Output row 0: first nreads1 read from detector row 0.
    assert stripe_array[0, 1024] == 0
    # Output row 1: start of first nreads2 block; head advanced by nreads1 + nskips1.
    assert stripe_array[1, 1024] == nreads1 + nskips1
    # Output row nreads1 + nreads2: start of second nreads2 block;
    # head advanced by an additional nreads2 + nskips2 (no reset).
    assert stripe_array[nreads1 + nreads2, 1024] == nreads1 + nskips1 + nreads2 + nskips2
    # Output row nreads1 + 2*nreads2: start of third nreads2 block;
    # head advanced by yet another nreads2 + nskips2
    assert stripe_array[nreads1 + 2 * nreads2, 1024] == nreads1 + nskips1 + 2 * (nreads2 + nskips2)


def test_generate_stripe_reference_superstripe():
    # Mock model for metadata
    model = datamodels.RampModel()
    sci_meta = model.meta
    sci_meta.exposure.ngroups = 5
    subarray_keys = _multistripe_subarray_keys()

    # Generate test array with pixel values equal to row number in detector frame.
    test_array = np.ones((2048, 256), dtype=int) * np.arange(2048)[:, None]

    # SUB204STRIPE_SOSS case
    stripe_params = (208, 256, 4, 204, 0, 0, 1, 1, 204, 10, -2, -1)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)

    # Function presumes input in science frame, so move test array to science frame
    # before supplying to function.
    stripe1_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe1_array.shape == (10, 256, 208)
    # detector row numbers count backwards in columns, but sections increase, left to right
    np.testing.assert_allclose(stripe1_array[0, 0], np.arange(207, -1, -1))
    np.testing.assert_allclose(stripe1_array[-1, 0, :-4], np.arange(2043, 1839, -1))
    # reference pixels are always 3, 2, 1, 0
    assert np.all(stripe1_array[:, :, -4:] == np.arange(3, -1, -1))

    # Test swapped axes: striped array is rotated
    stripe_params = (256, 208, 4, 204, 0, 0, 1, 1, 204, 10, -1, -2)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)
    stripe1swap_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    # Dimensions are swapped, values are the same
    assert stripe1swap_array.shape == (10, 208, 256)
    np.testing.assert_allclose(stripe1swap_array[0, :, 0], np.arange(207, -1, -1))
    np.testing.assert_allclose(stripe1swap_array[-1, :-4, 0], np.arange(2043, 1839, -1))
    assert np.all(stripe1swap_array[:, -4:, :] == np.arange(3, -1, -1)[None, :, None])

    # Test positive fast and slow axes: values read forward instead of backward
    stripe_params = (208, 256, 4, 204, 0, 0, 1, 1, 204, 10, 2, 1)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)

    stripe2_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe2_array.shape == (10, 256, 208)
    # detector row numbers count forward in columns
    np.testing.assert_allclose(stripe2_array[0, 0], np.arange(0, 208))
    np.testing.assert_allclose(stripe2_array[-1, 0, 4:], np.arange(1840, 2044))
    # reference pixels are always 0, 1, 2, 3
    assert np.all(stripe2_array[:, :, :4] == np.arange(0, 4))


def test_generate_stripe_reference_superstripe_repeat_stripe_zero():
    """
    Test generate_stripe_reference in superstripe mode when repeat_stripe=0.

    When repeat_stripe=0 the read head does not reset between stripes;
    successive nreads2 blocks are separated by nskips2 rows advancing
    monotonically through the detector.
    """
    # Mock model for metadata
    model = datamodels.RampModel()
    sci_meta = model.meta
    sci_meta.exposure.ngroups = 5
    subarray_keys = _multistripe_subarray_keys()

    # Generate test array with pixel values equal to row number in detector frame.
    test_array = np.ones((2048, 256), dtype=int) * np.arange(2048)[:, None]

    # SUB204STRIPE_SOSS case but with no repeat, skip 8 instead of reading them
    # without 8 pixel skip: (208, 256, 4, 204, 0, 0, 1, 1, 204, 10, -2, -1)
    stripe_params = (200, 256, 4, 196, 0, 8, 0, 1, 204, 10, -2, -1)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)

    stripe_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe_array.shape == (10, 256, 200)
    np.testing.assert_allclose(stripe_array[0, 0], np.arange(199, -1, -1))
    np.testing.assert_allclose(stripe_array[-1, 0, :-4], np.arange(2035, 1839, -1))
    assert np.all(stripe_array[:, :, -4:] == np.arange(3, -1, -1))


def test_generate_stripe_reference_super_sub_stripe():
    """Test generate_stripe_reference in superstripe mode with substripes."""
    # Mock model for metadata
    model = datamodels.RampModel()
    sci_meta = model.meta
    sci_meta.exposure.ngroups = 5
    subarray_keys = _multistripe_subarray_keys()

    # Generate test array with pixel values equal to row number in detector frame.
    test_array = np.ones((2048, 256), dtype=int) * np.arange(2048)[:, None]

    # SUB204STRIPE_SOSS case but reduce nreads2 to 50, set skip2 to 1
    # to do substripes within superstripes.
    # without substripes: (208, 256, 4, 204, 0, 0, 1, 1, 204, 10, -2, -1)
    stripe_params = (216, 256, 4, 50, 0, 1, 1, 1, 204, 10, -2, -1)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)

    stripe_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe_array.shape == (10, 256, 216)

    # check the values in the first stripe:
    # reference pixels are interleaved between 4 50 pixel stripes
    start = 206
    expected = []
    for i in range(4):
        expected.append(np.arange(start, start - 50, -1))
        expected.append(np.arange(3, -1, -1))
        start -= 51
    expected = np.hstack(expected)
    np.testing.assert_allclose(stripe_array[0, 0], expected)

    # Same, but don't interleave: same stripe is read repeatedly
    stripe_params = (216, 256, 4, 50, 0, 1, 1, 0, 204, 10, -2, -1)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)

    stripe_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe_array.shape == (10, 256, 216)

    # check the values in the first stripe:
    # reference pixels are interleaved between 4 50 pixel stripes,
    # all reading the same detector pixels
    expected = []
    for i in range(4):
        expected.append(np.arange(53, -1, -1))
    expected = np.hstack(expected)
    np.testing.assert_allclose(stripe_array[0, 0], expected)


def test_generate_stripe_reference_super_sub_stripe_repeat_0():
    """Test generate_stripe_reference in superstripe mode with substripes, no repeat."""
    # Mock model for metadata
    model = datamodels.RampModel()
    sci_meta = model.meta
    sci_meta.exposure.ngroups = 5
    subarray_keys = _multistripe_subarray_keys()

    # Generate test array with pixel values equal to row number in detector frame.
    test_array = np.ones((2048, 256), dtype=int) * np.arange(2048)[:, None]

    # SUB204STRIPE_SOSS case but reduce nreads2 to 50, nreads1 to 1,
    # set skip1 to 2, skip2 to 4, set repeat to 0.
    # without substripes: (208, 256, 4, 204, 0, 0, 1, 1, 204, 10, -2, -1)
    stripe_params = (204, 256, 4, 50, 0, 1, 0, 0, 204, 10, -2, -1)
    _assign_metadata(sci_meta, subarray_keys, stripe_params)

    stripe_array = stripe_utils.generate_stripe_reference(
        science_detector_frame_transform(
            test_array, sci_meta.subarray.fastaxis, sci_meta.subarray.slowaxis
        ),
        model,
    )
    assert stripe_array.shape == (10, 256, 204)

    # check the values in the first stripe: should be 50 pixel reads,
    # followed by 4 reference pixels at the end
    start = 206
    expected = []
    for i in range(4):
        expected.append(np.arange(start, start - 50, -1))
        start -= 51
    expected.append(np.arange(3, -1, -1))
    expected = np.hstack(expected)
    np.testing.assert_allclose(stripe_array[0, 0], expected)

    # check the values in the last stripe: same, but start value is different
    expected[:-4] += 2042 - 206
    np.testing.assert_allclose(stripe_array[-1, 0], expected)


def test_generate_stripe_reference_invalid_shape(superstripe_model):
    ref_array = np.ones((2, 3, 4, 5, 6))
    with pytest.raises(ValueError, match=r"Unsupported shape: len\(ref_shape\) == 5"):
        stripe_utils.generate_stripe_reference(ref_array, superstripe_model)


def test_stripe_read_substripe(substripe_model):
    mask_model = make_superstripe_mask_model()
    striped_mask_model = stripe_utils.stripe_read(substripe_model, mask_model, ["dq"])
    assert type(striped_mask_model) is datamodels.MaskModel

    expected_shape = substripe_model.data.shape[-2:]
    assert striped_mask_model.dq.shape == expected_shape


def test_stripe_read_superstripe(superstripe_model):
    mask_model = make_superstripe_mask_model()
    striped_mask_model = stripe_utils.stripe_read(superstripe_model, mask_model, ["dq"])
    assert type(striped_mask_model) is datamodels.ReferenceFileModel

    expected_shape = (
        superstripe_model.meta.subarray.num_superstripe,
        *superstripe_model.data.shape[-2:],
    )
    assert striped_mask_model.dq.shape == expected_shape
