import pytest

from jwst.assign_wcs.tests.helpers import make_mock_dhs_nrca1_rate
from jwst.dq_init.tests.helpers import make_superstripe_model
from jwst.lib import stripe_utils


@pytest.fixture(scope="module")
def substripe_model():
    return make_mock_dhs_nrca1_rate()


@pytest.fixture(scope="module")
def superstripe_model():
    return make_superstripe_model()


def test_generate_substripe_ranges(substripe_model):
    all_ranges = stripe_utils.generate_substripe_ranges(substripe_model)
    expected = {
        "full": {0: [1527, 1567], 1: [1652, 1692], 2: [1777, 1817], 3: [1902, 1942]},
        "subarray": {0: [1, 41], 1: [42, 82], 2: [83, 123], 3: [124, 164]},
        "reference_full": {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]},
        "reference_subarray": {0: [0, 1], 1: [41, 42], 2: [82, 83], 3: [123, 124]},
    }
    assert all_ranges == expected


def test_generate_substripe_ranges_invalid_input(substripe_model):
    model = substripe_model.copy()
    model.meta.subarray.multistripe_reads2 = 0

    with pytest.raises(ValueError, match="Invalid value for multistripe_reads2"):
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


def test_generate_superstripe_ranges_invalid_input(superstripe_model):
    model = superstripe_model.copy()
    model.meta.subarray.multistripe_reads2 = 0

    with pytest.raises(ValueError, match="Invalid value for multistripe_reads2"):
        stripe_utils.generate_superstripe_ranges(model)
