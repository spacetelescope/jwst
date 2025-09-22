import numpy as np
import pytest
from astropy.modeling.mappings import Identity
from stcal.alignment import sregion_to_footprint

from jwst.resample.combine_sregions import (
    _combine_footprints,
    _polygons_to_sregion,
    combine_sregions,
)


@pytest.fixture
def simple_footprint_pair():
    """A pair of overlapping footprints."""
    footprint1 = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]
    )
    footprint2 = np.array(
        [
            [0.5, 0.5],
            [1.5, 0.5],
            [1.5, 1.5],
            [0.5, 1.5],
        ]
    )
    return np.array([footprint1, footprint2])


def _rotate_to_start(arr, start):
    """
    Rotate a closed polygon array so that it starts at the given point.

    This is useful for testing because we don't care about the starting point of the polygon,
    but we do care about the order of the points.
    """
    idx = np.where(np.all(np.isclose(arr, start, rtol=1e-12), axis=1))[0]
    if len(idx) == 0:
        raise ValueError("Start point not found in array")
    idx = idx[0]
    return np.roll(arr, -idx, axis=0)


def test_combine_footprints_simple(simple_footprint_pair):
    combined = _combine_footprints(simple_footprint_pair)

    # this is always a list but in this case there is only one combined region
    assert isinstance(combined, list)
    assert len(combined) == 1
    combined = combined[0]
    assert combined.shape == (8, 2)  # no duplicate last point

    expected = np.array(
        [
            [0.0, 0.0],
            [0.0, 1.0],
            [0.5, 1.0],
            [0.5, 1.5],
            [1.5, 1.5],
            [1.5, 0.5],
            [1.0, 0.5],
            [1.0, 0.0],
        ]
    )
    combined = _rotate_to_start(combined, expected[0])

    np.testing.assert_allclose(combined, expected)


@pytest.fixture
def complex_footprint_set():
    """
    Set of nine footprints in three groups.

    * Five footprints of the same size that overlap and combine into one region.
    * One standalone footprint that does not overlap with any others.
    * Two footprints that are rotated by 45 degrees w.r.t. one another and combine into one region,
      plus a smaller footprint inside the combined region that does not affect the union.
    """
    offsets = np.array(
        [
            [-2.0, 0.0],
            [0.0, 0.0],
            [0.1, 0.1],
            [0.5, 0.2],
            [0.9, 0.4],
            [1.0, 1.0],
            [-2.0, 2.0],
            [-2.0, 2.0],
            [-2.0, 2.0],
        ]
    )
    offsets[:, 0] += 180

    square = np.array(
        [
            [0, 0],
            [1, 0],
            [1, 1],
            [0, 1],
        ]
    )
    rotated_square = np.array(
        [
            [0.5, -0.25],
            [1.25, 0.5],
            [0.5, 1.25],
            [-0.25, 0.5],
        ]
    )
    mini_square = np.array(
        [
            [0.3, 0.3],
            [0.7, 0.3],
            [0.7, 0.7],
            [0.3, 0.7],
        ]
    )
    footprints = np.empty((offsets.shape[0], square.shape[0], 2))
    for i, offset in enumerate(offsets):
        if i == offsets.shape[0] - 2:
            footprint = rotated_square + offset
        elif i == offsets.shape[0] - 1:
            footprint = mini_square + offset
        else:
            footprint = square + offset
        footprints[i] = footprint
    return footprints


def test_combine_footprints_complex(complex_footprint_set):
    combined = _combine_footprints(complex_footprint_set)

    # this is always a list but in this case there are three combined regions
    assert isinstance(combined, list)
    assert len(combined) == 3

    region_shapes = sorted([region.shape for region in combined])
    expected_shapes = [(4, 2), (16, 2), (20, 2)]
    assert region_shapes == expected_shapes


@pytest.fixture
def shared_edge_footprints():
    """Two footprints that share an edge."""
    footprint1 = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]
    )
    footprint2 = np.array(
        [
            [1.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
        ]
    )
    return np.array([footprint1, footprint2])


def test_combine_footprints_shared_edge(shared_edge_footprints):
    """
    Test combining footprints that share an edge.

    The expected result is a single region without duplicate points along the shared edge.
    """
    combined = _combine_footprints(shared_edge_footprints)

    # this is always a list but in this case there is only one combined region
    assert isinstance(combined, list)
    assert len(combined) == 1
    combined = combined[0]
    assert combined.shape == (4, 2)  # no duplicate last point

    expected = np.array(
        [
            [0.0, 0.0],
            [0.0, 1.0],
            [2.0, 1.0],
            [2.0, 0.0],
        ]
    )
    combined = _rotate_to_start(combined, expected[0])

    np.testing.assert_allclose(combined, expected)


@pytest.fixture
def shared_vertex_footprints():
    """Two footprints that share a vertex."""
    footprint1 = np.array(
        [
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]
    )
    footprint2 = np.array(
        [
            [1.0, 0.0],
            [2.0, 0.0],
            [2.0, -1.0],
            [1.0, -1.0],
        ]
    )
    return np.array([footprint1, footprint2])


def test_combine_footprints_shared_vertex(shared_vertex_footprints):
    """
    Test combining footprints that share a vertex.

    The expected result is two separate regions.
    """
    combined = _combine_footprints(shared_vertex_footprints)

    assert isinstance(combined, list)
    assert len(combined) == 2

    expected = [shared_vertex_footprints[0], shared_vertex_footprints[1]]

    # We want the two regions to be exactly as expected, but we do not care about the order
    matched = [False, False]
    for region in combined:
        for i, exp in enumerate(expected):
            try:
                region_rot = _rotate_to_start(region, exp[0])
                np.testing.assert_allclose(region_rot, exp)
                matched[i] = True
                break
            except (AssertionError, ValueError):
                continue
    assert all(matched)


@pytest.fixture
def hole_footprint():
    """Four footprints that, when combined, leave a hole in the middle."""
    left = np.array(
        [
            [0.0, 0.0],
            [0.0, 5.0],
            [1.0, 5.0],
            [1.0, 0.0],
        ]
    )
    top = np.array(
        [
            [0.0, 4.0],
            [0.0, 5.0],
            [5.0, 5.0],
            [5.0, 4.0],
        ]
    )
    right = np.array(
        [
            [4.0, 0.0],
            [4.0, 5.0],
            [5.0, 5.0],
            [5.0, 0.0],
        ]
    )
    bottom = np.array(
        [
            [0.0, 0.0],
            [5.0, 0.0],
            [5.0, 1.0],
            [0.0, 1.0],
        ]
    )
    return np.array([left, top, right, bottom])


def test_combine_footprints_hole(hole_footprint):
    """
    Test combining footprints where their combination has a hole in the middle.

    The expected result is a single region with NO hole, because there is no easy way to specify
    an S_REGION keyword with a hole in the middle,
    and MAST cannot handle that case.
    """
    combined = _combine_footprints(hole_footprint)

    # this is always a list but in this case there is only one combined region
    assert isinstance(combined, list)
    assert len(combined) == 1
    combined = combined[0]
    assert combined.shape == (4, 2)  # no duplicate last point

    expected = np.array(
        [
            [0.0, 0.0],
            [0.0, 5.0],
            [5.0, 5.0],
            [5.0, 0.0],
        ]
    )
    combined = _rotate_to_start(combined, expected[0])

    np.testing.assert_allclose(combined, expected)


@pytest.fixture
def det2world():
    """Simple WCS transform that is just an identity transform."""
    return Identity(2)


def _footprints_to_sregion_list(footprints):
    """Helper function to convert a list of footprints to a list of S_REGION strings."""
    return [_polygons_to_sregion([footprint]) for footprint in footprints]


def test_sregion_pair(simple_footprint_pair, det2world):
    """Test for a pair of overlapping footprints as S_REGION strings."""
    sregion_list = _footprints_to_sregion_list(simple_footprint_pair)
    # since we're constructing this with code, let's first make sure the list looks right
    assert (
        sregion_list[0]
        == "POLYGON ICRS  0.000000000 0.000000000 1.000000000 0.000000000 1.000000000 1.000000000 0.000000000 1.000000000"
    )
    assert (
        sregion_list[1]
        == "POLYGON ICRS  0.500000000 0.500000000 1.500000000 0.500000000 1.500000000 1.500000000 0.500000000 1.500000000"
    )

    # now combine them and test their combination
    combined_sregion = combine_sregions(sregion_list, det2world)
    expected = "POLYGON ICRS  1.000000000 0.000000000 0.000000000 0.000000000 0.000000000 1.000000000 0.500000000 1.000000000 0.500000000 1.500000000 1.500000000 1.500000000 1.500000000 0.500000000 1.000000000 0.500000000"
    assert combined_sregion == expected


def test_sregion_complex(complex_footprint_set, det2world):
    """Test for a complex set of footprints as S_REGION strings."""
    sregion_list = _footprints_to_sregion_list(complex_footprint_set)
    combined_sregion = combine_sregions(sregion_list, det2world)

    assert combined_sregion.count("POLYGON ICRS") == 3
    assert combined_sregion.startswith("POLYGON ICRS")
    sregions_out = [s.strip() for s in combined_sregion.split("POLYGON ICRS")][1:]
    footprints_out = [sregion_to_footprint(sregion) for sregion in sregions_out]

    # Too tedious to check exact coordinates, so just check the shapes of the regions.
    region_shapes = sorted([region.shape for region in footprints_out])
    expected_shapes = [(4, 2), (16, 2), (20, 2)]
    assert region_shapes == expected_shapes


def test_sregion_intersection(complex_footprint_set, det2world):
    """
    Test for a complex set of footprints intersected with a bounding box.

    The bounding box is chosen to intersect only two of the three combined regions,
    and to cut both of them.
    """
    bbox = np.array(
        [
            [178.5, 0.5],
            [181.5, 0.5],
            [181.5, 1.5],
            [178.5, 1.5],
        ]
    )
    sregion_list = _footprints_to_sregion_list(complex_footprint_set)
    combined_sregion = combine_sregions(sregion_list, det2world, intersect_footprint=bbox)

    assert combined_sregion.count("POLYGON ICRS") == 2
    assert combined_sregion.startswith("POLYGON ICRS")
    expected = "POLYGON ICRS  179.000000000 1.000000000 179.000000000 0.500000000 178.500000000 0.500000000 178.500000000 1.000000000  POLYGON ICRS  180.500000000 1.100000000 180.500000000 1.200000000 180.900000000 1.200000000 180.900000000 1.400000000 181.000000000 1.400000000 181.000000000 1.500000000 181.500000000 1.500000000 181.500000000 0.500000000 180.000000000 0.500000000 180.000000000 1.000000000 180.100000000 1.000000000 180.100000000 1.100000000"
    assert combined_sregion == expected


def test_sregion_no_overlap(complex_footprint_set, det2world):
    """Test error raise when intersection footprint does not overlap any input footprints."""
    bbox = np.array(
        [
            [10.0, 10.0],
            [11.0, 10.0],
            [11.0, 11.0],
            [10.0, 11.0],
        ]
    )
    sregion_list = _footprints_to_sregion_list(complex_footprint_set)
    with pytest.raises(ValueError) as excinfo:
        combine_sregions(sregion_list, det2world, intersect_footprint=bbox)
    assert str(excinfo.value) == "No overlap between input s_regions and intersection footprint"
