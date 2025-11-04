import pytest

from jwst.lib.multiprocessing import determine_ncores


@pytest.mark.parametrize(
    "max_cores, num_cores, expected",
    [
        ("none", 4, 1),
        ("quarter", 4, 1),
        ("half", 4, 2),
        ("all", 4, 4),
        ("none", 1, 1),
        (None, 1, 1),
        (3, 5, 3),
        (100, 5, 5),
    ],
)
def test_determine_ncores(max_cores, num_cores, expected):
    """Test determine_ncores with various inputs."""
    assert determine_ncores(max_cores, num_cores) == expected


def test_determine_ncores_invalid():
    """Test that invalid max_cores raises ValueError."""
    with pytest.raises(ValueError, match="Invalid value for max_cores"):
        determine_ncores("invalid", 4)
