import logging

from jwst.cube_build.ifu_cube import _check_cube_size


def test_typical_cube_does_not_warn(caplog):
    with caplog.at_level(logging.WARNING, logger="jwst.cube_build.ifu_cube"):
        total = _check_cube_size(60, 60, 8000)
    assert total == 28_800_000
    assert "exceeding" not in caplog.text


def test_large_cube_warns_before_allocation(caplog):
    with caplog.at_level(logging.WARNING, logger="jwst.cube_build.ifu_cube"):
        total = _check_cube_size(1001, 1001, 100)
    assert total == 100_200_100
    assert "100200100 voxels" in caplog.text
    assert "100000000-voxel warning threshold" in caplog.text
