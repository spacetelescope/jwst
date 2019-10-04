from jwst.datamodels import SlitModel

from jwst.resample.resample_spec import find_dispersion_axis


def test_find_dispersion_axis():
    """
    Test the find_dispersion_axis() function
    """
    dm = SlitModel()

    dm.meta.wcsinfo.dispersion_direction = 1    # horizontal
    assert find_dispersion_axis(dm) == 0        # X axis for wcs functions

    dm.meta.wcsinfo.dispersion_direction = 2    # vertical
    assert find_dispersion_axis(dm) == 1        # Y axis for wcs functions
