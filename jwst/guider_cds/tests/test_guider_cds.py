import numpy as np
import pytest

from stdatamodels.jwst import datamodels

from jwst.guider_cds.guider_cds import get_dataset_info, guider_cds


@pytest.fixture
def make_guider_image():
    """Generate science image"""

    image = datamodels.GuiderRawModel()

    image.meta.instrument.name = 'FGS'
    image.meta.exposure.frame_time = 234.3423235
    image.meta.exposure.ngroups = 4
    image.meta.exposure.group_time = 465.643643
    image.meta.exposure.type = "FGS_FINEGUIDE"

    image.data = np.random.rand(4, 10, 10, 10)

    return image


def test_get_dataset_info(make_guider_image):
    """Make sure information assigned to datamodel is retrieved correctly."""

    model = make_guider_image

    imshape, n_int, grp_time, exp_type = get_dataset_info(model)

    assert imshape == (10, 10)
    assert n_int == model.data.shape[0]
    assert grp_time == model.meta.exposure.group_time
    assert exp_type == model.meta.exposure.type


def test_guider_cds_fineguide_mode(make_guider_image):
    """Test the fine guiding mode"""

    model = make_guider_image

    truth = np.zeros(model.data.shape)

    result = guider_cds(model)

    n_int = model.data.shape[0]
    imshape = (model.data.shape[2], model.data.shape[3])

    slope_int_cube = np.zeros((n_int,) + imshape, dtype=np.float32)

    for num_int in range(0, n_int):
        data_sect = model.data[num_int, :, :, :]
        first_4 = data_sect[:4, :, :].mean(axis=0)
        last_4 = data_sect[-4:, :, :].mean(axis=0)
        slope_int_cube[num_int, :, :] = last_4 - first_4

        truth = slope_int_cube / model.meta.exposure.group_time

    assert np.allclose(result.data, truth)


@pytest.mark.parametrize("exptype", ['FGS_ACQ1', 'FGS_ACQ2', 'FGS_TRACK'])
def test_guider_cds_acq_track_modes(exptype, make_guider_image):
    """Test acq and track exptypes"""

    model = make_guider_image
    model.meta.exposure.type = exptype

    truth = np.zeros(model.data.shape)

    result = guider_cds(model)

    n_int = model.data.shape[0]
    imshape = (model.data.shape[2], model.data.shape[3])
    slope_int_cube = np.zeros((n_int,) + imshape, dtype=np.float32)

    for num_int in range(0, n_int):
        data_sect = model.data[num_int, :, :, :]
        grp_last = data_sect[1, :, :]
        grp_first = data_sect[0, :, :]
        slope_int_cube[num_int, :, :] = grp_last - grp_first

        truth = slope_int_cube / model.meta.exposure.group_time

    assert np.allclose(result.data, truth)


@pytest.mark.parametrize("exptype", ['FGS_ID-IMAGE', 'FGS_ID-STACK'])
def test_guider_cds_id_modes(exptype, make_guider_image):
    """Test fgs id exptypes"""

    model = make_guider_image
    model.meta.exposure.type = exptype

    result = guider_cds(model)

    n_int = model.data.shape[0]
    imshape = (model.data.shape[2], model.data.shape[3])

    truth = np.zeros((1,) + imshape)

    for num_int in range(0, n_int):
        data_sect = model.data[num_int, :, :, :]
        grp_last = data_sect[1, :, :]
        grp_first = data_sect[0, :, :]

        if num_int == 0:
            diff_int0 = grp_last - grp_first
        if num_int == 1:
            diff_int1 = grp_last - grp_first

    truth[0, :, :] = np.minimum(diff_int1, diff_int0) / model.meta.exposure.group_time

    assert np.allclose(result.data[0, :, :], truth[0, :, :])


def test_unit_assignment(make_guider_image):
    """Test that correct units are returned"""

    model = make_guider_image

    result = guider_cds(model)

    assert result.meta.bunit_data == 'DN/s'


def test_table_extensions(make_guider_image):
    """Test that tables are assigned to result of pipeline"""

    model = make_guider_image

    model.planned_star_table = np.arange(0, 11)
    model.flight_star_table = np.arange(0, 11)
    model.pointing_table = np.arange(0, 11)
    model.centroid_table = np.arange(0, 11)
    model.track_sub_table = np.arange(0, 11)

    result = guider_cds(model)

    assert 'planned_star_table' in result
    assert 'flight_star_table' in result
    assert 'pointing_table' in result
    assert 'centroid_table' in result
    assert 'track_sub_table' in result


def test_err_nonzero(make_guider_image):
    """Make sure that the ERR array in output are not all zero."""

    model = make_guider_image

    result = guider_cds(model)

    assert result.err.max() > 0
