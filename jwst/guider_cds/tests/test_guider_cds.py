import numpy as np
import pytest

from jwst import datamodels
from jwst.guider_cds.guider_cds import get_dataset_info, guider_cds



def test_get_dataset_info():
    """Make sure information assined to datamodel is retrieved correctly."""

    model = datamodels.GuiderRawModel()

    model.meta.instrument.name = 'FGS'
    model.meta.exposure.frame_time = 234.3423235
    model.meta.exposure.ngroups = 4
    model.meta.exposure.group_time = 465.643643
    model.meta.exposure.type = 'FGS_FINEGUIDE'

    model.data = np.random.rand(1,10,10,10)

    imshape, n_int, grp_time, exp_type = get_dataset_info(model)

    assert (imshape == (10,10) and
            n_int == model.data.shape[0] and
            grp_time == model.meta.exposure.group_time and
            exp_type == model.meta.exposure.type)

@pytest.mark.parametrize("exptype", ['FGS_FINEGUIDE','FGS_ACQ1', 'FGS_ACQ2', 'FGS_TRACK', 'FGS_ID-IMAGE', 'FGS_ID-STACK'])
def test_guider_cds(exptype):
    
    model = datamodels.GuiderRawModel()

    model.meta.instrument.name = 'FGS'
    model.meta.exposure.frame_time = 234.3423235
    model.meta.exposure.ngroups = 4
    model.meta.exposure.group_time = 465.643643
    model.meta.exposure.type = exptype

    model.data = np.random.rand(4,10,10,10)

    result = guider_cds(model)

    n_int = model.data.shape[0]
    imshape = (model.data.shape[2], model.data.shape[3])
    slope_int_cube = np.zeros((n_int,) + imshape, dtype=np.float32)
    

    for num_int in range(0, n_int):
        data_sect = model.data[num_int, :, :, :]
        if model.meta.exposure.type ==  'FGS_FINEGUIDE':
            first_4 = data_sect[:4, :, :].mean(axis=0)
            last_4 = data_sect[-4:, :, :].mean(axis=0)
            slope_int_cube[num_int, :, :] = last_4 - first_4
        elif exptype[:6] == 'FGS_ID':
            grp_last = data_sect[1, :, :]
            grp_first = data_sect[0, :, :]

            if num_int == 0:
                diff_int0 = grp_last - grp_first
            if num_int == 1:
                diff_int1 = grp_last - grp_first
        else:
            grp_last = data_sect[1, :, :]
            grp_first = data_sect[0, :, :]
            slope_int_cube[num_int, :, :] = grp_last - grp_first
    
    if exptype[:6] == 'FGS_ID':
        truth[0, :, :] = np.minimum(diff_int1, diff_int0) / model.meta.exposure.group_time
    else:
        truth = slope_int_cube /  model.meta.exposure.group_time

    assert np.array_equal(result.data, truth)