"""
Unit tests for master background subtraction
"""
import pytest

from jwst.master_background import MasterBackgroundStep
from jwst import datamodels
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs


@pytest.fixture(scope='module')
def user_background(tmpdir_factory):
    """Generate a user background spectrum"""

    filename = tmpdir_factory.mktemp('master_background_user_input')
    filename = str(filename.join('user_background.fits'))
    data = datamodels.SlitModel()
    data.save(filename)

    return filename


def _generate_data():
    """Generate data of each type of input for master background step"""

    image = datamodels.ImageModel((10, 10))
    multislit = datamodels.MultiSlitModel()
    multislit.slits.append(datamodels.SlitModel((10, 10)))
    ifu_image = datamodels.IFUImageModel((10, 10))
    container = datamodels.ModelContainer([image])
    cube = datamodels.CubeModel((2, 10, 10))

    # Return the data and a status dependent on whether the step can process it
    return [(image, 'COMPLETE'),
            (multislit, 'COMPLETE'),
            (ifu_image, 'COMPLETE'),
            (container, 'COMPLETE'),
            (cube, 'SKIPPED'),
            ]

@pytest.mark.parametrize('input_data, status', _generate_data())
def test_master_background_init(input_data, status, _jail, user_background):
    """Verify data can run through the step"""

    result = MasterBackgroundStep.call(input_data)
    result = MasterBackgroundStep.call(input_data, user_background=user_background)
    result = MasterBackgroundStep.call(input_data, save_background=True)
    collect_pipeline_cfgs('./config')
    result = MasterBackgroundStep.call(
        input_data,
        config_file='config/master_background.cfg'
        )

    assert type(input_data) is type(result)

    try:
        assert result.meta.cal_step.master_back_sub == status
    except AttributeError:
        for model in result:
            assert model.meta.cal_step.master_back_sub == status
