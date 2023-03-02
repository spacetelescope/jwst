import pytest

from stdatamodels.jwst import datamodels

from jwst.resample import ResampleSpecStep, ResampleStep


@pytest.mark.parametrize('resample_class', [ResampleSpecStep, ResampleStep])
def test_multi_integration_input(resample_class):
    cube = datamodels.CubeModel((5, 100, 100))
    cube.meta.instrument.name = 'MIRI'
    cube.meta.observation.date = '2018-09-07'
    cube.meta.observation.time = '10:32:20.181'

    # Resample can't handle cubes, so it should fail
    with pytest.raises(RuntimeError):
        resample_class().call(cube)
