import os
import stdatamodels.jwst.datamodels as dm
from jwst.extract_1d import Extract1dStep
from jwst.photom import PhotomStep

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
TEST_FILES = {
    'SOSS/FULL': os.path.join(data_dir, 'jw01209001001_03101_00001_nis_pathloss_dummy.fits'),
}


def test_expected_skip_niriss_soss_full():

    infile = TEST_FILES['SOSS/FULL']
    with dm.open(infile) as model:
        result = Extract1dStep().process(model)
        result2 = PhotomStep().process(result)
        assert result2.meta.cal_step.photom == 'SKIPPED'
        assert result2.meta.cal_step.extract_1d == 'SKIPPED'

    assert result2.data == model.data
