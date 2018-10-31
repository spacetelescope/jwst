import pytest

from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

from jwst.tests.base_classes import BaseJWSTTest, raw_from_asn


@pytest.mark.bigdata
class TestSpec2Pipeline(BaseJWSTTest):
    input_loc = 'niriss'
    ref_loc = ['test_spec2pipeline', 'truth']
    test_dir = 'test_spec2pipeline'

    @pytest.mark.xfail(reason='https://github.com/STScI-JWST/jwst/issues/2007')
    def test_nis_wfss_spec2(self):
        """
        Regression test of calwebb_spec2 pipeline performed on NIRISS WFSS data.
        """
        asn_file = self.get_data(self.test_dir,
                                 'jw87600-a3001_20171109T145456_spec2_001_asn.json')
        for file in raw_from_asn(asn_file):
            self.get_data(self.test_dir, file)

        pipe = Spec2Pipeline()
        pipe.save_bsub = True
        pipe.save_results = True
        pipe.resample_spec.save_results = True
        pipe.extract_1d.save_results = True
        pipe.run(asn_file)

        outputs = [('jw87600017001_02101_00002_nis_cal.fits',
                    'jw87600017001_02101_00002_nis_cal_ref.fits'),
                   ('jw87600017001_02101_00002_nis_x1d.fits',
                    'jw87600017001_02101_00002_nis_x1d_ref.fits')]
        self.compare_outputs(outputs)
