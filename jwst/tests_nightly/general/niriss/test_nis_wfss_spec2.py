import pytest

from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.tests.base_classes import BaseJWSTTest, raw_from_asn


@pytest.mark.bigdata
class TestSpec2Pipeline(BaseJWSTTest):
    input_loc = 'niriss'
    ref_loc = ['test_spec2pipeline', 'truth']
    test_dir = 'test_spec2pipeline'

    def test_nis_wfss_spec2(self):
        """
        Regression test of calwebb_spec2 pipeline performed on NIRISS WFSS data.
        """
        # Collect data
        asn_file = self.get_data(self.test_dir,
                                 'jw87600-a3001_20171109T145456_spec2_001_asn.json')
        for file in raw_from_asn(asn_file):
            self.get_data(self.test_dir, file)

        # Run the step
        collect_pipeline_cfgs('cfgs')
        Spec2Pipeline.call(asn_file, config_file='cfgs/calwebb_spec2.cfg', save_bsub=True)

        # Test results.
        outputs = [('jw87600017001_02101_00002_nis_cal.fits',
                    'jw87600017001_02101_00002_nis_cal_ref.fits'),
                   ('jw87600017001_02101_00002_nis_x1d.fits',
                    'jw87600017001_02101_00002_nis_x1d_ref.fits')]
        self.compare_outputs(outputs)
