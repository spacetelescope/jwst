import pytest

from jwst.pipeline import Image2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs

from jwst.tests.base_classes import BaseJWSTTest


@pytest.mark.bigdata
class TestImage2Pipeline(BaseJWSTTest):
    input_loc = 'nircam'
    ref_loc = ['test_image2pipeline', 'truth']

    def test_image2pipeline2_cal(self):
        """
        Regression test of calwebb_image2 pipeline performed on NIRCam data.
        """
        input_file = self.get_data('test_image2pipeline',
                                   'jw82500001003_02101_00001_NRCALONG_rate.fits')
        output_file = 'jw82500001003_02101_00001_NRCALONG_cal.fits'

        collect_pipeline_cfgs('cfgs')
        Image2Pipeline.call(input_file,
                            config_file='cfgs/calwebb_image2.cfg',
                            output_file=output_file)

        outputs = [(output_file,
                    'jw82500001003_02101_00001_NRCALONG_cal_ref.fits',
                    ['primary','sci','err','dq','area']),
                   ('jw82500001003_02101_00001_NRCALONG_i2d.fits',
                    'jw82500001003_02101_00001_NRCALONG_i2d_ref.fits',
                    ['primary','sci','con','wht'])
                   ]
        self.compare_outputs(outputs)
