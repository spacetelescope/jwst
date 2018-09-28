import pytest
from jwst.pipeline.calwebb_image2 import Image2Pipeline

from jwst.tests.base_test import BaseJWSTTest


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
        Image2Pipeline.call(input_file,
                            output_file=output_file)

        outputs = [(output_file,
                    'jw82500001003_02101_00001_NRCALONG_cal_ref.fits',
                    ['primary','sci','err','dq','area']),
                   ('jw82500001003_02101_00001_NRCALONG_i2d.fits',
                    'jw82500001003_02101_00001_NRCALONG_i2d_ref.fits',
                    ['primary','sci','con','wht'])
                   ]
        self.compare_outputs(outputs)
