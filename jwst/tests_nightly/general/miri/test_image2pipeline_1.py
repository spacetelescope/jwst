import pytest
from jwst.pipeline.calwebb_image2 import Image2Pipeline

from jwst.tests.base_test import BaseJWSTTest


@pytest.mark.bigdata
class TestImage2Pipeline(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_image2pipeline', 'truth']

    def test_image2pipeline1(self, _jail):
        """
        Regression test of calwebb_image2 pipeline performed on MIRI data.
        """
        input_file = self.get_data('test_image2pipeline',
                                   'jw00001001001_01101_00001_mirimage_rate.fits')
        Image2Pipeline.call(input_file)

        outputs = [('jw00001001001_01101_00001_mirimage_cal.fits',
                    'jw00001001001_01101_00001_mirimage_cal_ref.fits',
                    ['primary','sci','err','dq','area']),
                    ('jw00001001001_01101_00001_mirimage_i2d.fits',
                    'jw00001001001_01101_00001_mirimage_i2d_ref.fits',
                    ['primary','sci','wht','con'])
                  ]
        self.compare_outputs(outputs)
