import pytest
from jwst.pipeline.calwebb_image2 import Image2Pipeline

from ..resources import NIRCamTest

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

class TestImage2Pipeline(NIRCamTest):

    ref_loc = ['test_image2pipeline']

    def test_image2pipeline2_cal(self):
        """
        Regression test of calwebb_image2 pipeline performed on NIRCam data.
        """
        input_file = self.get_data('test_image2pipeline',
                                   'jw82500001003_02101_00001_NRCALONG_rate_ref.fits')
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
