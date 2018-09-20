import pytest
from jwst.pipeline.calwebb_image2 import Image2Pipeline

from ..resources import FGSTest

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


class TestImage2Pipeline(FGSTest):
    ref_loc = ['test_image2pipeline']#, 'truth']

    def test_fgs_image2pipeline1(self):
        """

        Regression test of calwebb_image2 pipeline performed on FGS imaging mode data.

        """
        input_file = self.get_data('test_image2pipeline',
                                   'jw86500007001_02101_00001_GUIDER2_rate.fits')
        output_file = 'jw86500007001_02101_00001_GUIDER2_cal.fits'

        Image2Pipeline.call(input_file,output_file = output_file)

        outputs = [(output_file, 'jw86500007001_02101_00001_GUIDER2_cal_ref.fits',
                    ['primary','sci','err','dq','area'])
                  ]
        self.compare_outputs(outputs)
