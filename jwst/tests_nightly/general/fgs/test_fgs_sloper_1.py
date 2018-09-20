import pytest
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

from ..resources import FGSTest

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


class TestSloperPipeline(FGSTest):
    ref_loc = ['test_sloperpipeline']#, 'truth']

    def test_fgs_detector1_1(self):
        """
        Regression test of calwebb_detector1 pipeline performed on FGS imaging mode data.
        """
        input_file = self.get_data('test_sloperpipeline',
                                    'jw86500007001_02101_00001_GUIDER2_uncal.fits')
        pipe = Detector1Pipeline()
        pipe.ipc.skip = True
        pipe.refpix.odd_even_columns = True
        pipe.refpix.use_side_ref_pixels = True
        pipe.refpix.side_smoothing_length = 11
        pipe.refpix.side_gain = 1.0
        pipe.refpix.odd_even_rows = True
        pipe.jump.rejection_threshold = 250.0
        pipe.persistence.skip = True
        pipe.ramp_fit.save_opt = False
        pipe.save_calibrated_ramp = True
        pipe.output_file = 'jw86500007001_02101_00001_GUIDER2_rate.fits'

        pipe.run(input_file)

        outputs = [('jw86500007001_02101_00001_GUIDER2_ramp.fits',
                    'jw86500007001_02101_00001_GUIDER2_ramp_ref.fits',
                    ['primary','sci','err','groupdq','pixeldq']),
                   ('jw86500007001_02101_00001_GUIDER2_rateints.fits',
                    'jw86500007001_02101_00001_GUIDER2_rateints_ref.fits',
                    ['primary','sci','err','dq']),
                   ('jw86500007001_02101_00001_GUIDER2_rate.fits',
                    'jw86500007001_02101_00001_GUIDER2_rate_ref.fits',
                    ['primary','sci','err','dq'])
                  ]
        self.compare_outputs(outputs)
