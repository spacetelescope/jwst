import pytest
from jwst.pipeline.calwebb_image2 import Image2Pipeline

from jwst.tests.base_test import NIRCamTest


@pytest.mark.bigdata
class TestImage2Pipeline(NIRCamTest):

    ref_loc = ['test_image2pipeline']

    def test_image2pipeline2b(self):
        """
        Regression test of calwebb_image2 pipeline performed on NIRCam data,
        using a multiple integration rate (rateints) file as input.
        """
        input_file = self.get_data('test_image2pipeline',
                                   'jw82500001003_02101_00001_NRCALONG_rateints.fits')
        output_file = 'jw82500001003_02101_00001_NRCALONG_calints.fits'

        Image2Pipeline.call(input_file,
                            output_file=output_file)

        outputs = [(output_file,
                    'jw82500001003_02101_00001_NRCALONG_calints_ref.fits',
                    ['primary','sci','err','dq','area']
                    )
                  ]
        self.compare_outputs(outputs)
