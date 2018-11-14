import pytest
from jwst.pipeline import Image2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs

from jwst.tests.base_classes import BaseJWSTTest


@pytest.mark.bigdata
class TestImage2Pipeline(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_image2pipeline', 'truth']

    def test_image2pipeline1(self):
        """
        Regression test of calwebb_image2 pipeline performed on MIRI data.
        """
        input_file = self.get_data('test_image2pipeline',
                                   'jw00001001001_01101_00001_mirimage_rate.fits')
        collect_pipeline_cfgs('cfgs')
        Image2Pipeline.call(input_file,
                            config_file='cfgs/calwebb_image2.cfg'
                           )

        outputs = [('jw00001001001_01101_00001_mirimage_cal.fits',
                    'jw00001001001_01101_00001_mirimage_cal_ref.fits',
                    ['primary','sci','err','dq','area']),
                    ('jw00001001001_01101_00001_mirimage_i2d.fits',
                    'jw00001001001_01101_00001_mirimage_i2d_ref.fits',
                    ['primary','sci','wht','con'])
                  ]
        self.compare_outputs(outputs)
