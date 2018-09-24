import pytest

from jwst.pipeline.calwebb_guider import GuiderPipeline

from jwst.tests.base_test import BaseJWSTTest


@pytest.mark.bigdata
class TestGuiderPipeline(BaseJWSTTest):
    input_loc = 'fgs'
    ref_loc = ['test_guiderpipeline', 'truth']
    test_dir = 'test_guiderpipeline'

    rtol = 0.000001

    def test_guider_pipeline1(self, _jail):
        """
        Regression test of calwebb_guider pipeline performed on ID-image data.
        """
        input_file = self.get_data(self.test_dir,
                                    'jw88600073001_gs-id_7_image-uncal.fits')

        GuiderPipeline.call(input_file,
                            output_file='jw88600073001_gs-id_7_image-cal.fits')

        # Compare calibrated ramp product
        outputs = [('jw88600073001_gs-id_7_image-cal.fits',
                    'jw88600073001_gs-id_7_image-cal_ref.fits',
                    ['primary','sci','dq'])
                  ]
        self.compare_outputs(outputs)


    def test_guider_pipeline2(self, _jail):
        """
        Regression test of calwebb_guider pipeline performed on ACQ-1 data.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw88600073001_gs-acq1_2016022183837_uncal.fits')

        GuiderPipeline.call(input_file,
                            output_file='jw88600073001_gs-acq1_2016022183837_cal.fits')

        # Compare calibrated ramp product
        outputs = [('jw88600073001_gs-acq1_2016022183837_cal.fits',
                    'jw88600073001_gs-acq1_2016022183837_cal_ref.fits',
                    ['primary','sci','dq'])

                  ]
        self.compare_outputs(outputs)


    def test_guider_pipeline3(self, _jail):
        """

        Regression test of calwebb_guider pipeline performed on ID STACKED data.

        """
        input_file = self.get_data(self.test_dir,
                                    'jw86600004001_gs-id_1_stacked-uncal.fits')

        GuiderPipeline.call(input_file,
                            output_file='jw86600004001_gs-id_1_stacked-cal.fits')

        # Compare calibrated ramp product
        outputs = [('jw86600004001_gs-id_1_stacked-cal.fits',
                    'jw86600004001_gs-id_1_stacked-cal_ref.fits',
                    ['primary','sci','dq'])

                   ]
        self.compare_outputs(outputs)
