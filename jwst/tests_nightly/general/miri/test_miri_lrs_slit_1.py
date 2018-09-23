import pytest
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

from jwst.tests.base_test import MIRITest


@pytest.mark.bigdata
class TestSpec2Pipeline(MIRITest):
    ref_loc = ['test_lrs_slit', 'truth']

    def test_miri_lrs_slit_1(self):
        """
        Regression test of calwebb_spec2 pipeline performed on a single
        MIRI LRS fixed-slit exposure.
        """
        input_file = self.get_data('test_lrs_slit',
                              'jw00035001001_01101_00001_MIRIMAGE_rate.fits')

        step = Spec2Pipeline()
        step.save_bsub=True,
        step.save_results=True
        step.resample_spec.save_results = True
        step.cube_build.save_results = True
        step.extract_1d.save_results = True
        step.run(input_file)

        outputs = [{'files':('jw00035001001_01101_00001_MIRIMAGE_cal.fits',
                            'jw00035001001_01101_00001_MIRIMAGE_cal_ref.fits',
                            ['primary','sci','err','dq','relsens']),
                    'pars': {}
                   },
                   {'files':('jw00035001001_01101_00001_MIRIMAGE_x1d.fits',
                             'jw00035001001_01101_00001_MIRIMAGE_x1d_ref.fits',
                             ['primary','extract1d']),
                    'pars': {}
                    }
                  ]
        self.compare_outputs(outputs)
