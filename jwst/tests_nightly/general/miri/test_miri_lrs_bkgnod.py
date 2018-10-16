import pytest
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

from jwst.tests.base_test import BaseJWSTTest


@pytest.mark.bigdata
class TestSpec2Pipeline(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_spec2pipeline', 'truth']

    test_dir = 'test_spec2pipeline'

    def test_miri_lrs_bkgnod(self):
        """

        Regression test of calwebb_spec2 pipeline performed on an association
        of nodded MIRI LRS fixed-slit exposures.

        """
        asn_file = self.get_data(self.test_dir,
                                   'lrs_bkgnod_asn.json')
        for file in self.raw_from_asn(asn_file):
            self.get_data(self.test_dir, file)

        step = Spec2Pipeline()
        step.save_bsub=True,
        step.save_results=True
        step.resample_spec.save_results = True
        step.cube_build.save_results = True
        step.extract_1d.save_results = True
        step.run(asn_file)

        outputs = [('test_lrs1_bsub.fits', 'test_lrs1_bsub_ref.fits',
                    ['primary','sci','err','dq']),
                   ('test_lrs2_bsub.fits','test_lrs2_bsub_ref.fits',
                    ['primary','sci','err','dq']),
                   ('test_lrs3_bsub.fits','test_lrs3_bsub_ref.fits',
                    ['primary','sci','err','dq']),
                   ('test_lrs4_bsub.fits','test_lrs4_bsub_ref.fits',
                    ['primary','sci','err','dq']),
                   ('test_lrs1_cal.fits', 'test_lrs1_cal_ref.fits',
                    ['primary','sci','err','dq','relsens']),
                   ('test_lrs2_cal.fits', 'test_lrs2_cal_ref.fits',
                    ['primary','sci','err','dq','relsens']),
                   ('test_lrs3_cal.fits', 'test_lrs3_cal_ref.fits',
                    ['primary','sci','err','dq','relsens']),
                   ('test_lrs4_cal.fits', 'test_lrs4_cal_ref.fits',
                    ['primary','sci','err','dq','relsens'])
        ]
        self.compare_outputs(outputs)

    def test_miri_lrs_slit_1(self):
        """

        Regression test of calwebb_spec2 pipeline performed on a single
        MIRI LRS fixed-slit exposure.

        """
        input_file = self.get_data(self.test_dir,
                                   'jw00035001001_01101_00001_MIRIMAGE_rate.fits')

        step = Spec2Pipeline()
        step.save_bsub=True,
        step.save_results=True
        step.resample_spec.save_results = True
        step.cube_build.save_results = True
        step.extract_1d.save_results = True
        step.run(input_file)

        outputs = [('jw00035001001_01101_00001_MIRIMAGE_cal.fits',
                    'jw00035001001_01101_00001_MIRIMAGE_cal_ref.fits',
                    ['primary','sci','err','dq','relsens']),
                   ('jw00035001001_01101_00001_MIRIMAGE_x1d.fits',
                    'jw00035001001_01101_00001_MIRIMAGE_x1d_ref.fits',
                    ['primary','extract1d'])
                   ]
        self.compare_outputs(outputs)

    def test_miri_lrs_slit_1b(self):
        """
        Regression test of calwebb_spec2 pipeline performed on a single
        MIRI LRS fixed-slit exposure with multiple integrations.  Compare _calints.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw00035001001_01101_00001_MIRIMAGE_rateints.fits')

        step = Spec2Pipeline()
        step.save_bsub=True,
        step.save_results=True
        step.extract_1d.save_results = True
        step.run(input_file)

        outputs = [('jw00035001001_01101_00001_MIRIMAGE_calints.fits',
                    'jw00035001001_01101_00001_MIRIMAGE_calints_ref.fits',
                    ['primary','sci','err','dq','relsens']),
                    ('jw00035001001_01101_00001_MIRIMAGE_x1dints.fits',
                     'jw00035001001_01101_00001_MIRIMAGE_x1dints_ref.fits',
                     ['primary', ('extract1d', 1), ('extract1d', 2), ('extract1d', 3), ('extract1d', 4)]
                    )
                   ]
        self.compare_outputs(outputs)

    def test_mrs2pipeline1(self):
        """

        Regression test of calwebb_spec2 pipeline performed on MIRI MRS data.

        """
        test_dir = 'test_mrs2pipeline'
        self.ref_loc = ['test_mrs2pipeline', 'truth']

        input_file = self.get_data(test_dir,
                                   'jw80500018001_02101_00002_MIRIFUSHORT_rate.fits')
        step = Spec2Pipeline()
        step.save_bsub=True,
        step.save_results=True
        step.resample_spec.save_results = True
        step.cube_build.save_results = True
        step.extract_1d.save_results = True
        step.run(input_file)

        outputs = [('jw80500018001_02101_00002_MIRIFUSHORT_cal.fits',
                    'jw80500018001_02101_00002_MIRIFUSHORT_cal_ref.fits',
                    ['primary','sci','err','dq']),
                   ('jw80500018001_02101_00002_MIRIFUSHORT_s3d.fits',
                    'jw80500018001_02101_00002_MIRIFUSHORT_s3d_ref.fits',
                    ['primary','sci','err','dq','wmap']),
                   ('jw80500018001_02101_00002_MIRIFUSHORT_x1d.fits',
                    'jw80500018001_02101_00002_MIRIFUSHORT_x1d_ref.fits',
                    ['primary','extract1d'])
                   ]
        self.compare_outputs(outputs)

    def test_mrs_spec2(self):
        """

        Regression test of calwebb_spec2 pipeline performed on MIRI MRS data.

        """
        self.rtol = 0.000001
        input_file = self.get_data(self.test_dir,
                                   'jw10001001001_01101_00001_mirifushort_rate.fits')
        step = Spec2Pipeline()
        step.save_bsub=True,
        step.save_results=True
        step.resample_spec.save_results = True
        step.cube_build.save_results = True
        step.extract_1d.save_results = True
        step.run(input_file)

        outputs = [('jw10001001001_01101_00001_mirifushort_cal.fits',
                    'jw10001001001_01101_00001_mirifushort_cal_ref.fits',
                    ['primary','sci','err','dq','relsens2d']),
                    ('jw10001001001_01101_00001_mirifushort_s3d.fits',
                     'jw10001001001_01101_00001_mirifushort_s3d_ref.fits',
                     ['primary','sci','err','dq','wmap']),
                    ('jw10001001001_01101_00001_mirifushort_x1d.fits',
                     'jw10001001001_01101_00001_mirifushort_x1d_ref.fits',
                     ['primary','extract1d'])
                  ]
        self.compare_outputs(outputs)
