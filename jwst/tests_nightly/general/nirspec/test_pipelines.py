import pytest

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step
from jwst.pipeline.calwebb_dark import DarkPipeline
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

from jwst.tests.base_test import BaseJWSTTest


@pytest.mark.bigdata
class TestNIRSpecPipelines(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_pipelines', 'truth']
    test_dir = 'test_pipelines'

    def test_nirspec_dark_pipeline(self, _jail):
        """
        Regression test of calwebb_dark pipeline performed on NIRSpec raw data.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw84500013001_02103_00003_NRS1_uncal.fits')

        pipe = DarkPipeline()
        pipe.suffix = 'dark'
        pipe.ipc.skip = True
        pipe.refpix.odd_even_columns = True
        pipe.refpix.use_side_ref_pixels = True
        pipe.refpix.side_smoothing_length = 11
        pipe.refpix.side_gain = 1.0
        pipe.refpix.odd_even_rows = True
        pipe.output_file = 'jw84500013001_02103_00003_NRS1_uncal.fits'

        pipe.run(input_file)

        outputs = [('jw84500013001_02103_00003_NRS1_dark.fits',
                    'jw84500013001_02103_00003_NRS1_dark_ref.fits',
                    ['primary','sci','err','pixeldq','groupdq'])]
        self.compare_outputs(outputs)

    def test_nrs_fs_brightobj_spec2(self, _jail):
        """
        Regression test of calwebb_spec2 pipeline performed on NIRSpec
        fixed-slit data that uses the NRS_BRIGHTOBJ mode (S1600A1 slit).
        """
        input_file = self.get_data(self.test_dir,
                                   'jw84600042001_02101_00001_nrs2_rateints.fits')
        collect_pipeline_cfgs()
        args = [
            'calwebb_tso_spec2.cfg',
            input_file
        ]
        Step.from_cmdline(args)

        outputs = [('jw84600042001_02101_00001_nrs2_calints.fits',
                    'jw84600042001_02101_00001_nrs2_calints_ref.fits'),
                   ('jw84600042001_02101_00001_nrs2_x1dints.fits',
                    'jw84600042001_02101_00001_nrs2_x1dints_ref.fits')
                  ]
        self.compare_outputs(outputs)

    def test_nrs_msa_spec2(self, _jail):
        """
        Regression test of calwebb_spec2 pipeline performed on NIRSpec MSA data.
        """
        input = 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod.fits'
        input_file = self.get_data(self.test_dir, input)
        self.get_data(self.test_dir, 'jw95065006001_0_short_msa.fits')

        # define step for use in test
        step = Spec2Pipeline()
        step.save_bsub = False
        step.output_use_model = True
        step.resample_spec.save_results = True
        step.extract_1d.save_results = True
        step.extract_1d.smoothing_length = 0
        step.extract_1d.bkg_order = 0
        step.run(input_file)

        outputs = [('F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod_cal.fits',
                    'f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod_cal_ref.fits'),
                   ('F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod_s2d.fits',
                    'f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod_s2d_ref.fits'),
                   ('F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod_x1d.fits',
                    'f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod_x1d_ref.fits')
                  ]
        self.compare_outputs(outputs)

    def test_nrs_msa_spec2b(self, _jail):
        """
        Regression test of calwebb_spec2 pipeline performed on NIRSpec MSA data,
        including barshadow correction.
        """
        input = 'jw95065_nrs_msaspec_barshadow.fits'
        input_file = self.get_data(self.test_dir, input)
        self.get_data(self.test_dir, 'jwst_nirspec_shutters_barshadow.fits')

        step = Spec2Pipeline()
        step.output_file='jw95065_nrs_msaspec_barshadow_cal.fits'
        step.save_bsub = False
        step.save_results = True
        step.resample_spec.save_results = True
        step.extract_1d.save_results = True
        step.run(input_file)

        outputs = [('jw95065_nrs_msaspec_barshadow_cal.fits',
                    'jw95065_nrs_msaspec_barshadow_cal_ref.fits'),
                   ('jw95065_nrs_msaspec_barshadow_s2d.fits',
                    'jw95065_nrs_msaspec_barshadow_s2d_ref.fits'),
                   ('jw95065_nrs_msaspec_barshadow_x1d.fits',
                    'jw95065_nrs_msaspec_barshadow_x1d_ref.fits')
                  ]
        self.compare_outputs(outputs)
