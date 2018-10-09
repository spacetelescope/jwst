import os
import glob

import pytest

from jwst.tests.base_test import BaseJWSTTest
from jwst.ramp_fitting.ramp_fit_step import RampFitStep
from jwst.wfs_combine.wfs_combine_step import WfsCombineStep
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

from jwst.stpipe import Step


@pytest.mark.bigdata
class TestWFSImage2(BaseJWSTTest):
    input_loc = 'nircam'
    ref_loc = ['test_datasets', 'sdp_jw82600_wfs', 'level2a', 'truth']
    test_dir = 'test_datasets'

    def test_wfs_image2(self):
        """
        Regression test of the WFS&C `calwebb_wfs-image2.cfg` pipeline
        """
        data_base = 'jw82600026001_02101_00001_nrca1_rate'
        ext = '.fits'
        input_name = '{}{}'.format(data_base, ext)
        input_file = self.get_data(self.test_dir, 'sdp_jw82600_wfs', 'level2a',
                                   input_name)

        self.get_custom_cfgs(self.test_dir, 'sdp_jw82600_wfs', 'cfgs')
        Step.from_cmdline([os.path.join('cfgs', 'calwebb_wfs-image2.cfg'), input_file])

        output_files = glob.glob('*')
        output_files.remove('cfgs')
        try:
            # TEST_BIGDATA set to artifactory will copy it to working dir
            output_files.remove(input_name)
        except ValueError:
            # File won't exist if local cache is used
            pass
        
        cal_name = input_name.replace('rate', 'cal')
        output_name = input_name.replace('rate','cal_ref')
        outputs = [(cal_name, output_name)]
        self.compare_outputs(outputs)

        assert cal_name in output_files
        output_files.remove(cal_name)
        assert not output_files, 'Unexpected output files {}'.format(output_files)


@pytest.mark.bigdata
class TestDetector1Pipeline(BaseJWSTTest):
    input_loc = 'nircam'
    ref_loc = ['test_detector1pipeline', 'truth']
    test_dir = 'test_detector1pipeline'

    def test_detector1pipeline3(self):
        """
        Regression test of calwebb_detector1 pipeline performed on NIRCam data.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw82500001003_02101_00001_NRCALONG_uncal.fits')
        step = Detector1Pipeline()
        step.save_calibrated_ramp = True
        step.ipc.skip = True
        step.refpix.odd_even_columns = True
        step.refpix.use_side_ref_pixels = False
        step.refpix.side_smoothing_length = 10
        step.refpix.side_gain = 1.0
        step.persistence.skip = True
        step.jump.rejection_threshold = 250.0
        step.ramp_fit.save_opt = True
        step.output_file = 'jw82500001003_02101_00001_NRCALONG_rate.fits'
        step.run(input_file)

        outputs = [('jw82500001003_02101_00001_NRCALONG_ramp.fits',
                    'jw82500001003_02101_00001_NRCALONG_ramp_ref.fits'),
                   ('jw82500001003_02101_00001_NRCALONG_rate.fits',
                    'jw82500001003_02101_00001_NRCALONG_rate_ref.fits'),
                   ('jw82500001003_02101_00001_NRCALONG_rateints.fits',
                    'jw82500001003_02101_00001_NRCALONG_rateints_ref.fits')
                  ]
        self.compare_outputs(outputs)


@pytest.mark.bigdata
class TestNIRCamRamp(BaseJWSTTest):
    input_loc = 'nircam'
    ref_loc = ['test_ramp_fit', 'truth']
    test_dir = 'test_ramp_fit'

    def test_ramp_fit_nircam(self):
        """
        Regression test of ramp_fit step performed on NIRCam data.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw00017001001_01101_00001_NRCA1_jump.fits')

        result, result_int = RampFitStep.call(input_file,
                                  save_opt=True,
                                  opt_name='rampfit_opt_out.fits'
                                  )

        optout_file = 'rampfit_opt_out_fitopt.fits'
        output_file = result.meta.filename
        result.save(output_file)
        result.close()

        outputs = [(output_file,
                    'jw00017001001_01101_00001_NRCA1_ramp_fit.fits'),
                    (optout_file,
                     'rampfit_opt_out.fits',
                     ['primary','slope','sigslope','yint','sigyint','pedestal','weights','crmag']) ]
        self.compare_outputs(outputs)


@pytest.mark.bigdata
class TestWFSCombine(BaseJWSTTest):
    input_loc = 'nircam'
    ref_loc = ['test_wfs_combine', 'truth']
    test_dir = 'test_wfs_combine'

    def test_wfs_combine(self):
        """
        Regression test of wfs_combine using do_refine=False (default)
        Association table has 3 (identical) pairs of input files to combine
        """
        asn_file = self.get_data(self.test_dir,
                                   'wfs_3sets_asn.json')
        for filename in self.raw_from_asn(asn_file):
            self.get_data(self.test_dir, filename)

        WfsCombineStep.call(asn_file)

        outputs = [('test_wfscom_wfscmb.fits',
                    'test_wfscom.fits'),
                   ('test_wfscoma_wfscmb.fits',
                    'test_wfscoma.fits'),
                   ('test_wfscomb_wfscmb.fits',
                    'test_wfscomb.fits')
                  ]
        self.compare_outputs(outputs)

    def test_wfs_combine1(self):
        """
        Regression test of wfs_combine using do_refine=True
        """
        asn_file = self.get_data(self.test_dir,
                                   'wfs_3sets_asn2.json')
        for filename in self.raw_from_asn(asn_file):
            self.get_data(self.test_dir, filename)

        WfsCombineStep.call(asn_file,
                            do_refine=True )

        outputs = [('test_wfscom2_wfscmb.fits',
                    'test_wfscom_do_ref.fits'),
                    ('test_wfscom2a_wfscmb.fits',
                     'test_wfscoma_do_ref.fits'),
                    ('test_wfscom2b_wfscmb.fits',
                     'test_wfscomb_do_ref.fits')
                   ]
        self.compare_outputs(outputs)

    def test_wfs_combine2(self):
        """
        Regression test of wfs_combine using do_refine=True
        """
        asn_file = self.get_data(self.test_dir,
                                   'wfs_3sets_asn3.json')
        for filename in self.raw_from_asn(asn_file):
            self.get_data(self.test_dir, filename)

        WfsCombineStep.call(asn_file,
                            do_refine=True)

        outputs = [('test_wfscom3_wfscmb.fits',
                    'test_wfscom_do_ref.fits'),
                   ('test_wfscom3a_wfscmb.fits',
                    'test_wfscoma_do_ref.fits'),
                   ('test_wfscom3b_wfscmb.fits',
                    'test_wfscomb_do_ref.fits')
                  ]
        self.compare_outputs(outputs)
