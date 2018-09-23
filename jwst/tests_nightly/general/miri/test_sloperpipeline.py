from glob import glob

from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
from jwst.tests.base_test import MIRITest


class TestMIRISloperPipeline(MIRITest):
    ref_loc = ['test_sloperpipeline','truth']
    test_dir = 'test_sloperpipeline'
    
    def test_gain_scale_naming(self):
        """
        Regression test for gain_scale naming when results are requested to
        be saved for the gain_scale step.
        """
        expfile = 'jw00001001001_01101_00001_MIRIMAGE'
        input_file = self.get_data(self.test_dir, expfile+'_uncal.fits')

        step = Detector1Pipeline()
        step.group_scale.skip = True
        step.dq_init.skip = True
        step.saturation.skip = True
        step.ipc.skip = True
        step.superbias.skip = True
        step.refpix.skip = True
        step.rscd.skip = True
        step.firstframe.skip = True
        step.lastframe.skip = True
        step.linearity.skip = True
        step.dark_current.skip = True
        step.persistence.skip = True
        step.jump.skip = True
        step.ramp_fit.skip = False

        step.gain_scale.skip = False
        step.gain_scale.save_results = True

        step.run(input_file)


        files = glob('*.fits')

        output_file = expfile + '_gain_scale.fits'
        assert output_file in files
        files.remove(output_file)

        output_file = expfile + '_gain_scaleints.fits'
        assert output_file in files
        files.remove(output_file)

        assert not len(files)     
        
    def test_detector1pipeline1(self):
        """
        Regression test of calwebb_detector1 pipeline performed on MIRI data.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw00001001001_01101_00001_MIRIMAGE_uncal.fits')
                                   
        step = Detector1Pipeline()
        step.save_calibrated_ramp = True
        step.ipc.skip = True
        step.refpix.odd_even_columns = True
        step.refpix.use_side_ref_pixels = True
        step.refpix.side_smoothing_length=11
        step.refpix.side_gain=1.0
        step.refpix.odd_even_rows = True
        step.persistence.skip = True
        step.jump.rejection_threshold = 250.0
        step.ramp_fit.save_opt = False
        step.output_file='jw00001001001_01101_00001_MIRIMAGE'
        step.suffix='rate'

        step.run(input_file)
   
        outputs = [('jw00001001001_01101_00001_MIRIMAGE_ramp.fits',
                    'jw00001001001_01101_00001_MIRIMAGE_uncal_jump.fits'),
                   ('jw00001001001_01101_00001_MIRIMAGE_rateints.fits',
                    'jw00001001001_01101_00001_MIRIMAGE_uncal_integ.fits'),
                   ('jw00001001001_01101_00001_MIRIMAGE_rate.fits',
                    'jw00001001001_01101_00001_MIRIMAGE_uncal_MiriSloperPipeline.fits')
                   ]
        self.compare_outputs(outputs)

    def test_detector1pipeline2(self):
        """
        Regression test of calwebb_detector1 pipeline performed on MIRI data.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw80600012001_02101_00003_mirimage_uncal.fits')
                                   
        step = Detector1Pipeline()
        step.save_calibrated_ramp = True
        step.ipc.skip = True
        step.refpix.odd_even_columns = True
        step.refpix.use_side_ref_pixels = True
        step.refpix.side_smoothing_length=11
        step.refpix.side_gain=1.0
        step.refpix.odd_even_rows = True
        step.persistence.skip = True
        step.jump.rejection_threshold = 250.0
        step.ramp_fit.save_opt = False
        step.output_file='jw80600012001_02101_00003_mirimage'
        step.suffix='rate'

        step.run(input_file)
        
        outputs = [('jw80600012001_02101_00003_mirimage_ramp.fits',
                    'jw80600012001_02101_00003_mirimage_ramp.fits'),
                   ('jw80600012001_02101_00003_mirimage_rateints.fits',
                    'jw80600012001_02101_00003_mirimage_rateints.fits'),
                   ('jw80600012001_02101_00003_mirimage_rate.fits',
                    'jw80600012001_02101_00003_mirimage_rate.fits')
                   ]
        self.compare_outputs(outputs)

