import pytest

from jwst.tests.base_test import MIRITestSteps
from jwst.tests.base_test import pytest_generate_tests

from jwst.refpix.refpix_step import RefPixStep
from jwst.dark_current.dark_current_step import DarkCurrentStep
from jwst.dq_init.dq_init_step import DQInitStep
from jwst.emission.emission_step import EmissionStep
from jwst.extract_1d.extract_1d_step import Extract1dStep
from jwst.flatfield.flat_field_step import FlatFieldStep
from jwst.fringe.fringe_step import FringeStep
from jwst.jump.jump_step import JumpStep
from jwst.lastframe.lastframe_step import LastFrameStep
from jwst.linearity.linearity_step import LinearityStep
from jwst.photom.photom_step import PhotomStep
from jwst.rscd.rscd_step import RSCD_Step
from jwst.saturation.saturation_step import SaturationStep
from jwst.srctype.srctype_step import SourceTypeStep
from jwst.straylight.straylight_step import StraylightStep

# Parameterized regression tests for MIRI processing
# All tests in this set run with 1 input file and 
#  only generate 1 output for comparison.
#
@pytest.mark.bigdata
class TestMIRISteps(MIRITestSteps):

    params = {'test_steps': 
                [
                # test_refpix_miri: refpix step performed on MIRI data
                 dict(input='jw00001001001_01101_00001_MIRIMAGE_saturation.fits',
                      test_dir='test_bias_drift',
                      step_class=RefPixStep, 
                      step_pars=dict(use_side_ref_pixels=False, 
                                     side_smoothing_length=10, 
                                     side_gain=1.0),
                      output_truth='jw00001001001_01101_00001_MIRIMAGE_bias_drift.fits',
                      output_hdus=[],
                      id='test_refpix_miri'
                     ),
                # test_refpix_miri2: refpix step performed on MIRI data
                 dict(input='jw00025001001_01107_00001_MIRIMAGE_saturation.fits',
                      test_dir='test_bias_drift',
                      step_class=RefPixStep, 
                      step_pars=dict(use_side_ref_pixels=False, 
                                     side_smoothing_length=10, 
                                     side_gain=1.0),
                      output_truth='jw00025001001_01107_00001_MIRIMAGE_bias_drift.fits',
                      output_hdus=[],
                      id='test_refpix_miri2'
                     ),
                # test_dark_current_miri: dark current step performed on MIRI data
                 dict(input='jw00001001001_01101_00001_MIRIMAGE_bias_drift.fits',
                      test_dir='test_dark_step',
                      step_class=DarkCurrentStep, 
                      step_pars=dict(),
                      output_truth='jw00001001001_01101_00001_MIRIMAGE_dark_current.fits',
                      output_hdus=[],
                      id='test_dark_current_miri'
                     ),
                # test_dark_current_miri2: dark current step performed on MIRI data
                 dict(input='jw80600012001_02101_00003_mirimage_lastframe.fits',
                      test_dir='test_dark_step',
                      step_class=DarkCurrentStep, 
                      step_pars=dict(),
                      output_truth='jw80600012001_02101_00003_mirimage_dark.fits',
                      output_hdus=[],
                      id='test_dark_current_miri2'
                     ),
                # test_dq_init_miri: dq_init step performed on uncalibrated MIRI data
                 dict(input='jw00001001001_01101_00001_MIRIMAGE_uncal.fits',
                      test_dir='test_dq_init',
                      step_class=DQInitStep, 
                      step_pars=dict(),
                      output_truth='jw00001001001_01101_00001_MIRIMAGE_dq_init.fits',
                      output_hdus=[],
                      id='test_dq_init_miri'
                     ),
                # test_dq_init_miri2: dq_init step performed on uncalibrated MIRI data
                 dict(input='jw80600012001_02101_00003_mirimage_uncal.fits',
                      test_dir='test_dq_init',
                      step_class=DQInitStep, 
                      step_pars=dict(),
                      output_truth='jw80600012001_02101_00003_mirimage_dqinit.fits',
                      output_hdus=[],
                      id='test_dq_init_miri2'
                     ),
                # test_emission_miri: emission step performed on calibrated miri data
                 dict(input='jw00001001001_01101_00001_MIRIMAGE_flat_field.fits',
                      test_dir='test_emission',
                      step_class=EmissionStep, 
                      step_pars=dict(),
                      output_truth='jw00001001001_01101_00001_MIRIMAGE_emission.fits',
                      output_hdus=[],
                      id='test_emission_miri'
                     ),
                # test_extract1d_miri: extract_1d step performed on MIRI LRS fixed-slit data
                 dict(input='jw00035001001_01101_00001_mirimage_photom.fits',
                      test_dir='test_extract1d',
                      step_class=Extract1dStep, 
                      step_pars=dict(smoothing_length=0),
                      output_truth='jw00035001001_01101_00001_mirimage_x1d.fits',
                      output_hdus=[],
                      id='test_extract1d_miri'
                     ),
                # test_extract1d_miri2: extract_1d step performed on MIRI LRS slitless data
                 dict(input='jw80600012001_02101_00003_mirimage_photom.fits',
                      test_dir='test_extract1d',
                      step_class=Extract1dStep, 
                      step_pars=dict(smoothing_length=0),
                      output_truth='jw80600012001_02101_00003_mirimage_x1d.fits',
                      output_hdus=['primary',('extract1d',1),('extract1d',2),('extract1d',3),('extract1d',4)],
                      id='test_extract1d_miri2'
                     ),                
                # test_flat_field_miri: flat_field step performed on MIRI data.
                 dict(input='jw00001001001_01101_00001_MIRIMAGE_assign_wcs.fits',
                      test_dir='test_flat_field',
                      step_class=FlatFieldStep, 
                      step_pars=dict(),
                      output_truth='jw00001001001_01101_00001_MIRIMAGE_flat_field.fits',
                      output_hdus=[],
                      id='test_flat_field_miri'
                     ),
                # test_flat_field_miri2: flat_field step performed on MIRI data.
                 dict(input='jw80600012001_02101_00003_mirimage_assign_wcs.fits',
                      test_dir='test_flat_field',
                      step_class=FlatFieldStep, 
                      step_pars=dict(),
                      output_truth='jw80600012001_02101_00003_mirimage_flat_field.fits',
                      output_hdus=[],
                      id='test_flat_field_miri2'
                     ),
                # test_fringe_miri: fringe performed on MIRI data.
                 dict(input='fringe1_input.fits',
                      test_dir='test_fringe',
                      step_class=FringeStep, 
                      step_pars=dict(),
                      output_truth='baseline_fringe1.fits',
                      output_hdus=['primary','sci','err','dq'],
                      id='test_fringe_miri'
                     ),
                # test_fringe_miri2: fringe performed on MIRI data.
                 dict(input='fringe2_input.fits',
                      test_dir='test_fringe',
                      step_class=FringeStep, 
                      step_pars=dict(),
                      output_truth='baseline_fringe2.fits',
                      output_hdus=['primary','sci','err','dq'],
                      id='test_fringe_miri2'
                     ),
                # test_fringe_miri3: fringe performed on MIRI data.
                 dict(input='fringe3_input.fits',
                      test_dir='test_fringe',
                      step_class=FringeStep, 
                      step_pars=dict(),
                      output_truth='baseline_fringe3.fits',
                      output_hdus=['primary','sci','err','dq'],
                      id='test_fringe_miri3'
                     ),
                # test_jump_miri: jump step performed on MIRI data.
                 dict(input='jw00001001001_01101_00001_MIRIMAGE_linearity.fits',
                      test_dir='test_jump',
                      step_class=JumpStep, 
                      step_pars=dict(rejection_threshold=200.0),
                      output_truth='jw00001001001_01101_00001_MIRIMAGE_jump.fits',
                      output_hdus=[],
                      id='test_jump_miri'
                     ),
                # test_jump_miri2: jump step performed on MIRI data.
                 dict(input='jw80600012001_02101_00003_mirimage_dark.fits',
                      test_dir='test_jump',
                      step_class=JumpStep, 
                      step_pars=dict(rejection_threshold=25.0),
                      output_truth='jw80600012001_02101_00003_mirimage_jump.fits',
                      output_hdus=[],
                      id='test_jump_miri2'
                     ),
                # test_lastframe_miri2: lastframe step performed on MIRI data
                 dict(input='jw80600012001_02101_00003_mirimage_rscd.fits',
                      test_dir='test_lastframe',
                      step_class=LastFrameStep, 
                      step_pars=dict(),
                      output_truth='jw80600012001_02101_00003_mirimage_lastframe.fits',
                      output_hdus=[],
                      id='test_lastframe_miri2'
                     ),
                # test_linearity_miri: linearity step performed on MIRI data
                 dict(input='jw00001001001_01101_00001_MIRIMAGE_dark_current.fits',
                      test_dir='test_linearity',
                      step_class=LinearityStep, 
                      step_pars=dict(),
                      output_truth='jw00001001001_01101_00001_MIRIMAGE_linearity.fits',
                      output_hdus=[],
                      id='test_linearity_miri'
                     ),
                # test_linearity_miri2: linearity step performed on MIRI data
                 dict(input='jw80600012001_02101_00003_mirimage_saturation.fits',
                      test_dir='test_linearity',
                      step_class=LinearityStep, 
                      step_pars=dict(),
                      output_truth='jw80600012001_02101_00003_mirimage_linearity.fits',
                      output_hdus=[],
                      id='test_linearity_miri2'
                     ),
                # test_photom_miri: photom step performed on MIRI imaging data
                 dict(input='jw00001001001_01101_00001_MIRIMAGE_emission.fits',
                      test_dir='test_photom',
                      step_class=PhotomStep, 
                      step_pars=dict(),
                      output_truth='jw00001001001_01101_00001_MIRIMAGE_photom.fits',
                      output_hdus=[],
                      id='test_photom_miri'
                     ),
                 # test_photom_miri2: photom step performed on MIRI LRS slitless data
                 dict(input='jw80600012001_02101_00003_mirimage_srctype.fits',
                      test_dir='test_photom',
                      step_class=PhotomStep, 
                      step_pars=dict(),
                      output_truth='jw80600012001_02101_00003_mirimage_photom.fits',
                      output_hdus=[],
                      id='test_photom_miri2'
                     ),
                 # test_rscd_miri2: RSCD step performed on MIRI data
                 dict(input='jw80600012001_02101_00003_mirimage_linearity.fits',
                      test_dir='test_rscd',
                      step_class=RSCD_Step, 
                      step_pars=dict(),
                      output_truth='jw80600012001_02101_00003_mirimage_rscd.fits',
                      output_hdus=[],
                      id='test_rscd_miri2'
                     ),
                 # test_saturation_miri: saturation step performed on uncalibrated MIRI data
                 dict(input='jw00001001001_01101_00001_MIRIMAGE_dq_init.fits',
                      test_dir='test_saturation',
                      step_class=SaturationStep, 
                      step_pars=dict(),
                      output_truth='jw00001001001_01101_00001_MIRIMAGE_saturation.fits',
                      output_hdus=['primary','sci','err','pixeldq','groupdq'],
                      id='test_saturation_miri'
                     ),
                 # test_saturation_miri2: saturation step performed on uncalibrated MIRI data
                 dict(input='jw80600012001_02101_00003_mirimage_dqinit.fits',
                      test_dir='test_saturation',
                      step_class=SaturationStep, 
                      step_pars=dict(),
                      output_truth='jw80600012001_02101_00003_mirimage_saturation.fits',
                      output_hdus=[],
                      id='test_saturation_miri2'
                     ),
                 # test_srctype2: srctype step performed on MIRI LRS slitless data
                 dict(input='jw80600012001_02101_00003_mirimage_flat_field.fits',
                      test_dir='test_srctype',
                      step_class=SourceTypeStep, 
                      step_pars=dict(),
                      output_truth='jw80600012001_02101_00003_mirimage_srctype.fits',
                      output_hdus=[],
                      id='test_srctype2'
                     ),
                 # test_straylight1_miri: straylight performed on MIRI IFUSHORT data
                 dict(input='jw80500018001_02101_00002_MIRIFUSHORT_flatfield.fits',
                      test_dir='test_straylight',
                      step_class=StraylightStep, 
                      step_pars=dict(),
                      output_truth='jw80500018001_02101_00002_MIRIFUSHORT_straylight.fits',
                      output_hdus=['primary','sci','err','dq'],
                      id='test_straylight1_miri'
                     ),
                 # test_straylight2_miri: straylight performed on MIRI IFULONG data
                 dict(input='jw80500018001_02101_00002_MIRIFULONG_flatfield.fits',
                      test_dir='test_straylight',
                      step_class=StraylightStep, 
                      step_pars=dict(),
                      output_truth='jw80500018001_02101_00002_MIRIFULONG_straylight.fits',
                      output_hdus=[],
                      id='test_straylight2_miri'
                     ),
               ]
             }


