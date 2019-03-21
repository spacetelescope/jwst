import pytest

from jwst.tests.base_classes import BaseJWSTTestSteps
from jwst.tests.base_classes import pytest_generate_tests # noqa: F401

from jwst.refpix import RefPixStep
from jwst.dark_current import DarkCurrentStep
from jwst.dq_init import DQInitStep
from jwst.extract_1d import Extract1dStep
from jwst.extract_2d import Extract2dStep
from jwst.flatfield import FlatFieldStep
from jwst.group_scale import GroupScaleStep
from jwst.jump import JumpStep
from jwst.linearity import LinearityStep
from jwst.photom import PhotomStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep


# Parameterized regression tests for NIRISS processing
# All tests in this set run with 1 input file and
#  only generate 1 output for comparison.
#
@pytest.mark.bigdata
class TestNIRSpecSteps(BaseJWSTTestSteps):
    input_loc = 'nirspec'

    params = {'test_steps':
                [dict(input='jw00023001001_01101_00001_NRS1_dq_init.fits',
                      test_dir='test_bias_drift',
                      step_class=RefPixStep,
                      step_pars=dict(odd_even_columns=True,
                                     use_side_ref_pixels=False,
                                     side_smoothing_length=10,
                                     side_gain=1.0),
                      output_truth='jw00023001001_01101_00001_NRS1_bias_drift.fits',
                      output_hdus=[],
                      id='refpix_nirspec'
                      ),
                 dict(input='jw00023001001_01101_00001_NRS1_saturation.fits',
                      test_dir='test_dark_step',
                      step_class=DarkCurrentStep,
                      step_pars=dict(),
                      output_truth='jw00023001001_01101_00001_NRS1_dark_current.fits',
                      output_hdus=[],
                      id='dark_current_nirspec'
                      ),
                 dict(input='jw00023001001_01101_00001_NRS1_uncal.fits',
                      test_dir='test_dq_init',
                      step_class=DQInitStep,
                      step_pars=dict(),
                      output_truth='jw00023001001_01101_00001_NRS1_dq_init.fits',
                      output_hdus=[],
                      id='dq_init_nirspec'
                      ),
                 dict(input='jw00023001001_01101_00001_NRS1_cal.fits',
                      test_dir='test_extract_1d',
                      step_class=Extract1dStep,
                      step_pars=dict(),
                      output_truth='jw00023001001_01101_00001_NRS1_spec.fits',
                      output_hdus=['primary',('extract1d',1),('extract1d',2),
                                   ('extract1d',3),('extract1d',4)],
                      id='extract1d_nirspec'
                      ),
                 dict(input='jw00023001001_01101_00001_NRS1_assign_wcs.fits',
                      test_dir='test_extract_2d',
                      step_class=Extract2dStep,
                      step_pars=dict(),
                      output_truth='jw00023001001_01101_00001_NRS1_extract_2d.fits',
                      output_hdus=['primary', ('sci', 1), ('err', 1), ('dq', 1),
                                         ('sci', 2), ('err', 2), ('dq', 2),
                                         ('sci', 3), ('err', 3), ('dq', 3),
                                         ('sci', 4), ('err', 4), ('dq', 4),
                                         ('sci', 5), ('err', 5), ('dq', 5)],
                      id='extract2d_nirspec'
                      ),
                 dict(input='jw00023001001_01101_00001_NRS1_extract_2d.fits',
                      test_dir='test_flat_field',
                      step_class=FlatFieldStep,
                      step_pars=dict(),
                      output_truth='jw00023001001_01101_00001_NRS1_flat_field.fits',
                      output_hdus=['primary', ('sci', 1), ('err', 1), ('dq', 1),
                                         ('sci', 2), ('err', 2), ('dq', 2),
                                         ('sci', 3), ('err', 3), ('dq', 3),
                                         ('sci', 4), ('err', 4), ('dq', 4)],
                      id='flat_field_nirspec'
                      ),
                 dict(input='NRSIRS2_230_491_uncal.fits',
                      test_dir='test_group_scale',
                      step_class=GroupScaleStep,
                      step_pars=dict(),
                      output_truth='NRSIRS2_230_491_groupscale.fits',
                      output_hdus=[],
                      id='group_scale_nirspec'
                      ),
                 dict(input='jw00023001001_01101_00001_NRS1_linearity.fits',
                      test_dir='test_jump',
                      step_class=JumpStep,
                      step_pars=dict(rejection_threshold=50.0),
                      output_truth='jw00023001001_01101_00001_NRS1_jump.fits',
                      output_hdus=[],
                      id='jump_nirspec'
                      ),
                 dict(input='jw00023001001_01101_00001_NRS1_dark_current.fits',
                      test_dir='test_linearity',
                      step_class=LinearityStep,
                      step_pars=dict(),
                      output_truth='jw00023001001_01101_00001_NRS1_linearity.fits',
                      output_hdus=[],
                      id='linearity_nirspec'
                      ),
                 dict(input='jw00023001001_01101_00001_NRS1_flat_field.fits',
                      test_dir='test_photom',
                      step_class=PhotomStep,
                      step_pars=dict(),
                      output_truth='jw00023001001_01101_00001_NRS1_photom.fits',
                      output_hdus=['primary',('sci',1),('err',1),('dq',1),('relsens',1),
                                              ('sci',2),('err',2),('dq',2),('relsens',2),
                                              ('sci',3),('err',3),('dq',3),('relsens',3),
                                              ('sci',4),('err',4),('dq',4),('relsens',4)],
                      id='photom_nirspec'
                      ),
                 dict(input='jw84600007001_02101_00001_nrs1_superbias.fits',
                      test_dir='test_bias_drift',
                      step_class=RefPixStep,
                      step_pars=dict(),
                      output_truth='jw84600007001_02101_00001_nrs1_refpix.fits',
                      output_hdus=[],
                      id='refpix_nirspec_irs2'
                      ),
                 dict(input='jw00023001001_01101_00001_NRS1_bias_drift.fits',
                      test_dir='test_saturation',
                      step_class=SaturationStep,
                      step_pars=dict(),
                      output_truth='jw00023001001_01101_00001_NRS1_saturation.fits',
                      output_hdus=[],
                      id='saturation_nirspec'
                      ),
                 dict(input='jw00011001001_01106_00001_NRS2_saturation.fits',
                      test_dir='test_superbias',
                      step_class=SuperBiasStep,
                      step_pars=dict(),
                      output_truth='jw00011001001_01106_00001_NRS2_superbias.fits',
                      output_hdus=[],
                      id='superbias_nirspec'
                      )
                ]
              }
