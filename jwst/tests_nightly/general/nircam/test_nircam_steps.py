import pytest

from jwst.tests.base_classes import BaseJWSTTestSteps
from jwst.tests.base_classes import pytest_generate_tests # noqa: F401

from jwst.refpix import RefPixStep
from jwst.dark_current import DarkCurrentStep
from jwst.dq_init import DQInitStep
from jwst.emission import EmissionStep
from jwst.flatfield import FlatFieldStep
from jwst.ipc import IPCStep
from jwst.jump import JumpStep
from jwst.linearity import LinearityStep
from jwst.persistence import PersistenceStep
from jwst.photom import PhotomStep
from jwst.saturation import SaturationStep


# Parameterized regression tests for NIRCAM processing
# All tests in this set run with 1 input file and
#  only generate 1 output for comparison.
#
@pytest.mark.bigdata
class TestNIRCamSteps(BaseJWSTTestSteps):
    input_loc = 'nircam'

    params = {'test_steps':
                [dict(input='jw00017001001_01101_00001_NRCA1_dq_init.fits',
                      test_dir='test_bias_drift',
                      step_class=RefPixStep,
                      step_pars=dict(odd_even_columns=True,
                                     use_side_ref_pixels=False,
                                     side_smoothing_length=10,
                                     side_gain=1.0),
                      output_truth='jw00017001001_01101_00001_NRCA1_bias_drift.fits',
                      output_hdus=[],
                      id='refpix_nircam'

                      ),
                 dict(input='jw00017001001_01101_00001_NRCA1_saturation.fits',
                      test_dir='test_dark_step',
                      step_class=DarkCurrentStep,
                      step_pars=dict(),
                      output_truth='jw00017001001_01101_00001_NRCA1_dark_current.fits',
                      output_hdus=[],
                      id='dark_current_nircam'

                      ),
                 dict(input='jw00017001001_01101_00001_NRCA1_uncal.fits',
                      test_dir='test_dq_init',
                      step_class=DQInitStep,
                      step_pars=dict(),
                      output_truth='jw00017001001_01101_00001_NRCA1_dq_init.fits',
                      output_hdus=[],
                      id='dq_init_nircam'
                      ),
                 dict(input='jw00017001001_01101_00001_NRCA1_persistence.fits',
                      test_dir='test_emission',
                      step_class=EmissionStep,
                      step_pars=dict(),
                      output_truth='jw00017001001_01101_00001_NRCA1_emission.fits',
                      output_hdus=[],
                      id='emission_nircam'
                      ),
                 dict(input='jw00017001001_01101_00001_NRCA1_ramp_fit.fits',
                      test_dir='test_flat_field',
                      step_class=FlatFieldStep,
                      step_pars=dict(),
                      output_truth='jw00017001001_01101_00001_NRCA1_flat_field.fits',
                      output_hdus=[],
                      id='flat_field_nircam'
                      ),
                 dict(input='jw00017001001_01101_00001_NRCA3_uncal.fits',
                      test_dir='test_ipc_step',
                      step_class=IPCStep,
                      step_pars=dict(),
                      output_truth='jw00017001001_01101_00001_NRCA3_ipc.fits',
                      output_hdus=['primary', 'sci'],
                      id='ipc_nircam'
                      ),
                 dict(input='jw00017001001_01101_00001_NRCA1_linearity.fits',
                      test_dir='test_jump',
                      step_class=JumpStep,
                      step_pars=dict(rejection_threshold=25.0),
                      output_truth='jw00017001001_01101_00001_NRCA1_jump.fits',
                      output_hdus=[],
                      id='jump_nircam'
                      ),
                 dict(input='jw00017001001_01101_00001_NRCA1_dark_current.fits',
                      test_dir='test_linearity',
                      step_class=LinearityStep,
                      step_pars=dict(),
                      output_truth='jw00017001001_01101_00001_NRCA1_linearity.fits',
                      output_hdus=[],
                      id='linearity_nircam'
                      ),
                 dict(input='jw00017001001_01101_00001_NRCA1_ramp.fits',
                      test_dir='test_persistence',
                      step_class=PersistenceStep,
                      step_pars=dict(),
                      output_truth='jw00017001001_01101_00001_NRCA1_persistence.fits',
                      output_hdus=[],
                      id='persistence_nircam'
                      ),
                 dict(input='jw00017001001_01101_00001_NRCA1_emission.fits',
                      test_dir='test_photom',
                      step_class=PhotomStep,
                      step_pars=dict(),
                      output_truth='jw00017001001_01101_00001_NRCA1_photom.fits',
                      output_hdus=[],
                      id='photom_nircam'
                      ),
                 dict(input='jw00017001001_01101_00001_NRCA1_bias_drift.fits',
                      test_dir='test_saturation',
                      step_class=SaturationStep,
                      step_pars=dict(),
                      output_truth='jw00017001001_01101_00001_NRCA1_saturation.fits',
                      output_hdus=[],
                      id='saturation_nircam'
                      ),
                ]
              }
