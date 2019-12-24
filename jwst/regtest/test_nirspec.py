"""Regression tests for NIRSpec"""
import os

import pytest

from . import regtestdata as rt


@pytest.fixture(scope='module')
def run_nrs_ifu(jail, rtdata_module):
    """Run the pipelines"""
    step_params = {
        'input_path': 'nirspec/spectroscopic/nrs_ifu_uncal.fits',
        'step': 'calwebb_detector1.cfg',
        'args': [
            '--steps.group_scale.save_results=true',
            '--steps.dq_init.save_results=true',
            '--steps.saturation.save_results=true',
            '--steps.ipc.save_results=true',
            '--steps.superbias.save_results=true',
            '--steps.refpix.save_results=true',
            '--steps.rscd.save_results=true',
            '--steps.firstframe.save_results=true',
            '--steps.lastframe.save_results=true',
            '--steps.linearity.save_results=true',
            '--steps.dark_current.save_results=true',
            '--steps.persistence.save_results=true',
            '--steps.jump.save_results=true',
            '--steps.ramp_fit.save_results=true',
            '--steps.gain_scale.save_results=true',
        ]
    }

    return rt.run_step_from_dict(rtdata_module, **step_params)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ['0_ramp_fit', 'dark_current', 'dq_init', 'gain_scale', 'jump', 'linearity', 'rate', 'refpix',
     'saturation', 'superbias']
)
def test_nrs_ifu(run_nrs_ifu, fitsdiff_default_kwargs, suffix):
    """Regression test matching output files"""
    rt.is_like_truth(run_nrs_ifu, fitsdiff_default_kwargs, suffix,
                     truth_path='truth/nirspec/test_nirspec')
