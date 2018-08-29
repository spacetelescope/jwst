from glob import glob
import os
import pytest

from astropy.io import fits as pf

from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_detector1pipeline1(_bigdata):
    """
    Regression test for gain_scale naming when results are requested to
    be saved for the gain_scale step.
    """

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

    expfile = 'jw00001001001_01101_00001_MIRIMAGE'
    step.run(_bigdata+'/miri/test_sloperpipeline/' + expfile + '_uncal.fits')


    files = glob('*.fits')

    output_file = expfile + '_gain_scale.fits'
    assert output_file in files
    files.remove(output_file)

    output_file = expfile + '_gain_scaleints.fits'
    assert output_file in files
    files.remove(output_file)

    assert not len(files)
