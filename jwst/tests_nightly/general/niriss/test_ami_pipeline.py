import os
import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_ami3 import Ami3Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_ami_pipeline(_bigdata):
    """
    Regression test of the AMI pipeline performed on NIRISS AMI data.
    """
    pipe = Ami3Pipeline()
    pipe.save_averages = True
    pipe.ami_analyze.oversample = 3
    pipe.ami_analyze.rotation = 1.49
    pipe.run(_bigdata + '/niriss/test_ami_pipeline/test_lg1_asn.json')

    h = fits.open('test_targ_aminorm.fits')
    href = fits.open(_bigdata+'/niriss/test_ami_pipeline/ami_pipeline_targ_lgnorm.fits')

    result = fits.diff.FITSDiff(h, href,
                              ignore_hdus=['ASDF', 'HDRTAB'],
                              ignore_keywords=['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.001
    )

    assert result.identical, result.report()
