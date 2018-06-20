import os
import pytest
from astropy.io import fits as pf
from jwst.ami.ami_analyze_step import AmiAnalyzeStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_ami_analyze(_bigdata):
    """

    Regression test of ami_analyze step performed on NIRISS AMI data.

    """
    suffix = 'ami_analyze'
    output_file_base, output_file = add_suffix('ami_analyze_output_16.fits', suffix)

    AmiAnalyzeStep.call(_bigdata+'/niriss/test_ami_analyze/ami_analyze_input_16.fits',
                        oversample=3, rotation=1.49,
                        output_file=output_file_base, suffix=suffix
                        )

    h = pf.open(output_file)
    href = pf.open(_bigdata+'/niriss/test_ami_analyze/ami_analyze_ref_output_16.fits')
    newh = pf.HDUList([h['primary'],h['fit'],h['resid'],h['closure_amp'],
                       h['closure_pha'],h['fringe_amp'],h['fringe_pha'],
                       h['pupil_pha'],h['solns']])
    newhref = pf.HDUList([href['primary'],href['fit'],href['resid'],href['closure_amp'],
                          href['closure_pha'],href['fringe_amp'],href['fringe_pha'],
                          href['pupil_pha'],href['solns']])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.001
    )
    assert result.identical, result.report()
