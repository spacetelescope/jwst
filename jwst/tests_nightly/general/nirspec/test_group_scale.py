import os
import pytest
from astropy.io import fits as pf
from jwst.group_scale.group_scale_step import GroupScaleStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_group_scale_nirspec(_bigdata):
    """

    Regression test of group scale step performed on NIRSpec data.

    """
    output_file_base, output_file = add_suffix('groupscale1_output.fits', 'group_scale')

    try:
        os.remove(output_file)
    except:
        pass

    GroupScaleStep.call(_bigdata+'/nirspec/test_group_scale/NRSIRS2_230_491_uncal.fits',
                        output_file=output_file_base, name='group_scale'
                        )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/nirspec/test_group_scale/NRSIRS2_230_491_groupscale.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
