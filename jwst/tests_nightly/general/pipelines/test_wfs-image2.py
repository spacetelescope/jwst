"""Test WFS&C calibrations"""

from glob import glob
from os import path as op
import pytest

from astropy.io import fits

from jwst.lib.suffix import replace_suffix
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_wfs_image2(_bigdata):
    """
    Regression test of the WFS&C `calwebb_wfs-image2.cfg` pipeline
    """

    data_path = op.join(_bigdata, 'nircam', 'test_datasets', 'sdp_jw82600_wfs', 'level2a')
    data_base = 'jw82600026001_02101_00001_nrca1_rate'
    ext = '.fits'

    collect_pipeline_cfgs('cfgs')
    Step.from_cmdline([op.join('cfgs', 'calwebb_wfs-image2.cfg'), op.join(data_path, data_base + ext)])

    output_files = glob('*')
    output_files.remove('cfgs')

    cal_name = replace_suffix(data_base, 'cal') + ext
    assert cal_name in output_files
    output_files.remove(cal_name)
    assert not output_files, 'Unexpected output files {}'.format(output_files)

    na = cal_name
    h = fits.open(na)
    nb = op.join(data_path, replace_suffix(data_base, 'cal_ref')) + ext
    href = fits.open(nb)
    newh = fits.HDUList([h['primary'], h['sci'], h['err'], h['dq'], h['area']])
    newhref = fits.HDUList([href['primary'], href['sci'], href['err'], href['dq'], href['area']])
    result = fits.diff.FITSDiff(
        newh,
        newhref,
        ignore_keywords=['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX'],
        rtol=0.00001
    )
    assert result.identical, result.report()
