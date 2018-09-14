"""Test calwebb_spec3 against NIRSpec Fixed-slit science (FSS)"""
from glob import glob
from os import path
import pytest

from astropy.io import fits

from jwst.associations import load_asn
from jwst.pipeline import Spec3Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


@pytest.mark.xfail(
    reason='Input data not available',
    run=False
)
def test_save_source_only(_bigdata):
    """Test saving the source-based files only"""
    datapath = path.join(
        _bigdata, 'nirspec', 'test_datasets', 'fss', '93045', 'level2b'
    )
    pipe = Spec3Pipeline()
    pipe.mrs_imatch.skip = True
    pipe.outlier_detection.skip = True
    pipe.resample_spec.skip = True
    pipe.cube_build.skip = True
    pipe.extract_1d.skip = True

    asn_path = path.join(datapath, 'jw93045-o010_20180725t035735_spec3_001_asn.json') 
    pipe.run(asn_path)

    # Check resulting product
    with open(asn_path) as fh:
        asn = load_asn(fh)
    base_name = asn['products'][0]['name']
    product_name = base_name.format(source_id='ss400a1') + '_cal.fits'
    output_files = glob('*')

    if product_name in output_files:
        output_files.remove(product_name)
    else:
        assert False

    assert len(output_files) == 0


@pytest.mark.xfail(
    reason='Dataset fails at resample',
    run=False
)
def test_nrs_fs_spec3(_bigdata):
    """
    Regression test of calwebb_spec3 pipeline performed on
    NIRSpec fixed-slit data.
    """
    cfg_dir = './cfgs'
    collect_pipeline_cfgs(cfg_dir)
    datapath = path.join(
        _bigdata, 'nirspec', 'test_datasets', 'fss', '93045', 'level2b'
    )
    args = [
        path.join(cfg_dir, 'calwebb_spec3.cfg'),
        path.join(datapath, 'jw93045-o010_20180725t035735_spec3_001_asn.json')
    ]

    Step.from_cmdline(args)

    # Compare results
    assert False

    ignore_keywords = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX']

    na = 'jw00023001001_01101_00001_NRS1_cal.fits'
    nb = _bigdata+'/pipelines/jw00023001001_01101_00001_NRS1_cal_ref.fits'
    h = fits.open(na)
    href = fits.open(nb)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol = 0.00001)
    assert result.identical, result.report()

    na = 'jw00023001001_01101_00001_NRS1_s2d.fits'
    nb = _bigdata+'/pipelines/jw00023001001_01101_00001_NRS1_s2d_ref.fits'
    h = fits.open(na)
    href = fits.open(nb)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol = 0.00001)
    assert result.identical, result.report()

    na = 'jw00023001001_01101_00001_NRS1_x1d.fits'
    nb = _bigdata+'/pipelines/jw00023001001_01101_00001_NRS1_x1d_ref.fits'
    h = fits.open(na)
    href = fits.open(nb)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol = 0.00001)
    assert result.identical, result.report()
