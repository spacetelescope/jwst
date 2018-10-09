import os
import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_coron3 import Coron3Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_coron3_pipeline1(_bigdata):
    """Regression test of calwebb_coron3 pipeline.

    Test will be performed on NIRCam simulated data.
    """
    subdir = os.path.join(_bigdata, 'nircam', 'test_coron3')
    asn_name = 'jw99999-a3001_20170327t121212_coron3_001_asn.json'
    asn_file = os.path.join(subdir, asn_name)
    override_psfmask = os.path.join(subdir, 'jwst_nircam_psfmask_somb.fits')

    ignore_keywords = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX']

    pipe = Coron3Pipeline()
    pipe.align_refs.override_psfmask = override_psfmask
    pipe.outlier_detection.resample_data = False
    pipe.run(asn_file)

    # Compare psfstack product
    n_cur = 'jw99999-a3001_t1_nircam_f140m-maskbar_psfstack.fits'
    n_ref_name = 'jw99999-a3001_t1_nircam_f140m-maskbar_psfstack_ref.fits'
    n_ref = os.path.join(subdir, n_ref_name)
    h = fits.open(n_cur)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.001)
    assert result.identical, result.report()

    # Compare psfalign product
    n_cur = 'jw9999947001_02102_00001_nrcb3_a3001_psfalign.fits'
    n_ref_name = 'jw99999-a3001_t1_nircam_f140m-maskbar_psfalign_ref.fits'
    n_ref = os.path.join(subdir, n_ref_name)
    h = fits.open(n_cur)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.001)
    assert result.identical, result.report()

    # Compare psfsub product
    n_cur = 'jw9999947001_02102_00001_nrcb3_a3001_psfsub.fits'
    n_ref_name = 'jw9999947001_02102_00001_nrcb3_psfsub_ref.fits'
    n_ref = os.path.join(subdir, n_ref_name)
    h = fits.open(n_cur)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.001)
    assert result.identical, result.report()

    # Compare level-2c product
    n_cur = 'jw9999947001_02102_00001_nrcb3_a3001_crfints.fits'
    n_ref_name = 'jw9999947001_02102_00001_nrcb3_a3001_crfints_ref.fits'
    n_ref = os.path.join(subdir, n_ref_name)
    h = fits.open(n_cur)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.001)
    assert result.identical, result.report()

    n_cur = 'jw9999947001_02102_00002_nrcb3_a3001_crfints.fits'
    n_ref = os.path.join(subdir,
                        'jw9999947001_02102_00002_nrcb3_a3001_crfints_ref.fits')
    h = fits.open(n_cur)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.001)
    assert result.identical, result.report()

    # Compare i2d product
    n_cur = 'jw99999-a3001_t1_nircam_f140m-maskbar_i2d.fits'
    n_ref = os.path.join(subdir,
                         'jw99999-a3001_t1_nircam_f140m-maskbar_i2d_ref.fits')
    h = fits.open(n_cur)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF', 'HDRTAB'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.001)
    assert result.identical, result.report()

    # Compare the HDRTAB in the i2d product
    result = fits.diff.HDUDiff(h['HDRTAB'],
                              href['HDRTAB'],
                              ignore_keywords=ignore_keywords+['NAXIS1', 'TFORM*'],
                              ignore_fields=ignore_keywords,
                              rtol=0.001)
    assert result.identical, result.report()
