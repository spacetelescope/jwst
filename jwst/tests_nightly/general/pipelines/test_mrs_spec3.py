import os
import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_spec3_pipeline1(_bigdata):
    """
    Regression test of calwebb_spec3 pipeline on simulated
    MIRI MRS dithered data.
    """

    subdir = os.path.join(_bigdata, 'pipelines', 'mrs_calspec3')
    asn_file = os.path.join(subdir, "test_asn4.json")

    step = Spec3Pipeline()
    step.save_bsub = False
    step.mrs_imatch.suffix = 'mrs_imatch'
    step.mrs_imatch.bkg_degree = 1
    step.mrs_imatch.subtract = False
    step.outlier_detection.skip = True
    step.output_use_model = True
    step.resample_spec.save_results = True
    step.resample_spec.suffix = 's2d'
    step.cube_build.save_results = True
    step.cube_build.suffix = 's3d'
    step.extract_1d.save_results = True
    step.extract_1d.suffix = 'x1d'
    step.run(asn_file)

    ignore_kws = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX']

    # Compare cube product 1
    n_cur = 'det_image_ch1-short_s3d.fits'
    n_ref = os.path.join(subdir, 'det_image_ch1-short_s3d_ref.fits')


    h = fits.open(n_cur)
    href = fits.open(n_ref)
    newh = fits.HDUList([h['primary'], h['sci'], h['err'], h['dq'], h['wmap']])
    newhref = fits.HDUList([href['primary'], href['sci'], href['err'],
                          href['dq'], href['wmap']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_kws,
                              rtol=0.000001)
    assert result.identical, result.report()

    # Compare cube product 2
    n_cur = 'det_image_ch2-short_s3d.fits'
    n_ref = os.path.join(subdir, 'det_image_ch2-short_s3d_ref.fits')


    h = fits.open(n_cur)
    href = fits.open(n_ref)
    newh = fits.HDUList([h['primary'], h['sci'], h['err'], h['dq'], h['wmap']])
    newhref = fits.HDUList([href['primary'], href['sci'], href['err'],
                          href['dq'], href['wmap']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_kws,
                              rtol=0.000001)
    assert result.identical, result.report()

    # Compare x1d product 1
    n_cur = 'det_image_ch1-short_x1d.fits'
    n_ref = os.path.join(subdir, 'det_image_ch1-short_x1d_ref.fits')


    h = fits.open(n_cur)
    href = fits.open(n_ref)
    newh = fits.HDUList([h['primary'], h['extract1d']])
    newhref = fits.HDUList([href['primary'], href['extract1d']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_kws,
                              rtol=0.000001)
    assert result.identical, result.report()

    # Compare x1d product 2
    n_cur = 'det_image_ch2-short_x1d.fits'
    n_ref = os.path.join(subdir, 'det_image_ch2-short_x1d_ref.fits')


    h = fits.open(n_cur)
    href = fits.open(n_ref)
    newh = fits.HDUList([h['primary'], h['extract1d']])
    newhref = fits.HDUList([href['primary'], href['extract1d']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_kws,
                              rtol=0.000001)
    assert result.identical, result.report()
