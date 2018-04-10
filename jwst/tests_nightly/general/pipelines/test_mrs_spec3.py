import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline

BIGDATA = os.environ['TEST_BIGDATA']


def test_spec3_pipeline1():
    """Regression test definitions for CALSPEC3 pipeline.

    Regression test of calwebb_spec3 pipeline on simulated
    MIRI MRS dithered data.

    """
    subdir = os.path.join(BIGDATA, 'pipelines', 'mrs_calspec3')
    ignore_kws = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX']

    try:
        os.remove("det_image_ch1-short_s3d.fits")
        os.remove("det_image_ch2-short_s3d.fits")
        os.remove("det_image_ch1-short_x1d.fits")
        os.remove("det_image_ch2-short_x1d.fits")
    except:
        pass

    asn_file = os.path.join(subdir, "test_asn17.json")
    step = Spec3Pipeline()
    step.save_bsub = False
    step.output_use_model = True
    step.resample_spec.save_results = True
    step.resample_spec.skip = True
    step.extract_1d.save_results = True
    step.run(asn_file)

    # Compare cube product 1
    n_cur = 'det_image_ch1-short_s3d.fits'
    n_ref = os.path.join(subdir, 'det_image_ch1-short_s3d_ref.fits')

    print(' Fitsdiff comparison between the resampled file - a:', n_cur)
    print(' ... and the reference file - b:', n_ref)

    h = pf.open(n_cur)
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'], h['sci'], h['err'], h['dq'], h['wmap']])
    newhref = pf.HDUList([href['primary'], href['sci'], href['err'],
                          href['dq'], href['wmap']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_kws,
                              rtol=0.000001)
    result.report()
    try:
        assert result.identical is True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # Compare cube product 2
    n_cur = 'det_image_ch2-short_s3d.fits'
    n_ref = os.path.join(subdir, 'det_image_ch2-short_s3d_ref.fits')

    print(' Fitsdiff comparison between the resampled file - a:', n_cur)
    print(' ... and the reference file - b:', n_ref)

    h = pf.open(n_cur)
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'], h['sci'], h['err'], h['dq'], h['wmap']])
    newhref = pf.HDUList([href['primary'], href['sci'], href['err'],
                          href['dq'], href['wmap']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_kws,
                              rtol=0.000001)
    result.report()
    try:
        assert result.identical is True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # Compare x1d product 1
    n_cur = 'det_image_ch1-short_x1d.fits'
    n_ref = os.path.join(subdir, 'det_image_ch1-short_x1d_ref.fits')

    print(' Fitsdiff comparison between the resampled file - a:', n_cur)
    print(' ... and the reference file - b:', n_ref)

    h = pf.open(n_cur)
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'], h['extract1d']])
    newhref = pf.HDUList([href['primary'], href['extract1d']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_kws,
                              rtol=0.000001)
    result.report()
    try:
        assert result.identical is True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # Compare x1d product 2
    n_cur = 'det_image_ch2-short_x1d.fits'
    n_ref = os.path.join(subdir, 'det_image_ch2-short_x1d_ref.fits')

    print(' Fitsdiff comparison between the resampled file - a:', n_cur)
    print(' ... and the reference file - b:', n_ref)

    h = pf.open(n_cur)
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'], h['extract1d']])
    newhref = pf.HDUList([href['primary'], href['extract1d']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_kws,
                              rtol=0.000001)
    result.report()
    try:
        assert result.identical is True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)
