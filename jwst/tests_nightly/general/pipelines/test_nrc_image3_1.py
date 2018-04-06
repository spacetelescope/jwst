import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_image3 import Image3Pipeline

BIGDATA = os.environ['TEST_BIGDATA']


def test_image3_pipeline1():
    """Regression test definitions for CALIMAGE3 pipeline.

    Regression test of calwebb_image3 pipeline on NIRCam
    simulated long-wave data.

    """
    subdir = os.path.join(BIGDATA, 'pipelines', 'nircam_calimage3')
    ignore_kws = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX']

    try:
        os.remove("nrca5_47Tuc_subpix_dither1_newpos_a3001_crf.fits")
        os.remove("mosaic_long_i2d.fits")
        os.remove("mosaic_long_cat.ecsv")
    except:
        pass

    asn_file = os.path.join(subdir, "mosaic_long_asn.json")
    Image3Pipeline.call(asn_file, config_file='calwebb_image3.cfg')

    # Compare level-2c product
    n_cur = 'nrca5_47Tuc_subpix_dither1_newpos_a3001_crf.fits'
    ref_filename = 'nrca5_47Tuc_subpix_dither1_newpos_cal-a3001_ref.fits'
    n_ref = os.path.join(subdir, ref_filename)
    print(' Fitsdiff comparison between the level-2c file - a:', n_cur)
    print(' ... and the reference file - b:', n_ref)

    h = pf.open(n_cur)
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'], h['sci'], h['err'], h['dq']])
    newhref = pf.HDUList([href['primary'], href['sci'],
                          href['err'], href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              rtol=0.00001)

    result.report()
    try:
        assert result.identical is True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # Compare resampled product
    n_cur = 'mosaic_long_i2d.fits'
    n_ref = os.path.join(subdir, 'mosaic_long_i2d_ref.fits')

    print(' Fitsdiff comparison between the resampled file - a:', n_cur)
    print(' ... and the reference file - b:', n_ref)

    h = pf.open(n_cur)
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'], h['sci'], h['con'],
                       h['wht'], h['hdrtab']])
    newhref = pf.HDUList([href['primary'], href['sci'],
                         href['con'], href['wht'], href['hdrtab']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_kws,
                              rtol=0.00001)
    result.report()
    try:
        assert result.identical is True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)
