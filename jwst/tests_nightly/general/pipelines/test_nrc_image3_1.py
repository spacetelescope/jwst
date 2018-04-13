import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_image3 import Image3Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_image3_pipeline1():
    """Regression test definitions for CALIMAGE3 pipeline.

    Regression test of calwebb_image3 pipeline on NIRCam
    simulated long-wave data.

    """
    subdir = os.path.join(_bigdata, 'pipelines', 'nircam_calimage3')
    ignore_kws = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX']

    try:
        os.remove("nrca5_47Tuc_subpix_dither1_newpos_a3001_crf.fits")
        os.remove("mosaic_long_i2d.fits")
        os.remove("mosaic_long_cat.ecsv")
    except:
        pass

    asn_file = os.path.join(subdir, "mosaic_long_asn.json")
    step = Image3Pipeline()
    step.tweakreg.skip = True
    skymethod = 'global+match'
    step.skymatch.match_down =True
    step.skymatch.subtract = False
    step.skymatch.skystat = 'mode'
    step.skymatch.nclip = 5
    step.skymatch.lsigma = 4.0
    step.skymatch.usigma = 4.0
    step.skymatch.binwidth = 0.1
    step.outlier_detection.wht_type = 'exptime'
    step.outlier_detection.pixfrac = 1.0
    step.outlier_detection.kernel = 'square'
    step.outlier_detection.fillval = 'INDEF'
    step.outlier_detection.nlow = 0
    step.outlier_detection.nhigh = 0
    step.outlier_detection.maskpt = 0.7
    step.outlier_detection.grow = 1
    step.outlier_detection.snr = '4.0 3.0'
    step.outlier_detection.scale = '0.5 0.4'
    step.outlier_detection.backg = 0.0
    step.outlier_detection.save_intermediate_results = False
    step.outlier_detection.resample_data = True
    step.outlier_detection.good_bits = 4
    step.resample.single = False
    step.resample.wht_type = 'exptime'
    step.resample.pixfrac = 1.0
    step.resample.kernel = 'square'
    step.resample.fillval = 'INDEF'
    step.resample.good_bits = 4
    step.resample.blendheaders = True
    step.source_catalog.kernel_fwhm = 3.
    step.source_catalog.kernel_xsize = 5.
    step.source_catalog.kernel_ysize = 5.
    step.source_catalog.snr_threshold = 3.
    step.source_catalog.npixels = 50
    step.source_catalog.deblend = False

    step.run(asn_file)

    # Compare level-2c product
    n_cur = 'nrca5_47Tuc_subpix_dither1_newpos_a3001_crf.fits'
    ref_filename = 'nrca5_47Tuc_subpix_dither1_newpos_cal-a3001_ref.fits'
    n_ref = os.path.join(subdir, ref_filename)

    h = pf.open(n_cur)
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'], h['sci'], h['err'], h['dq']])
    newhref = pf.HDUList([href['primary'], href['sci'],
                          href['err'], href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              rtol=0.00001)
    assert result.identical, result.report()

    # Compare resampled product
    n_cur = 'mosaic_long_i2d.fits'
    n_ref = os.path.join(subdir, 'mosaic_long_i2d_ref.fits')


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
    assert result.identical, result.report()
