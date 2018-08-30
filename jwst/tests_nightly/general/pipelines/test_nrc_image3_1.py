import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_image3 import Image3Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_image3_pipeline1(_bigdata):
    """Regression test definitions for CALIMAGE3 pipeline.

    Regression test of calwebb_image3 pipeline on NIRCam
    simulated long-wave data.
    """

    subdir = os.path.join(_bigdata, 'pipelines', 'nircam_calimage3')

    ignore_keywords = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX']

    asn_file = os.path.join(subdir, "mosaic_long_asn.json")

    pipe = Image3Pipeline()
    pipe.tweakreg.skip = True
    pipe.skymethod = 'global+match'
    pipe.skymatch.match_down = True
    pipe.skymatch.subtract = False
    pipe.skymatch.skystat = 'mode'
    pipe.skymatch.nclip = 5
    pipe.skymatch.lsigma = 4.0
    pipe.skymatch.usigma = 4.0
    pipe.skymatch.binwidth = 0.1
    pipe.outlier_detection.weight_type = 'exptime'
    pipe.outlier_detection.pixfrac = 1.0
    pipe.outlier_detection.kernel = 'square'
    pipe.outlier_detection.fillval = 'INDEF'
    pipe.outlier_detection.nlow = 0
    pipe.outlier_detection.nhigh = 0
    pipe.outlier_detection.maskpt = 0.7
    pipe.outlier_detection.grow = 1
    pipe.outlier_detection.snr = '4.0 3.0'
    pipe.outlier_detection.scale = '0.5 0.4'
    pipe.outlier_detection.backg = 0.0
    pipe.outlier_detection.save_intermediate_results = False
    pipe.outlier_detection.resample_data = True
    pipe.outlier_detection.good_bits = 4
    pipe.resample.single = False
    pipe.resample.weight_type = 'exptime'
    pipe.resample.pixfrac = 1.0
    pipe.resample.kernel = 'square'
    pipe.resample.fillval = 'INDEF'
    pipe.resample.good_bits = 4
    pipe.resample.blendheaders = True
    pipe.source_catalog.kernel_fwhm = 3.
    pipe.source_catalog.kernel_xsize = 5.
    pipe.source_catalog.kernel_ysize = 5.
    pipe.source_catalog.snr_threshold = 3.
    pipe.source_catalog.npixels = 50
    pipe.source_catalog.deblend = False

    pipe.run(asn_file)

    # Compare level-2c crf product
    n_cur = 'nrca5_47Tuc_subpix_dither1_newpos_a3001_crf.fits'
    ref_filename = 'nrca5_47Tuc_subpix_dither1_newpos_cal-a3001_ref.fits'
    n_ref = os.path.join(subdir, ref_filename)
    h = pf.open(n_cur)
    href = pf.open(n_ref)
    result = pf.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.00001)
    assert result.identical, result.report()

    # Compare i2d product
    n_cur = 'mosaic_long_i2d.fits'
    n_ref = os.path.join(subdir, 'mosaic_long_i2d_ref.fits')
    h = pf.open(n_cur)
    href = pf.open(n_ref)
    result = pf.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF', 'HDRTAB'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.0001)
    assert result.identical, result.report()

    # Compare the HDRTAB in the i2d product
    result = pf.diff.HDUDiff(h['HDRTAB'],
                              href['HDRTAB'],
                              ignore_keywords=ignore_keywords+['NAXIS1', 'TFORM*'],
                              ignore_fields=ignore_keywords,
                              rtol=0.0001)
    assert result.identical, result.report()
