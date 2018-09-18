import os
import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_image3 import Image3Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_image3_pipeline2(_bigdata):
    """Regression test definitions for CALIMAGE3 pipeline.

    Regression test of calwebb_image3 pipeline on NIRCam
    simulated long-wave data with a 6-point dither.
    """

    subdir = os.path.join(_bigdata, 'pipelines', 'nircam_calimage3')

    ignore_keywords = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX']

    asn_file = os.path.join(subdir, "jw10002-o001_20171116t191235_image3_002_asn.json")

    pipe = Image3Pipeline()
    pipe.tweakreg.save_catalogs = False
    pipe.tweakreg.catalog_format = 'ecsv'
    pipe.tweakreg.kernel_fwhm = 2.
    pipe.tweakreg.snr_threshold = 5.
    pipe.tweakreg.enforce_user_order = True
    pipe.tweakreg.expand_refcat = False
    pipe.tweakreg.minobj = 15
    pipe.tweakreg.searchrad = 10.0
    pipe.tweakreg.use2dhist = True
    pipe.tweakreg.separation = 0.5
    pipe.tweakreg.tolerance = 1.0
    pipe.tweakreg.xoffset = 0.0
    pipe.tweakreg.yoffset = 0.0
    pipe.tweakreg.fitgeometry = 'rscale'
    pipe.tweakreg.nclip = 3
    pipe.tweakreg.sigma = 3.0
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
    pipe.resample.suffix = 'i2d'
    pipe.source_catalog.kernel_fwhm = 3.
    pipe.source_catalog.kernel_xsize = 5.
    pipe.source_catalog.kernel_ysize = 5.
    pipe.source_catalog.snr_threshold = 3.
    pipe.source_catalog.npixels = 50
    pipe.source_catalog.deblend = False

    pipe.run(asn_file)

    # Compare one of the level-2c crf products
    n_cur = 'jw10002001001_01101_00004_nrcblong_o001_crf.fits'
    ref_filename = 'jw10002001001_01101_00004_nrcblong_o001_crf_ref.fits'
    n_ref = os.path.join(subdir, ref_filename)
    h = fits.open(n_cur)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.00001)
    assert result.identical, result.report()

    # Compare the i2d product
    n_cur = 'jw10002-o001_t002_nircam_f444w_i2d.fits'
    n_ref = os.path.join(subdir, 'jw10002-o001_t002_nircam_f444w_i2d_ref.fits')
    h = fits.open(n_cur)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF', 'HDRTAB'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.0001)
    assert result.identical, result.report()

    # Compare the HDRTAB in the i2d product
    result = fits.diff.HDUDiff(h['HDRTAB'],
                              href['HDRTAB'],
                              ignore_keywords=ignore_keywords+['NAXIS1', 'TFORM*'],
                              ignore_fields=ignore_keywords,
                              rtol=0.0001)
    assert result.identical, result.report()
