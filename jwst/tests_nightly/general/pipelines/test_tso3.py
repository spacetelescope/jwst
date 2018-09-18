import os
import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_tso3 import Tso3Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_tso3_pipeline_nrc1(_bigdata):
    """Regression test of calwebb_tso3 pipeline on NIRCam simulated data.

    Default imaging mode outlier_detection will be tested here.
    """
    subdir = 'pipelines/nircam_caltso3/'
    asn_file = os.path.join(_bigdata, subdir,
                            "jw93065-a3001_20170511t111213_tso3_001_asn.json")
    step = Tso3Pipeline()
    step.scale_detection = False
    step.outlier_detection.weight_type = 'exptime'
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
    step.outlier_detection.resample_data = False
    step.outlier_detection.good_bits = 4
    step.extract_1d.smoothing_length = 0
    step.extract_1d.bkg_order = 0

    step.run(asn_file)

    # Compare level-2c product
    fname = 'jw93065002001_02101_00001_nrca1_a3001_crfints.fits'
    reffile = 'jw93065002001_02101_00001_nrca1_a3001_crfints_ref.fits'
    extn_list = ['primary', 'sci', 'dq', 'err']

    refname = os.path.join(_bigdata, subdir, reffile)

    result = perform_fits_comparison(fname, refname, extn_list=extn_list)
    assert result.identical, result.report()


def test_tso3_pipeline_nrc2(_bigdata):
    """Regression test of calwebb_tso3 pipeline on NIRCam simulated data.

    Scaled imaging mode outlier_detection will be tested here.
    """
    subdir = 'pipelines/nircam_caltso3/'
    asn_file = os.path.join(_bigdata, subdir,
                            "jw93065-a3002_20170511t111213_tso3_001_asn.json")
    step = Tso3Pipeline()
    step.scale_detection = True
    step.outlier_detection.weight_type = 'exptime'
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
    step.outlier_detection.resample_data = False
    step.outlier_detection.good_bits = 4
    step.extract_1d.smoothing_length = 0
    step.extract_1d.bkg_order = 0

    step.run(asn_file)

    # Compare level-2c product
    fname = 'jw93065002002_02101_00001_nrca1_a3002_crfints.fits'
    reffile = 'jw93065002002_02101_00001_nrca1_a3002_crfints_ref.fits'
    extn_list = ['primary', 'sci', 'dq', 'err']

    refname = os.path.join(_bigdata, subdir, reffile)

    result = perform_fits_comparison(fname, refname, extn_list=extn_list)
    assert result.identical, result.report()


def test_tso3_pipeline_nis(_bigdata):
    """Regression test of calwebb_tso3 on NIRISS SOSS simulated data.
    """
    subdir3 = 'pipelines/niriss_caltso3/'
    asn_file = os.path.join(_bigdata, subdir3,
                            "jw87600-a3001_20170527t111213_tso3_001_asn.json")
    Tso3Pipeline.call(asn_file)

    # Compare level-2c product
    fname = 'jw87600024001_02101_00001_nis_a3001_crfints.fits'
    reffile = 'jw87600-a3001_t1_niriss_clear-gr700xd_crfints_ref.fits'
    extn_list = ['primary', 'sci', 'dq', 'err']

    refname = os.path.join(_bigdata, subdir3, reffile)

    result = perform_fits_comparison(fname, refname, extn_list=extn_list)
    assert result.identical, result.report()

    # Compare level-3 product
    fname = 'jw87600-a3001_t1_niriss_clear-gr700xd_x1dints.fits'
    reffile = 'jw87600-a3001_t1_niriss_clear-gr700xd_x1dints_ref.fits'
    extn_list = ['primary', 'extract1d']

    refname = os.path.join(_bigdata, subdir3, reffile)

    result = perform_fits_comparison(fname, refname, extn_list=extn_list)
    assert result.identical, result.report()




# Utility function to simplify FITS comparisons
def perform_fits_comparison(filename, refname, **pars):
    """Perform FITSDiff comparison.

    Comparision will be done on `filename` in current directory with
        file that has same filename in `_bigdata` directory.

        Parameters
        ----------
        filename : str
            Filename (no path) for file to be compared

        refname : str
            Full filename (with path) of file to be used as reference

        extn_list : list
            List of FITS extensions to include in comparison.
            Default: ['primary', 'sci', 'con', 'wht', 'hdrtab']

        ignore_kws : list
            List of header keywords to ignore during comparison
            Default: ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX']

        ignore_fields : list
            List of header keywords to ignore during comparison
            Default: ['DATE']

        rtol : float
            Level of difference below which constitutes agreement in values
            Default: 0.00001

        Returns
        -------
        result : obj
            astropy.io.fits.diff.FITSDIFF object with results of comparison

        filenames : list
            List of input and reference filenames used for comparison

    """
    extn_list = pars.get('extn_list',
                         ['primary', 'sci', 'con', 'wht', 'hdrtab'])
    ignore_kws = pars.get('ignore_kws',
                          ['DATE', 'CAL_VER', 'CAL_VCS',
                           'CRDS_VER', 'CRDS_CTX'])
    ignore_fields = pars.get('ignore_fields', ['DATE'])
    rtol = pars.get('rtol', 0.00001)

    # Compare resampled product
    h = fits.open(filename)
    href = fits.open(refname)
    newh = fits.HDUList([h[extn] for extn in extn_list])
    newhref = fits.HDUList([href[extn] for extn in extn_list])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_fields,
                              rtol=rtol)
    return result
