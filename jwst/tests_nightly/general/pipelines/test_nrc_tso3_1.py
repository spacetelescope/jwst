import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_tso3 import Tso3Pipeline
# DISABLED
# import pandokia.helpers.filecomp as filecomp

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_tso3_pipeline1(_bigdata):
    """Regression test of calwebb_tso3 pipeline on NIRCam simulated data.

    Default imaging mode outlier_detection will be tested here.
    """
    testname = "test_tso3_pipeline1"
    # Define where this test data resides in testing tree
    subdir = 'pipelines/nircam_caltso3/'  # Can NOT start with path separator

    # You need to have a tda dict for:
    #  - recording information to make FlagOK work
    #  - recording parameters to the task as attributes
    global tda
    tda = {}

    output = [
        # one dict for each output file to compare (i.e. each <val>)
        {
            'file'      : 'jw93065-a3001_t1_nircam_f150w-wlp8_phot.ecsv',
            'reference' : os.path.join(_bigdata,subdir,'jw93065-a3001_t1_nircam_f150w-wlp8_phot_ref.ecsv'),
            'comparator': 'diff',
            'args'      : {},
        }
    ]

    try:
        os.remove("jw93065002001_02101_00001_nrca1_a3001_crfints.fits")
        os.remove("jw93065-a3001_t1_nircam_f150w-wlp8_phot.ecsv")
    except Exception:
        pass

    # Run pipeline step...
    asn_file = os.path.join(_bigdata, subdir,
                            "jw93065-a3001_20170511t111213_tso3_001_asn.json")
    step = Tso3Pipeline()
    step.scale_detection = False
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

    result = perform_FITS_comparison(fname, refname, extn_list=extn_list)
    assert result.identical, result.report()

    # compare the output files - use this exact command
    # DISABLED
    # filecomp.compare_files(output, (__file__, testname), tda=tda,)


def test_tso3_pipeline2(_bigdata):
    """Regression test of calwebb_tso3 pipeline on NIRCam simulated data.

    Scaled imaging mode outlier_detection will be tested here.
    """
    testname = "test_tso3_pipeline2"

    # Define where this test data resides in testing tree
    subdir = 'pipelines/nircam_caltso3/'  # Can NOT start with path separator

    # You need to have a tda dict for:
    #  - recording information to make FlagOK work
    #  - recording parameters to the task as attributes
    global tda
    tda = {}

    output = [
        # one dict for each output file to compare (i.e. each <val>)
        {
            'file'      : 'jw93065-a3002_t1_nircam_f150w-wlp8_phot.ecsv',
            'reference' : os.path.join(_bigdata,subdir,'jw93065-a3002_t1_nircam_f150w-wlp8_phot_ref.ecsv'),
            'comparator': 'diff',
            'args'      : {},
        }
    ]

    try:
        os.remove("jw93065002002_02101_00001_nrca1_a3002_crfints.fits")
        os.remove("jw93065-a3002_t1_nircam_f150w-wlp8_phot.ecsv")
    except Exception:
        pass

    # Run pipeline step...
    asn_file = os.path.join(_bigdata, subdir,
                            "jw93065-a3002_20170511t111213_tso3_001_asn.json")
    Tso3Pipeline.call(asn_file,
                      config_file='calwebb_tso3_2.cfg')

    # Compare level-2c product
    fname = 'jw93065002002_02101_00001_nrca1_a3002_crfints.fits'
    reffile = 'jw93065002002_02101_00001_nrca1_a3002_crfints_ref.fits'
    extn_list = ['primary', 'sci', 'dq', 'err']

    refname = os.path.join(_bigdata, subdir, reffile)

    result = perform_FITS_comparison(fname, refname, extn_list=extn_list)
    assert result.identical, result.report()

    # compare the output files - use this exact command
    # DISABLED
    # filecomp.compare_files(output, (__file__, testname), tda=tda,)


# Utility function to simplify FITS comparisons
def perform_FITS_comparison(filename, refname, **pars):
    """Perform FITSDIFF comparison.

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
    h = pf.open(filename)
    href = pf.open(refname)
    newh = pf.HDUList([h[extn] for extn in extn_list])
    newhref = pf.HDUList([href[extn] for extn in extn_list])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords=ignore_kws,
                              ignore_fields=ignore_fields,
                              rtol=rtol)
    return result
