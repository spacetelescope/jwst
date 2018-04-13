import os
import pytest
from jwst.pipeline.calwebb_tso3 import Tso3Pipeline
# DISABLED
# import pandokia.helpers.filecomp as filecomp
from .test_nrc_tso3_1 import perform_FITS_comparison

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_tso3_pipeline3(_bigdata):
    """Regression test of calwebb_tso3 on NIRISS SOSS simulated data."""
    testname = "test_tso3_pipeline2"

    testname = 'test_tso3_pipeline3'
    # Define where this test data resides in testing tree
    subdir3 = 'pipelines/niriss_caltso3/'  # Can NOT start with path separator

    # You need to have a tda dict for:
    #  - recording information to make FlagOK work
    #  - recording parameters to the task as attributes
    global tda
    tda = {}

    output = [
        # one dict for each output file to compare (i.e. each <val>)
        {
            'file'      : 'jw87600-a3001_t1_niriss_clear-gr700xd_whtlt.ecsv',
            'reference' : os.path.join(_bigdata,subdir3,'jw87600-a3001_t1_niriss_clear-gr700xd_whtlt_ref.ecsv'),
            'comparator': 'diff',
            'args'      : {},
        }
    ]

    try:
        os.remove("jw87600024001_02101_00001_nis_a3001_crfints.fits")
        os.remove("jw87600-a3001_t1_niriss_clear-gr700xd_x1dints.fits")
        os.remove("jw87600-a3001_t1_niriss_clear-gr700xd_whtlt.ecsv")
    except Exception:
        pass

    # Run pipeline step...
    asn_file = os.path.join(_bigdata, subdir3,
                            "jw87600-a3001_20170527t111213_tso3_001_asn.json")
    Tso3Pipeline.call(asn_file,
                      config_file='calwebb_tso3_1.cfg')

    # Compare level-2c product
    fname = 'jw87600024001_02101_00001_nis_a3001_crfints.fits'
    reffile = 'jw87600-a3001_t1_niriss_clear-gr700xd_crfints_ref.fits'
    extn_list = ['primary', 'sci', 'dq', 'err']

    refname = os.path.join(_bigdata, subdir3, reffile)

    result = perform_FITS_comparison(fname, refname, extn_list=extn_list)
    assert result.identical, result.report()

    # Compare level-3 product
    fname = 'jw87600-a3001_t1_niriss_clear-gr700xd_x1dints.fits'
    reffile = 'jw87600-a3001_t1_niriss_clear-gr700xd_x1dints_ref.fits'
    extn_list = ['primary', 'extract1d']

    refname = os.path.join(_bigdata, subdir3, reffile)

    result = perform_FITS_comparison(fname, refname, extn_list=extn_list)
    assert result.identical, result.report()

    # compare the output files - use this exact command
    # DISABLED
    # filecomp.compare_files(output, (__file__, testname), tda=tda,)
