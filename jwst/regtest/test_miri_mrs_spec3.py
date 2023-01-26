"""Regression tests for MIRI MRS modes"""
import os
import pytest

from jwst.regtest import regtestdata as rt

from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step
from jwst.associations.asn_from_list import asn_from_list

# Define artifactory source and truth
INPUT_PATH = 'miri/mrs'
TRUTH_PATH = 'truth/test_miri_mrs'


@pytest.fixture(scope='module')
def run_spec3_single(jail, rtdata_module):
    """Run the Spec3Pipeline on single cal from Spec2Pipeline"""
    
    # Running Spec2 on MIRI data creates a single cube with 2 channels
    # Running Spec3 on MIRI data creates 2 channel cubes
    # This test will perform a regression tests on the 2 channel cubes


    rtdata = rtdata_module

    # Note that we use the truth file from spec2 processing as the input to spec3
    rtdata.get_data(TRUTH_PATH + '/' + 'jw01024001001_04101_00001_mirifulong_cal.fits')

    asn = asn_from_list([rtdata.input], product_name='jw01024001001_04101_00001_mirifulong')
    asn.data["program"] = "01024"
    asn.data["asn_type"] = "spec3"
    asn.sequence = 1
    asn_name, serialized = asn.dump(format="json")
    with open(asn_name, "w") as f:
        f.write(serialized)

    rtdata.input = asn_name

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.master_background.save_results=true',
        '--steps.outlier_detection.save_results=true',
        '--steps.resample_spec.save_results=true',
        '--steps.cube_build.save_results=true',
        '--steps.extract_1d.save_results=true',
        '--steps.combine_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope='module')
def run_spec3_multi(jail, rtdata_module):
    """Run the Spec3Pipeline dithered flight data """
    rtdata = rtdata_module

    # Use an association that has the full 12 band so we can create channel cubes
    # from the two detectors. We will run tests on a cube from IFUSHORT and one from
    # IFULONG
    
    rtdata.get_asn('miri/mrs/jw01024-c10000_spec3_00001_asn.json')

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.master_background.save_results=true',
        '--steps.outlier_detection.save_results=true',
        '--steps.cube_build.save_results=true',
        '--steps.extract_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata



@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw01024-c1000_t002_miri_ch1-short_s3d.fits',
        'jw01024-c1000_t002_miri_ch1-short_x1d.fits',
        'jw01024-c1000_t002_miri_ch3-short_s3d.fits',
        'jw01024-c1000_t002_miri_ch3-short_x1d.fits'
    ],
    #ids=["ch1-short-s3d", "ch1-short-x1d", "ch3-", "ch2-x1d"]
)


def test_spec3(run_spec3_single, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""

    rtdata = run_spec3_single
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, rtdata.output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw01024001001_04101_00001_mirifushort_c1000_crf.fits',
        'jw01024001001_04101_00002_mirifushort_c1000_crf.fits',
        'jw01024001001_04101_00001_mirifulong_c1000_crf.fits',
        'jw01024001001_04101_00002_mirifulong_c1000_crf.fits',
        'jw01024-c1000_miri_ch1-shortmedium_s3d.fits',
        'jw01024-c1000_miri_ch2-shortmedium_s3d.fits',
        'jw01024-c1000_miri_ch3-shortmedium_s3d.fits',
        'jw01024-c1000_miri_ch4-shortmedium_s3d.fits',		
        'jw01024-c1000_miri_ch1-shortmedium_x1d.fits',
        'jw01024-c1000_miri_ch2-shortmedium_x1d.fits',	
        'jw01024-c1000_miri_ch3-shortmedium_x1d.fits',
        'jw01024-c1000_miri_ch4-shortmedium_x1d.fits'
    ],
    #ids=["ch1-mrs_imatch", "ch2-mrs_imatch", "ch1-crf", "ch2-crf",
    #     "ch1-s3d", "ch2-s3d", "ch1-x1d", "ch2-x1d"]
)


def test_spec3_multi(run_spec3_multi, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rt.is_like_truth(
        run_spec3_multi, fitsdiff_default_kwargs, output,
        truth_path=TRUTH_PATH,
        is_suffix=False
    )
