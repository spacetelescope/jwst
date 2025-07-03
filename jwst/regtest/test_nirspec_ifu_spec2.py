"""Regression tests for NIRSpec IFU"""

import warnings

import pytest

from jwst.regtest import regtestdata as rt
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

# Define artifactory source and truth
INPUT_PATH = "nirspec/ifu"
TRUTH_PATH = "truth/test_nirspec_ifu"

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_spec2(rtdata_module, resource_tracker):
    """Run the Spec2Pipeline on a spec2 ASN containing a single exposure"""
    rtdata = rtdata_module

    # Setup the inputs
    asn_name = "jw01251-o004_20221105t023552_spec2_031_asn.json"
    asn_path = INPUT_PATH + "/" + asn_name

    # Run the pipeline
    step_params = {
        "input_path": asn_path,
        "step": "calwebb_spec2",
        "args": [
            "--steps.assign_wcs.save_results=true",
            "--steps.msa_flagging.save_results=true",
            "--steps.srctype.save_results=true",
            "--steps.flat_field.save_results=true",
            "--steps.pathloss.save_results=true",
        ],
    }

    with resource_tracker.track():
        rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


@pytest.fixture(scope="module")
def run_spec2_nsclean(rtdata_module, resource_tracker):
    """Run the Spec2Pipeline on a spec2 ASN containing a single exposure"""
    rtdata = rtdata_module

    # Set up the inputs
    asn_name = "jw01251-o004_20221105t023552_spec2_031_asn.json"
    asn_path = INPUT_PATH + "/" + asn_name

    # Run the pipeline
    step_params = {
        "input_path": asn_path,
        "step": "calwebb_spec2",
        "args": [
            "--output_file=jw01251004001_03107_00002_nrs1_nsc",
            "--steps.nsclean.skip=False",
            "--steps.nsclean.save_results=true",
        ],
    }

    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


@pytest.mark.slow
def test_log_tracked_resources_spec2(log_tracked_resources, run_spec2):
    log_tracked_resources()


@pytest.mark.slow
@pytest.mark.parametrize(
    "suffix",
    ["assign_wcs", "cal", "flat_field", "msa_flagging", "pathloss", "s3d", "srctype", "x1d"],
)
def test_spec2(run_spec2, fitsdiff_default_kwargs, suffix):
    """Regression test matching output files"""
    rt.is_like_truth(run_spec2, fitsdiff_default_kwargs, suffix, truth_path=TRUTH_PATH)


@pytest.mark.slow
@pytest.mark.parametrize("suffix", ["nsclean", "cal", "s3d", "x1d"])
def test_spec2_nsclean(run_spec2_nsclean, fitsdiff_default_kwargs, suffix):
    """Regression test matching output files"""

    # Run the pipeline and retrieve outputs
    rtdata = run_spec2_nsclean
    output = f"jw01251004001_03107_00002_nrs1_nsc_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(f"{TRUTH_PATH}/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.fixture
def run_photom(rtdata):
    """Run the photom step on an NRS IFU exposure with SRCTYPE=POINT"""

    # Setup the inputs
    rate_name = "jw01251004001_03107_00002_nrs1_pathloss.fits"
    rate_path = INPUT_PATH + "/" + rate_name

    # Run the step
    step_params = {
        "input_path": rate_path,
        "step": "jwst.photom.PhotomStep",
        "args": [
            "--save_results=True",
        ],
    }

    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


def test_photom(run_photom, fitsdiff_default_kwargs):
    """Test the photom step on an NRS IFU exposure with a point source"""
    rt.is_like_truth(run_photom, fitsdiff_default_kwargs, "photomstep", truth_path=TRUTH_PATH)


@pytest.fixture
def run_extract1d(rtdata):
    """Run the extract_1d step on an NRS IFU cube with SRCTYPE=POINT"""

    # Setup the inputs
    cube_name = "jw01251004001_03107_00002_nrs1_s3d.fits"
    cube_path = INPUT_PATH + "/" + cube_name

    # Run the step
    step_params = {
        "input_path": cube_path,
        "step": "jwst.extract_1d.Extract1dStep",
        "args": [
            "--save_results=True",
        ],
    }

    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


def test_extract1d(run_extract1d, fitsdiff_default_kwargs):
    """Test the extract_1d step on an NRS IFU cube with a point source"""
    rt.is_like_truth(run_extract1d, fitsdiff_default_kwargs, "extract1dstep", truth_path=TRUTH_PATH)
