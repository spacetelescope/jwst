"""Regression test of MIRI conronagraphic data through the image2 pipeline,
including multiple background exposures that have a mixture of NINTS values"""

import pytest
from jwst.regtest import regtestdata as rt

from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_image2(rtdata_module, resource_tracker):
    """Run the calwebb_image2 pipeline"""

    rtdata = rtdata_module
    rtdata.get_asn("miri/coron/jw03254-c1009_20240118t185403_image2_00001_asn.json")

    args = [
        "calwebb_image2",
        rtdata.input,
        "--steps.bkg_subtract.save_results=true",
        "--steps.bkg_subtract.save_combined_background=true",
        "--steps.assign_wcs.save_results=true",
        "--steps.flat_field.save_results=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_image2(log_tracked_resources, run_image2):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix", ["bsubints", "combinedbackground", "assign_wcs", "flat_field", "calints"]
)
def test_miri_coron_image2(run_image2, fitsdiff_default_kwargs, suffix):
    """Regression test for image2 processing of MIRI coronagraphic data with background exposures"""

    rt.is_like_truth(run_image2, fitsdiff_default_kwargs, suffix, "truth/test_miri_coron_image2")
