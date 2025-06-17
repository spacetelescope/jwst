import warnings

import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """
    Run the calwebb_spec3 pipeline on a NIRSpec IFU moving target.
    """
    rtdata = rtdata_module

    # Get the ASN file and input exposures
    # Note: ASN contains a background spectrum
    rtdata.get_asn("nirspec/ifu/jw01252-o001_spec3_00003_asn_with_bg.json")

    # Run the calwebb_spec3 pipeline; save results from intermediate steps
    args = ["calwebb_spec3", rtdata.input]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["cal", "o001_crf", "s3d", "x1d"])
def test_nirspec_ifu_spec3_moving_target(
    run_pipeline, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Test spec3 pipeline on a NIRSpec IFU moving target."""
    rtdata = rtdata_module

    if suffix in {"cal", "o001_crf"}:
        output = f"jw01252001001_03101_00001_nrs1_{suffix}.fits"
    else:
        output = f"jw01252-o001_t004_nirspec_prism-clear_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_ifu_spec3_moving_target/{output}")

    # Adjust tolerance for machine precision for resampled data
    if suffix == "s3d":
        fitsdiff_default_kwargs["rtol"] = 1e-2
        fitsdiff_default_kwargs["atol"] = 2e-4

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
