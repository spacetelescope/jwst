import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module, resource_tracker):
    """Run the calwebb_spec2 pipeline on a single NIRSpec MOS exposure."""

    rtdata = rtdata_module

    # Get the MSA metadata file referenced in the input exposure
    rtdata.get_data("nirspec/mos/jw01345066001_01_msa.fits")

    # Get the input ASN file and exposures
    rtdata.get_asn("nirspec/mos/jw01345-o066_20230831t181155_spec2_00010_asn.json")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.msa_flagging.save_results=true",
        "--steps.bkg_subtract.save_results=true",
        "--steps.extract_2d.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.wavecorr.save_results=true",
        "--steps.flat_field.save_results=true",
        "--steps.pathloss.save_results=true",
        "--steps.barshadow.save_results=true",
    ]

    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


@pytest.fixture(scope="module")
def run_pipeline_nsclean(rtdata_module, resource_tracker):
    """Run the calwebb_spec2 pipeline on NIRSpec MOS with nsclean."""

    rtdata = rtdata_module

    # Get the MSA metadata file referenced in the input exposure
    rtdata.get_data("nirspec/mos/jw01345066001_01_msa.fits")

    # Get the input ASN file and exposures
    rtdata.get_asn("nirspec/mos/jw01345-o066_20230831t181155_spec2_00010_asn.json")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--output_file=jw01345066001_05101_00003_nrs1_nsc",
        "--steps.nsclean.skip=False",
        "--steps.nsclean.save_results=True",
    ]

    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources(log_tracked_resources, run_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix",
    [
        "assign_wcs",
        "msa_flagging",
        "bsub",
        "extract_2d",
        "wavecorr",
        "flat_field",
        "srctype",
        "pathloss",
        "barshadow",
        "cal",
        "s2d",
        "x1d",
    ],
)
def test_nirspec_mos_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline on a
    NIRSpec MOS exposure. Using an exposure that's part of a
    3-shutter nod sequence, so there are nodded exposures available
    with which to do background subtraction."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = f"jw01345066001_05101_00003_nrs1_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nirspec_mos_spec2/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("suffix", ["nsclean", "cal", "s2d", "x1d"])
def test_nirspec_mos_spec2_nsclean(run_pipeline_nsclean, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline with nsclean."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline_nsclean
    output = f"jw01345066001_05101_00003_nrs1_nsc_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nirspec_mos_spec2/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
