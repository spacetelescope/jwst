import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run the calwebb_detector pipeline on a NIRSpec MOS exposure with picture frame correction."""

    rtdata = rtdata_module

    # Get the MSA metadata file referenced in the input exposure
    rtdata.get_data("nirspec/mos/jw04291001001_01_msa.fits")

    # Get an override pictureframe file
    rtdata.get_data("nirspec/mos/jwst_nirspec_pictureframe_nrs2.fits")

    # Get the input uncal file
    rtdata.get_data("nirspec/mos/jw04291001001_03101_00001_nrs2_uncal.fits")

    # Run the calwebb_detector pipeline
    args = [
        "calwebb_detector1",
        rtdata.input,
        "--steps.picture_frame.skip=false",
        "--steps.picture_frame.override_pictureframe=jwst_nirspec_pictureframe_nrs2.fits",
        "--steps.picture_frame.save_results=true",
        "--steps.picture_frame.save_mask=true",
        "--steps.picture_frame.save_correction=true",
    ]

    Step.from_cmdline(args)

    return rtdata


@pytest.mark.parametrize(
    "suffix",
    [
        "picture_frame",
        "pctfrm_mask",
        "pctfrm_correction",
        "rate",
    ],
)
def test_nirspec_picture_frame(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of the picture_frame step in the detector1 pipeline."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = f"jw04291001001_03101_00001_nrs2_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nirspec_picture_frame/" + output)

    # Ignore the custom picture frame file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_PCTFRM")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
