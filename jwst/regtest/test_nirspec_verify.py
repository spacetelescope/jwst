import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

# Input is jw01970001001_02102_00001_nrs1_uncal.fits, a verify-only image
# preceding an IFU observation, with EXP_TYPE=NRS_IMAGE.
# There is currently no data with EXP_TYPE=NRS_VERIFY available in
# the archive.
ROOT = "jw01970002001_03101_00001_nrs1"


@pytest.fixture(scope="module")
def run_detector1(rtdata_module):
    """Run a verify-only NIRSpec image through detector1."""
    rtdata = rtdata_module
    rtdata.get_data(f"nirspec/imaging/{ROOT}_uncal.fits")

    args = [
        "calwebb_detector1",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.superbias.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.jump.save_results=True",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_image2(run_detector1, rtdata_module):
    """Run a verify-only NIRSpec image through image2."""
    rtdata = rtdata_module
    rtdata.input = f"{ROOT}_rate.fits"

    args = [
        "calwebb_image2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.flat_field.save_results=true",
        "--steps.photom.save_results=true",
    ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "suffix",
    [
        "dq_init",
        "saturation",
        "superbias",
        "refpix",
        "linearity",
        "dark_current",
        "jump",
        "rate",
    ],
)
def test_verify_detector1(run_detector1, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test results of the detector1 processing."""
    rtdata = rtdata_module
    output = f"{ROOT}_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_verify/{output}")

    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 1e-3

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "suffix",
    [
        "assign_wcs",
        "flat_field",
        "cal",
    ],
)
def test_verify_image2(run_image2, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test results of the image2 processing."""
    rtdata = rtdata_module
    output = f"{ROOT}_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_verify/{output}")

    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 1e-3

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
