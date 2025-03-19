import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_detector1pipeline(rtdata_module):
    """Run calwebb_detector1 on NIRCam imaging long data"""
    rtdata = rtdata_module
    rtdata.get_data("nircam/dhs/sub164stripe4_dhs_mock_dark.fits")
    rtdata.get_data("nircam/dhs/jw04453010001_02106_00001_nrca1_genheader_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files
    args = ["calwebb_detector1", rtdata.input,
            "--steps.dq_init.save_results=True",
            "--steps.saturation.save_results=True",
            "--steps.superbias.save_results=True",
            "--steps.refpix.save_results=True",
            "--steps.linearity.save_results=True",
            "--steps.dark_current.save_results=True",
            "--steps.dark_current.override_dark=sub164stripe4_dhs_mock_dark.fits",
            "--steps.jump.save_results=True",
            "--steps.jump.rejection_threshold=50.0",
            ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["dq_init", "saturation", "superbias",
                                    "refpix", "linearity",
                                    "dark_current", "jump", "rate", "rateints"])
def test_nircam_dhs_detector1(run_detector1pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 pipeline performed on NIRCam DHS mock data."""
    rtdata = rtdata_module
    rtdata.input = "jw04453010001_02106_00001_nrca1_genheader_uncal.fits"
    output = "jw04453010001_02106_00001_nrca1_genheader_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_dhs/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
