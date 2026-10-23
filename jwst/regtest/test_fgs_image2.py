import pytest

from jwst.regtest.regtestdata import RTFile
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

INPUT_DATA_PATH = "rtdata/fgs/image2"
INPUT_DATA = {
    "image2_rate": RTFile(
        file_name="jw01029001001_04201_00001_guider2_rate.fits",
        path=INPUT_DATA_PATH,
    )
}


@pytest.fixture(scope="module")
def run_fgs_image2(rtdata_module):
    rtdata = rtdata_module
    rtdata.get_data(INPUT_DATA["image2_rate"].full_path)

    args = [
        "calwebb_image2",
        rtdata.input,
        "--steps.flat_field.save_results=True",
        "--steps.resample.save_results=True",
    ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["flat_field", "cal", "i2d"])
def test_fgs_image2(run_fgs_image2, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test for FGS data in the image2 pipeline"""
    rtdata = rtdata_module
    output = f"jw01029001001_04201_00001_guider2_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"rtdata/truth/test_fgs_image2/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "i2d":
        fitsdiff_default_kwargs["rtol"] = 3e-3
        fitsdiff_default_kwargs["atol"] = 2e-2

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
