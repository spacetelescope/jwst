import pytest

from jwst.regtest.regtestdata import RTData
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

INPUT_DATA_PATH = "fgs/image3"
ASN_FILE = "jw01029-o001_20240716t172128_image3_00001_asn.json"
INPUT_DATA = {
    ASN_FILE: RTData(
        file_name=ASN_FILE,
        path=INPUT_DATA_PATH,
        from_mast=True,
        asn=[
            "jw01029001001_06201_00001_guider2_cal.fits",
            "jw01029001001_06201_00002_guider2_cal.fits",
            "jw01029001001_06201_00003_guider2_cal.fits",
            "jw01029001001_06201_00004_guider2_cal.fits",
        ],
    )
}


@pytest.fixture(scope="module")
def run_fgs_image3(rtdata_module):
    rtdata = rtdata_module
    rtdata.get_asn(INPUT_DATA[ASN_FILE])

    args = ["calwebb_image3", rtdata.input]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["i2d", "segm"])
def test_fgs_image3(run_fgs_image3, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test for FGS data in the image3 pipeline"""
    rtdata = rtdata_module
    output = f"jw01029-o001_t009_fgs_clear_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_fgs_image3/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_fgs_image3_catalog(run_fgs_image3, rtdata_module, diff_astropy_tables):
    rtdata = rtdata_module
    rtdata.output = "jw01029-o001_t009_fgs_clear_cat.ecsv"
    rtdata.get_truth("truth/test_fgs_image3/jw01029-o001_t009_fgs_clear_cat.ecsv")

    assert diff_astropy_tables(rtdata.output, rtdata.truth, rtol=1e-3, atol=1e-4)
