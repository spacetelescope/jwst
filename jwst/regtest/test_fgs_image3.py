import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_fgs_image3(rtdata_module):
    rtdata = rtdata_module
    rtdata.get_asn("fgs/image3/jw01029-o001_20240716t172128_image3_00001_asn.json")

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
