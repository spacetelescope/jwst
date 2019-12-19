import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.mark.bigdata
def test_miri_dark_cal( _jail, rtdata, fitsdiff_default_kwargs ):
    rtdata.get_data("miri/image/jw00001001001_01101_00001_MIRIMAGE_uncal.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_dark.cfg", rtdata.input,
            "--steps.rscd.save_results=True"]
    Step.from_cmdline( args )
    rtdata.output = "jw00001001001_01101_00001_MIRIMAGE_dark.fits"

    rtdata.get_truth("truth/test_miri_dark_pipeline/\
        jw00001001001_01101_00001_MIRIMAGE_dark.fits")

    diff = FITSDiff( rtdata.output, rtdata.truth, **fitsdiff_default_kwargs )
    assert diff.identical, diff.report
