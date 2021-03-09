import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'exposure',
    ['jw00001001001_01101_00001_mirimage', 'jw02201001001_01101_00001_MIRIMAGE']
)
def test_miri_dark_pipeline(exposure, _jail, rtdata, fitsdiff_default_kwargs):
    rtdata.get_data(f"miri/image/{exposure}_uncal.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_dark.cfg", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = f"{exposure}_dark.fits"

    rtdata.get_truth(f"truth/test_miri_dark_pipeline/{exposure}_dark.fits")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
