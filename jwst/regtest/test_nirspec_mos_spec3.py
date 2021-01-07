import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_spec3 on NIRSpec MOS data."""
    rtdata = rtdata_module
    rtdata.get_asn("nirspec/mos/jw00626-o030_20191210t193826_spec3_001_asn.json")

    # Run the calwebb_spec3 pipeline on the association
    collect_pipeline_cfgs("config")
    args = ["config/calwebb_spec3.cfg", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
@pytest.mark.parametrize("source_id", ["s00000", "s00227", "s00279", "s00443",
                                       "s00482", "s02315"])
def test_nirspec_mos_spec3(run_pipeline, suffix, source_id, fitsdiff_default_kwargs):
    """Check results of calwebb_spec3"""
    rtdata = run_pipeline

    output = "jw00626-o030_" + source_id + "_nirspec_f170lp-g235m_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nirspec_mos_spec3/" + output)

    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
