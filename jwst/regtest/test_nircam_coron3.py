import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_coron3 on coronographic data."""
    rtdata = rtdata_module
    psfmask = rtdata.get_data("nircam/coron/jwst_nircam_psfmask_somb.fits")
    rtdata.get_asn("nircam/coron/jw99999-a3001_20170327t121212_coron3_001_asn.json")

    # Run the calwebb_coron3 pipeline on the association
    collect_pipeline_cfgs("config")
    args = [
        "config/calwebb_coron3.cfg", rtdata.input,
        f"--steps.align_refs.override_psfmask={psfmask}"
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["psfalign", "psfsub", "crfints"])
@pytest.mark.parametrize("exposure", ["00001", "00002"])
def test_nircam_coron3_sci_exp(run_pipeline, suffix, exposure, fitsdiff_default_kwargs):
    """Check intermediate results of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw9999947001_02102_" + exposure + "_nrcb3_a3001_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nircam_coron3/" + output)

    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["crfints"])
@pytest.mark.parametrize("exposure", ["00003", "00004", "00005"])
def test_nircam_coron3_psf_exp(run_pipeline, suffix, exposure, fitsdiff_default_kwargs):
    """Check intermediate results of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw9999947001_02102_" + exposure + "_nrcb3_a3001_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nircam_coron3/" + output)

    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["psfstack", "i2d"])
def test_nircam_coron3_product(run_pipeline, suffix, fitsdiff_default_kwargs):
    """Check final products of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw99999-a3001_t1_nircam_f140m-maskbar_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nircam_coron3/" + output)

    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
