import os

import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.lib.suffix import replace_suffix
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run the calwebb_spec2 pipeline on a NIRSpec lamp exposure."""
    rtdata = rtdata_module

    # Get the input exposure
    rtdata.get_data('nirspec/lamp/jw00626030001_02103_00005_nrs1_rate.fits')

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = ["jwst.pipeline.Spec2Pipeline", rtdata.input,
            "--steps.assign_wcs.save_results=true",
            "--steps.msa_flagging.save_results=true",
            "--steps.extract_2d.save_results=true",
            "--steps.flat_field.save_results=true"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", [
    "assign_wcs", "msa_flagging", "flat_field", "cal"])
def test_nirspec_lamp_ifu_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline on a
       NIRSpec lamp exposure in IFU mode."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = replace_suffix(
            os.path.splitext(os.path.basename(rtdata.input))[0], suffix) + '.fits'
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_lamp_spec2", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
