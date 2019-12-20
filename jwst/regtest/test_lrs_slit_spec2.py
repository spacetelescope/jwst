import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run the calwebb_spec2 pipeline on an ASN of nodded MIRI LRS
       fixedslit exposures."""

    rtdata = rtdata_module

    # Get the cfg files
    collect_pipeline_cfgs("config")

    # Get the spec2 ASN and its members
    rtdata.get_asn("miri/lrs/jw00623-o032_20191210t195246_spec2_001_asn.json")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    # NOTE: THE RESAMPLE_SPEC STEP IS SKIPPED FOR NOW, BECAUSE IT HAS A BUG
    # (the s2d image is all zeros)
    args = ["config/calwebb_spec2.cfg", rtdata.input,
            "--steps.resample_spec.skip=true",       # remove when bug fixed
            "--save_bsub=true",
            "--steps.flat_field.save_results=true",
            "--steps.srctype.save_results=true"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("output",[
    "bsub", "flat_field", "srctype", "cal", "x1d"])
def test_miri_lrs_slit_spec2(run_pipeline, fitsdiff_default_kwargs, output):
    """Regression test of the calwebb_spec2 pipeline on MIRI
       LRS fixedslit data using along-slit-nod pattern for
       background subtraction."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    rtdata.output = "jw00623032001_03102_00001_mirimage_" + output + ".fits"

    # Get the truth files
    rtdata.get_truth(os.path.join("truth", "test_miri_lrs_slit_spec2", rtdata.output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
