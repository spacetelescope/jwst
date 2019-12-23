import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.lib.suffix import replace_suffix
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

file_roots = ['jw00023001001_01101_00001_nrs1_', 'jw93045010001_02101_00001_nrs2_', 'jwtest1013001_01101_00001_nrs1_']
@pytest.fixture(scope="module", params=file_roots, ids=file_roots)
def run_pipeline(jail, rtdata_module, request):
    """Run the calwebb_spec2 pipeline on NIRSpec Fixed-Slit exposures.
       We currently test the following types of inputs:
         1) Full-frame exposure (all slits will be extracted)
         2) ALLSLITS subarray exposure (all slits will be extracted)
         3) S400A1 subarray exposure (1 slit extracted)"""

    rtdata = rtdata_module

    # Get the cfg files
    collect_pipeline_cfgs("config")

    # Get the input exposure
    rtdata.get_data('nirspec/spectroscopic/' + request.param + 'rate.fits')

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = ["config/calwebb_spec2.cfg", rtdata.input,
            "--steps.assign_wcs.save_results=true",
            "--steps.extract_2d.save_results=true",
            "--steps.flat_field.save_results=true",
            "--steps.srctype.save_results=true",
            "--steps.pathloss.save_results=true"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix",[
    "assign_wcs", "extract_2d", "flat_field", "pathloss", "srctype",
    "cal", "s2d", "x1d"])
def test_nirspec_fs_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline on a
       NIRSpec FS exposures."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = replace_suffix(
            os.path.splitext(os.path.basename(rtdata.input))[0], suffix) + '.fits'
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_fs_spec2", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
