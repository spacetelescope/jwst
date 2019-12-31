import os
import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.lib.suffix import replace_suffix
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.fixture(scope="module")
def run_tso_spec2_pipeline(jail, rtdata_module, request):
    """Run the calwebb_spec2 pipeline performed on NIRSpec
        fixed-slit data that uses the NRS_BRIGHTOBJ mode (S1600A1 slit)
    """

    rtdata = rtdata_module

    # Get the input exposure
    rtdata.get_data('nirspec/tso/jw84600042001_02101_00001_nrs2_rateints.fits')

    # Run the calwebb_spec2 pipeline;
    collect_pipeline_cfgs("config")
    args = ["config/calwebb_tso-spec2.cfg", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix",['calints', 'x1dints'])
def test_nirspec_brightobj_spec2(run_tso_spec2_pipeline, fitsdiff_default_kwargs, suffix):
    """
        Regression test of calwebb_spec2 pipeline performed on NIRSpec
        fixed-slit data that uses the NRS_BRIGHTOBJ mode (S1600A1 slit).
    """
    rtdata = run_tso_spec2_pipeline
    output = replace_suffix(
            os.path.splitext(os.path.basename(rtdata.input))[0], suffix) + '.fits'
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_brightobj_spec2", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
