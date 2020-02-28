import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.fixture(scope="module")
def run_pipelines(jail, rtdata_module):
    """Run the calwebb_wfs-image2 and calwebb_wfs-image3 pipelines
       on NIRCam WFS&C images. The calwebb_wfs-image3 pipeline is
       run twice: once with default params and a second time with
       the do_refine option turned on."""

    rtdata = rtdata_module

    # Get the cfg files
    collect_pipeline_cfgs("config")

    # Get the image2 ASN and members
    rtdata.get_asn("nircam/wavefront/jw00632-o003_20191210t194815_wfs-image2_001_asn.json")

    # Run the calwebb_wfs-image2 pipeline
    args = ["config/calwebb_wfs-image2.cfg", rtdata.input]
    Step.from_cmdline(args)

    # Get the image3 ASN file only, not the members. The members were just
    # produced by the calwebb_wfs-image2 pipeline run.
    rtdata.get_data("nircam/wavefront/jw00632-o003_20191210t194815_wfs-image3_001_asn.json")

    # Run the calwebb_wfs-image3 pipeline with no refinement (default)
    args = ["config/calwebb_wfs-image3.cfg", rtdata.input]
    Step.from_cmdline(args)

    # Get the second version of the image3 ASN, which has a unique output
    # product name. Only get the ASN file, not the members.
    rtdata.get_data("nircam/wavefront/jw00632-o003_20191210t194815_wfs-image3_002_asn.json")

    # Run the calwebb_wfs-image3 pipeline with refinement
    args = ["config/calwebb_wfs-image3.cfg", rtdata.input, "--do_refine=True"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("output",[
    "jw00632003002_03105_00001_nrca4_cal.fits",
    "jw00632003002_03105_00002_nrca4_cal.fits",
    "jw00632-o003_t001_nircam_f212n-wlp8-nrca4_wfscmb-05.fits",
    "jw00632-o003_t001_nircam_f212n-wlp8-nrca4_wfscmb-05_refine.fits"],
    ids=["cal1", "cal2", "wfscmb1", "wfscmb2"])
def test_nircam_wfsimage(run_pipelines, fitsdiff_default_kwargs, output):
    """Regression test of the calwebb_wfs-image2 and calwebb_wfs-image3
       pipelines on a dithered pair of NIRCam images."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipelines
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/nircam/test_wfs-image2and3", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
