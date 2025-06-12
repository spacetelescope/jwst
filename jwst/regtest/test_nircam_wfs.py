import os
from pathlib import Path
import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipelines(rtdata_module):
    """Run the calwebb_wfs-image2 and calwebb_wfs-image3 pipelines
    on NIRCam WFS&C images. The calwebb_wfs-image3 pipeline is
    run twice: once with default params and a second time with
    the do_refine option turned on."""

    rtdata = rtdata_module

    # Get the rate files and run the image2 pipeline on each
    rate_files = [
        "nircam/wavefront/jw02586173001_03104_00001_nrca3_rate.fits",
        "nircam/wavefront/jw02586173001_03104_00002_nrca3_rate.fits",
    ]
    for rate_file in rate_files:
        rtdata.get_data(rate_file)
        args = ["calwebb_image2", rtdata.input]
        Step.from_cmdline(args)

    # Get the image3 ASN file only, not the members. The members were just
    # produced by the calwebb_image2 pipeline run.
    rtdata.get_data("nircam/wavefront/jw02586-o173_20240923t165435_wfs-image3_00008_asn.json")

    # Run the calwebb_wfs-image3 pipeline with no refinement (default)
    args = ["calwebb_wfs-image3", rtdata.input]
    Step.from_cmdline(args)

    # Get the second version of the image3 ASN, which has a unique output
    # product name. Only get the ASN file, not the members.
    rtdata.get_data("nircam/wavefront/jw02586-o173_20240923t165435_wfs-image3_00008_asn_v2.json")

    # Run the calwebb_wfs-image3 pipeline with refinement
    args = ["calwebb_wfs-image3", rtdata.input, "--do_refine=True"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
def test_nicam_wfsimage_noextras(run_pipelines):
    """Ensure that i2d files are not created"""
    i2ds = list(Path(".").glob("*i2d*"))

    assert not i2ds


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "output",
    [
        "jw02586173001_03104_00001_nrca3_cal.fits",
        "jw02586173001_03104_00002_nrca3_cal.fits",
        "jw02586-o173_t035_nircam_f212n-wlm8-nrca3_wfscmb-04.fits",
        "jw02586-o173_t035_nircam_f212n-wlm8-nrca3_wfscmb-04_refine.fits",
    ],
    ids=["cal1", "cal2", "wfscmb1", "wfscmb2"],
)
def test_nircam_wfsimage(run_pipelines, fitsdiff_default_kwargs, output):
    """Regression test of the calwebb_wfs-image2 and calwebb_wfs-image3
    pipelines on a dithered pair of NIRCam images."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipelines
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nircam_wfs", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
