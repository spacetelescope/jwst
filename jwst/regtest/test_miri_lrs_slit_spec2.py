import os
import pytest
from astropy.io.fits.diff import FITSDiff

from numpy.testing import assert_allclose
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step
from jwst import datamodels
from gwcs.wcstools import grid_from_bounding_box

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
            "--steps.assign_wcs.save_results=true",
            "--steps.flat_field.save_results=true",
            "--steps.srctype.save_results=true"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("output",[
    "bsub", "flat_field", "assign_wcs", "srctype", "cal", "x1d"])
def test_miri_lrs_slit_spec2(run_pipeline, fitsdiff_default_kwargs, output):
    """Regression test of the calwebb_spec2 pipeline on MIRI
       LRS fixedslit data using along-slit-nod pattern for
       background subtraction."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    rtdata.output = "jw00623032001_03102_00001_mirimage_" + output + ".fits"

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_miri_lrs_slit_spec2",
                                  "jw00623032001_03102_00001_mirimage_" + output + ".fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_lrs_slit_wcs(run_pipeline, fitsdiff_default_kwargs):
    rtdata = run_pipeline

    # get input assign_wcs and truth file
    output = "jw00623032001_03102_00001_mirimage_assign_wcs.fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_miri_lrs_slit_spec2/" + output)

    # Compare the output and truth file
    with datamodels.open(output) as im, datamodels.open(rtdata.truth) as im_truth:
        x, y = grid_from_bounding_box(im.meta.wcs.bounding_box)
        ra, dec, lam = im.meta.wcs(x, y)
        ratruth, dectruth, lamtruth = im_truth.meta.wcs(x, y)
        assert_allclose(ra, ratruth)
        assert_allclose(dec, dectruth)
        assert_allclose(lam, lamtruth)

        # Test the inverse transform
        xtest, ytest = im.meta.wcs.backward_transform(ra, dec, lam)
        xtruth, ytruth = im_truth.meta.wcs.backward_transform(ratruth, dectruth, lamtruth)
        assert_allclose(xtest, xtruth)
        assert_allclose(ytest, ytruth)
