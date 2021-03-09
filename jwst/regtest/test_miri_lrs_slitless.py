import pytest
from astropy.io.fits.diff import FITSDiff
from numpy.testing import assert_allclose

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step
from gwcs.wcstools import grid_from_bounding_box
from jwst.associations.asn_from_list import asn_from_list
from jwst import datamodels

DATASET_ID = "jw00623026001_03106_00005_mirimage"
PRODUCT_NAME = "jw00623-a3001_t001_miri_p750l-slitlessprism"
PROGRAM = "00623"


@pytest.fixture(scope="module")
def run_tso1_pipeline(jail, rtdata_module):
    """Run the calwebb_tso1 pipeline on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module
    collect_pipeline_cfgs("config")
    rtdata.get_data(f"miri/lrs/{DATASET_ID}_uncal.fits")

    args = [
        "config/calwebb_tso1.cfg",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.jump.rejection_threshold=150",
        "--steps.jump.save_results=True"
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso_spec2_pipeline(run_tso1_pipeline, jail, rtdata_module):
    """Run the calwebb_tso-spec2 pipeline on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module

    rtdata.input = f"{DATASET_ID}_rateints.fits"

    collect_pipeline_cfgs("config")

    args = [
        "config/calwebb_tso-spec2.cfg",
        rtdata.input,
        "--steps.flat_field.save_results=true",
        "--steps.assign_wcs.save_results=true",
        "--steps.srctype.save_results=true",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def generate_tso3_asn():
    """Generate an association file that references the output of run_spec2_pipeline."""
    asn = asn_from_list([f"{DATASET_ID}_calints.fits"], product_name=PRODUCT_NAME)
    asn.data["program"] = PROGRAM
    asn.data["asn_type"] = "tso3"
    asn.sequence = 1

    name, serialized = asn.dump(format="json")
    with open(name, "w") as f:
        f.write(serialized)

    return name, asn["asn_id"]


@pytest.fixture(scope="module")
def run_tso3_pipeline(run_tso_spec2_pipeline, generate_tso3_asn):
    """Run the calwebb_tso3 pipeline on the output of run_spec2_pipeline."""
    asn_filename, _ = generate_tso3_asn

    args = [
        "config/calwebb_tso3.cfg",
        asn_filename,
        "--steps.outlier_detection.save_results=true",
    ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("step_suffix", ['dq_init', 'saturation', 'refpix',
    'linearity', 'dark_current', 'jump', 'rate', 'rateints'])
def test_miri_lrs_slitless_tso1(run_tso1_pipeline, rtdata_module, fitsdiff_default_kwargs, step_suffix):
    """
    Regression test of tso1 pipeline performed on MIRI LRS slitless TSO data.
    """
    rtdata = rtdata_module
    output_filename = f"{DATASET_ID}_{step_suffix}.fits"
    rtdata.output = output_filename

    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso1/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("step_suffix", ["flat_field", "srctype", "calints", "assign_wcs", "x1dints"])
def test_miri_lrs_slitless_tso_spec2(run_tso_spec2_pipeline, rtdata_module, fitsdiff_default_kwargs,
    step_suffix):
    """Compare the output of a MIRI LRS slitless calwebb_tso-spec2 pipeline."""
    rtdata = rtdata_module

    output_filename = f"{DATASET_ID}_{step_suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso_spec2/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("step_suffix", ["outlier_detection", "crfints"])
def test_miri_lrs_slitless_tso3(run_tso3_pipeline, generate_tso3_asn, rtdata_module,
    fitsdiff_default_kwargs, step_suffix):
    """Compare the output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module
    _, asn_id = generate_tso3_asn

    output_filename = f"{DATASET_ID}_{asn_id}_{step_suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso3/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_lrs_slitless_tso3_x1dints(run_tso3_pipeline, rtdata_module,
    fitsdiff_default_kwargs):
    """Compare the output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module

    output_filename = f"{PRODUCT_NAME}_x1dints.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso3/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_lrs_slitless_tso3_whtlt(run_tso3_pipeline, generate_tso3_asn,
    rtdata_module, diff_astropy_tables):
    """Compare the whitelight output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module
    _, asn_id = generate_tso3_asn

    output_filename = f"{PRODUCT_NAME}_whtlt.ecsv"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso3/{output_filename}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)


@pytest.mark.bigdata
def test_miri_lrs_slitless_wcs(run_tso_spec2_pipeline, fitsdiff_default_kwargs,
                               rtdata_module):
    """Compare the assign_wcs output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module
    output = f"{DATASET_ID}_assign_wcs.fits"
    # get input assign_wcs and truth file
    rtdata.output = output
    rtdata.get_truth("truth/test_miri_lrs_slitless_tso_spec2/"+output)

    # Compare the output and truth file
    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        x, y = grid_from_bounding_box(im.meta.wcs.bounding_box)
        ra, dec, lam = im.meta.wcs(x, y)
        ratruth, dectruth, lamtruth = im_truth.meta.wcs(x, y)
        assert_allclose(ra, ratruth)
        assert_allclose(dec, dectruth)
        assert_allclose(lam, lamtruth)
