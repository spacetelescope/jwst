import os

import pytest
from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose
from stdatamodels.jwst import datamodels

from jwst.regtest.regtestdata import RTData, trim_tso_data
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

INPUT_DATA_PATH = "miri/lrs"
RTDATA_TESTING_PATH = "rtdata_testing"
DATASET1_ID = "jw01536028001_03103_00001-seg001_mirimage"
DATASET2_ID = "jw01536028001_03103_00001-seg002_mirimage"
# DATASET3_ID = "jw01281001001_04103_00001-seg002_trim_mirimage"
DATASET3_ID = "jw01281001001_04103_00001-seg002_mirimage_mod"
ASN3_FILENAME = "jw01536-o028_20221202t215749_tso3_00001_asn.json"
PRODUCT_NAME = "jw01536-o028_t008_miri_p750l-slitlessprism"
ASN_ID = "o028"
# DATASET4_ID = "jw04496004001_03103_00001-seg001_mirimage_truncated"
DATASET4_ID = "jw04496004001_03103_00001-seg001_mirimage_mod"
TARG_DATASET4 = "jw04496004001_03102_00001-seg001_mirimage_rate.fits"


INPUT_DATA = {
    DATASET1_ID: RTData(
        file_name=DATASET1_ID + "_uncal.fits",
        path=INPUT_DATA_PATH,
    ),
    DATASET2_ID: RTData(
        file_name=DATASET2_ID + "_calints.fits",
        path=INPUT_DATA_PATH,
    ),
    DATASET3_ID: RTData(
        file_name=DATASET3_ID + "_uncal.fits",
        path=RTDATA_TESTING_PATH,
        from_mast=True,
        mod_code="trim_tso",
    ),
    ASN3_FILENAME: RTData(
        file_name=ASN3_FILENAME,
        path=INPUT_DATA_PATH,
        from_mast=False,
        asn_files=[
            "jw01536028001_03103_00001-seg001_mirimage_calints.fits",
            "jw01536028001_03103_00001-seg002_mirimage_calints.fits",
        ],
        asn_files_from_mast=True,
        mod_code="N/A",
        comment="N/A",
    ),
    DATASET4_ID: RTData(
        file_name=DATASET4_ID + "_rateints.fits",
        path=RTDATA_TESTING_PATH,
        from_mast=True,
        mod_code="trim_tso",
    ),
    "jw04496004001_03102_00001-seg001_mirimage_rate.fits": RTData(
        file_name="jw04496004001_03102_00001-seg001_mirimage_rate.fits",
        path=INPUT_DATA_PATH,
        from_mast=True,
    ),
}

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


def trim_tso():
    input_data = {
        INPUT_DATA[DATASET3_ID].file_name: {
            "ints_to_keep": 20,
            "intstart": 133,
            "ints_offset": 14,
        },
        INPUT_DATA[DATASET4_ID].file_name: {
            "ints_to_keep": 10,
            "intstart": 1,
            "ints_offset": 0,
        },
    }

    for file, fdict in input_data.items():
        trim_tso_data(file, fdict["ints_to_keep"], fdict["intstart"], fdict["ints_offset"])


@pytest.fixture(scope="module")
def run_tso1_pipeline(rtdata_module):
    """Run the calwebb_detector1 pipeline on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module
    rtdata.get_data(INPUT_DATA[DATASET1_ID].path + "/" + INPUT_DATA[DATASET1_ID].file_name)

    args = [
        "calwebb_detector1",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.lastframe.save_results=True",
        "--steps.reset.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_detector1_pipeline(rtdata_module):
    """Run calwebb_detector pipeline on a MIRI LRS slitless exposure for Segment 2 data.
    Focusing on the steps that depend on integration # and not covered by run_tso1_pipeline.
    Also test running RSC step"""
    rtdata = rtdata_module
    rtdata.get_data(INPUT_DATA[DATASET3_ID].path + "/" + INPUT_DATA[DATASET3_ID].file_name)

    args = [
        "calwebb_detector1",
        rtdata.input,
        "--steps.emicorr.save_results=True",
        "--steps.emicorr.algorithm=sequential",
        "--steps.rscd.skip=False",
        "--steps.rscd.save_results=True",
        "--steps.dark_current.save_results=True",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_detector1_pipeline_emicorr_joint(rtdata_module):
    """Run detector1 with an alternate emicorr algorithm."""
    rtdata = rtdata_module
    rtdata.get_data(INPUT_DATA[DATASET3_ID].path + "/" + INPUT_DATA[DATASET3_ID].file_name)

    args = [
        "calwebb_detector1",
        rtdata.input,
        f"--output_file={DATASET3_ID}_emijoint",
        "--steps.emicorr.algorithm=joint",
        "--steps.emicorr.save_results=True",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso_spec2_pipeline(run_tso1_pipeline, rtdata_module, resource_tracker):
    """Run the calwebb_tso-spec2 pipeline on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module

    rtdata.input = f"{DATASET1_ID}_rateints.fits"

    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.flat_field.save_results=true",
        "--steps.pixel_replace.save_results=true",
        "--steps.pixel_replace.skip=false",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso3_pipeline(run_tso_spec2_pipeline, rtdata_module, resource_tracker):
    """Run the calwebb_tso3 pipeline on the output of run_spec2_pipeline."""
    rtdata = rtdata_module
    rtdata.get_data(INPUT_DATA[DATASET2_ID].path + "/" + INPUT_DATA[DATASET2_ID].file_name)
    rtdata.get_data(INPUT_DATA[ASN3_FILENAME].path + "/" + INPUT_DATA[ASN3_FILENAME].file_name)

    args = [
        "calwebb_tso3",
        ASN3_FILENAME,
        "--steps.outlier_detection.save_results=true",
        "--steps.outlier_detection.save_intermediate_results=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_spec2(log_tracked_resources, run_tso_spec2_pipeline):
    log_tracked_resources()


def test_log_tracked_resources_spec3(log_tracked_resources, run_tso3_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize(
    "step_suffix",
    [
        "dq_init",
        "saturation",
        "lastframe",
        "reset",
        "linearity",
        "dark_current",
        "ramp",
        "rate",
        "rateints",
    ],
)
def test_miri_lrs_slitless_tso1(
    run_tso1_pipeline, rtdata_module, fitsdiff_default_kwargs, step_suffix
):
    """Regression test of tso1 pipeline performed on MIRI LRS slitless TSO data."""
    rtdata = rtdata_module
    output_filename = f"{DATASET1_ID}_{step_suffix}.fits"
    rtdata.output = output_filename

    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso1/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "step_suffix", ["rscd", "emicorr", "dark_current", "ramp", "rate", "rateints"]
)
def test_miri_lrs_slitless_detector1(
    run_detector1_pipeline, rtdata_module, fitsdiff_default_kwargs, step_suffix
):
    """
    Regression test of  detector1 pipeline.

    Performed on MIRI LRS slitless TSO data.
    Testing segment 2 data for RSCD, emicorr and dark_current.
    """
    rtdata = rtdata_module
    output_filename = f"{DATASET3_ID}_{step_suffix}.fits"
    rtdata.output = output_filename

    rtdata.get_truth(f"truth/test_miri_lrs_slitless_detector1/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("step_suffix", ["emicorr", "rate", "rateints"])
def test_miri_lrs_slitless_detector1_emicorr_joint(
    run_detector1_pipeline_emicorr_joint, rtdata_module, fitsdiff_default_kwargs, step_suffix
):
    """Regression test of  detector1 pipeline with an alternate emicorr algorithm."""
    rtdata = rtdata_module
    output_filename = f"{DATASET3_ID}_emijoint_{step_suffix}.fits"
    rtdata.output = output_filename

    rtdata.get_truth(f"truth/test_miri_lrs_slitless_detector1_emicorr_joint/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "step_suffix", ["assign_wcs", "srctype", "flat_field", "pixel_replace", "calints", "x1dints"]
)
def test_miri_lrs_slitless_tso_spec2(
    run_tso_spec2_pipeline, rtdata_module, fitsdiff_default_kwargs, step_suffix
):
    """Compare the output of a MIRI LRS slitless calwebb_tso-spec2 pipeline."""
    rtdata = rtdata_module

    output_filename = f"{DATASET1_ID}_{step_suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso_spec2/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("step_suffix", ["outlier_detection", "crfints"])
def test_miri_lrs_slitless_tso3(
    run_tso3_pipeline, rtdata_module, fitsdiff_default_kwargs, step_suffix
):
    """Compare the output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module

    median_filename = f"{DATASET1_ID}_{ASN_ID}_median.fits"
    assert os.path.isfile(median_filename)

    output_filename = f"{DATASET1_ID}_{ASN_ID}_{step_suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso3/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_miri_lrs_slitless_tso3_x1dints(run_tso3_pipeline, rtdata_module, fitsdiff_default_kwargs):
    """Compare the output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module

    output_filename = f"{PRODUCT_NAME}_x1dints.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso3/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_miri_lrs_slitless_tso3_whtlt(run_tso3_pipeline, rtdata_module, diff_astropy_tables):
    """Compare the whitelight output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module

    output_filename = f"{PRODUCT_NAME}_whtlt.ecsv"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso3/{output_filename}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)


def test_miri_lrs_slitless_wcs(run_tso_spec2_pipeline, fitsdiff_default_kwargs, rtdata_module):
    """Compare the assign_wcs output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module
    output = f"{DATASET1_ID}_assign_wcs.fits"
    # get input assign_wcs and truth file
    rtdata.output = output
    rtdata.get_truth("truth/test_miri_lrs_slitless_tso_spec2/" + output)

    # Compare the output and truth file
    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        x, y = grid_from_bounding_box(im.meta.wcs.bounding_box)
        ra, dec, lam = im.meta.wcs(x, y)
        ratruth, dectruth, lamtruth = im_truth.meta.wcs(x, y)
        assert_allclose(ra, ratruth)
        assert_allclose(dec, dectruth)
        assert_allclose(lam, lamtruth)


@pytest.fixture(scope="module")
def run_spec2_slitless_targ_centroid(rtdata_module):
    """Run calwebb_spec2 including targ_centroid step on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module

    # science exposure
    sci = rtdata.get_data(INPUT_DATA[DATASET4_ID].path + "/" + INPUT_DATA[DATASET4_ID].file_name)

    # target acquisition verification image
    taq = rtdata.get_data(
        INPUT_DATA[TARG_DATASET4].path + "/" + INPUT_DATA[TARG_DATASET4].file_name
    )

    args = [
        "calwebb_spec2",
        sci,
        "--steps.targ_centroid.skip=false",
        "--steps.targ_centroid.save_results=true",
        f"--steps.targ_centroid.ta_file={taq}",
    ]
    Step.from_cmdline(args)


@pytest.mark.parametrize("step_suffix", ["targ_centroid", "calints", "x1dints"])
def test_miri_lrs_slitless_spec2_targ_centroid(
    run_spec2_slitless_targ_centroid, rtdata_module, fitsdiff_default_kwargs, step_suffix
):
    """Compare the output of a MIRI LRS slitless calwebb_spec2 pipeline including targ_centroid step."""
    rtdata = rtdata_module

    output_filename = f"{DATASET4_ID}_{step_suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso_spec2/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
