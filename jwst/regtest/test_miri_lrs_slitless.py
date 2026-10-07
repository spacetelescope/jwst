import os

import numpy as np
import pytest
from astropy.table import Table
from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose
from stdatamodels.jwst import datamodels

from jwst.regtest.regtestdata import RTFile, trim_tso_data
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

INPUT_DATA_PATH = "rtdata/miri/lrs"
PRODUCT_NAME = "jw01536-o028_t008_miri_p750l-slitlessprism"
ASN_ID = "o028"


INPUT_DATA = {
    "jw01536028001_03103_00001-seg001_mirimage": RTFile(
        file_name="jw01536028001_03103_00001-seg001_mirimage_uncal.fits",
        path=INPUT_DATA_PATH,
    ),
    "jw01536028001_03103_00001-seg002_mirimage": RTFile(
        file_name="jw01536028001_03103_00001-seg002_mirimage_calints.fits",
        path=INPUT_DATA_PATH,
    ),
    "jw01281001001_04103_00001-seg002_mirimage_mod": RTFile(
        file_name="jw01281001001_04103_00001-seg002_mirimage_mod_uncal.fits",
        path=INPUT_DATA_PATH,
        from_mast=True,
        mod_code="trim_tso",
    ),
    "jw01536-o028_mod_tso3_00001_asn": RTFile(
        file_name="jw01536-o028_mod_tso3_00001_asn.json",
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
    "jw04496004001_03103_00001-seg001_mirimage_mod": RTFile(
        file_name="jw04496004001_03103_00001-seg001_mirimage_mod_rateints.fits",
        path=INPUT_DATA_PATH,
        from_mast=True,
        mod_code="trim_tso",
    ),
    "jw04496004001_03102_00001-seg001_mirimage": RTFile(
        file_name="jw04496004001_03102_00001-seg001_mirimage_rate.fits",
        path=INPUT_DATA_PATH,
        from_mast=True,
    ),
}

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


def trim_tso():
    input_data = {
        INPUT_DATA["jw01281001001_04103_00001-seg002_mirimage_mod"].file_name: {
            "ints_to_keep": 20,
            "intstart": 133,
            "ints_offset": 14,
        },
        INPUT_DATA["jw04496004001_03103_00001-seg001_mirimage_mod"].file_name: {
            "ints_to_keep": 10,
            "intstart": 1,
            "ints_offset": 0,
        },
    }

    for file, fdict in input_data.items():
        trim_tso_data(file, fdict["ints_to_keep"], fdict["intstart"], fdict["ints_offset"])


def test_trim_tso_data(tmp_cwd):
    """Test trim_tso_data function."""
    # mock uncal tso data
    n_integrations = 10
    n_groups = 10
    rows = 20
    cols = 20
    data = np.ones([n_integrations, n_groups, rows, cols])
    table = Table()
    fake_mjd = np.arange(6000.0, 6010.0).tolist()
    table["integration_number"] = list(range(n_integrations))
    table["int_start_MJD_UTC"] = fake_mjd
    table["int_mid_MJD_UTC"] = fake_mjd
    table["int_end_MJD_UTC"] = fake_mjd
    table["int_start_BJD_TDB"] = fake_mjd
    table["int_mid_BJD_TDB"] = fake_mjd
    table["int_end_BJD_TDB"] = fake_mjd
    mock_file = tmp_cwd / "mock_file_uncal.fits"
    dm = datamodels.Level1bModel(data=data, refout=data, int_times=table.as_array())
    dm.meta.exposure.integration_start = 30
    dm.meta.exposure.integration_end = 40
    dm.save(mock_file)
    dm.close()

    os.chdir(tmp_cwd)
    expected_file = tmp_cwd / "mock_file_mod_uncal.fits"
    trim_tso_data(str(mock_file), 5, 33, 3)
    with datamodels.open(expected_file) as dm:
        assert dm.meta.filename == expected_file.name
        assert dm.meta.exposure.integration_start == 33
        assert dm.meta.exposure.integration_end == 38
        assert np.shape(dm.data) == (5, 10, 20, 20)
        assert np.shape(dm.refout) == (5, 10, 20, 20)
        assert len(dm.int_times) == 5


@pytest.fixture(scope="module")
def run_tso1_pipeline(rtdata_module):
    """Run the calwebb_detector1 pipeline on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module
    rtdata.get_data(INPUT_DATA["jw01536028001_03103_00001-seg001_mirimage"].full_path)

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
    rtdata.get_data(INPUT_DATA["jw01281001001_04103_00001-seg002_mirimage_mod"].full_path)

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
    rtdata.get_data(INPUT_DATA["jw01281001001_04103_00001-seg002_mirimage_mod"].full_path)

    args = [
        "calwebb_detector1",
        rtdata.input,
        f"--output_file={INPUT_DATA['jw01281001001_04103_00001-seg002_mirimage_mod'].root_name}_emijoint",
        "--steps.emicorr.algorithm=joint",
        "--steps.emicorr.save_results=True",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso_spec2_pipeline(run_tso1_pipeline, rtdata_module, resource_tracker):
    """Run the calwebb_tso-spec2 pipeline on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module

    rtdata.input = (
        f"{INPUT_DATA['jw01536028001_03103_00001-seg001_mirimage'].root_name}_rateints.fits"
    )

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
    rtdata.get_data(INPUT_DATA["jw01536028001_03103_00001-seg002_mirimage"].full_path)
    rtdata.get_data(INPUT_DATA["jw01536-o028_mod_tso3_00001_asn"].full_path)

    args = [
        "calwebb_tso3",
        INPUT_DATA["jw01536-o028_mod_tso3_00001_asn"].file_name,
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
    output_filename = (
        f"{INPUT_DATA['jw01536028001_03103_00001-seg001_mirimage'].root_name}_{step_suffix}.fits"
    )
    rtdata.output = output_filename

    rtdata.get_truth(f"rtdata/truth/test_miri_lrs_slitless_tso1/{output_filename}")

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
    file_root = INPUT_DATA["jw01281001001_04103_00001-seg002_mirimage_mod"].root_name
    output_filename = f"{file_root}_{step_suffix}.fits"
    rtdata.output = output_filename

    rtdata.get_truth(f"rtdata/truth/test_miri_lrs_slitless_detector1/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("step_suffix", ["emicorr", "rate", "rateints"])
def test_miri_lrs_slitless_detector1_emicorr_joint(
    run_detector1_pipeline_emicorr_joint, rtdata_module, fitsdiff_default_kwargs, step_suffix
):
    """Regression test of  detector1 pipeline with an alternate emicorr algorithm."""
    rtdata = rtdata_module
    file_root = INPUT_DATA["jw01281001001_04103_00001-seg002_mirimage_mod"].root_name
    output_filename = f"{file_root}_emijoint_{step_suffix}.fits"
    rtdata.output = output_filename

    rtdata.get_truth(
        f"rtdata/truth/test_miri_lrs_slitless_detector1_emicorr_joint/{output_filename}"
    )

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

    file_root = INPUT_DATA["jw01536028001_03103_00001-seg001_mirimage"].root_name
    output_filename = f"{file_root}_{step_suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"rtdata/truth/test_miri_lrs_slitless_tso_spec2/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("step_suffix", ["outlier_detection", "crfints"])
def test_miri_lrs_slitless_tso3(
    run_tso3_pipeline, rtdata_module, fitsdiff_default_kwargs, step_suffix
):
    """Compare the output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module

    file_root = INPUT_DATA["jw01536028001_03103_00001-seg001_mirimage"].root_name
    median_filename = f"{file_root}_{ASN_ID}_median.fits"
    assert os.path.isfile(median_filename)

    output_filename = f"{file_root}_{ASN_ID}_{step_suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"rtdata/truth/test_miri_lrs_slitless_tso3/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_miri_lrs_slitless_tso3_x1dints(run_tso3_pipeline, rtdata_module, fitsdiff_default_kwargs):
    """Compare the output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module

    output_filename = f"{PRODUCT_NAME}_x1dints.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"rtdata/truth/test_miri_lrs_slitless_tso3/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_miri_lrs_slitless_tso3_whtlt(run_tso3_pipeline, rtdata_module, diff_astropy_tables):
    """Compare the whitelight output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module

    output_filename = f"{PRODUCT_NAME}_whtlt.ecsv"
    rtdata.output = output_filename
    rtdata.get_truth(f"rtdata/truth/test_miri_lrs_slitless_tso3/{output_filename}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)


def test_miri_lrs_slitless_wcs(run_tso_spec2_pipeline, fitsdiff_default_kwargs, rtdata_module):
    """Compare the assign_wcs output of a MIRI LRS slitless calwebb_tso3 pipeline."""
    rtdata = rtdata_module
    output = f"{INPUT_DATA['jw01536028001_03103_00001-seg001_mirimage'].root_name}_assign_wcs.fits"
    # get input assign_wcs and truth file
    rtdata.output = output
    rtdata.get_truth("rtdata/truth/test_miri_lrs_slitless_tso_spec2/" + output)

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
    sci = rtdata.get_data(INPUT_DATA["jw04496004001_03103_00001-seg001_mirimage_mod"].full_path)

    # target acquisition verification image
    taq = rtdata.get_data(INPUT_DATA["jw04496004001_03102_00001-seg001_mirimage"].full_path)

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

    file_root = INPUT_DATA["jw04496004001_03103_00001-seg001_mirimage_mod"].root_name
    output_filename = f"{file_root}_{step_suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"rtdata/truth/test_miri_lrs_slitless_tso_spec2/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
