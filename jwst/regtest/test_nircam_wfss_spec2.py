import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module, resource_tracker):
    """Run the calwebb_spec2 pipeline on a single NIRSpec WFSS exposure."""

    rtdata = rtdata_module

    # Get the input data; load individual data files first, load ASN file last
    rtdata.get_data("nircam/wfss/jw01076-o103_t001_nircam_clear-f250m_cat.ecsv")
    rtdata.get_data("nircam/wfss/jw01076103001_02104_00001_nrcblong_rate.fits")
    rtdata.get_data("nircam/wfss/jw01076-o103_20231023t160322_spec2_00003_asn.json")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.bkg_subtract.skip=False",
        "--steps.bkg_subtract.wfss_mmag_extract=20.0",
        "--steps.bkg_subtract.save_results=true",
        "--steps.extract_2d.save_results=true",
        "--steps.extract_2d.wfss_mmag_extract=19.0",
        "--steps.extract_2d.wfss_nbright=20",
        "--steps.srctype.save_results=true",
        "--steps.flat_field.save_results=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_spec2(log_tracked_resources, run_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix", ["assign_wcs", "bsub", "extract_2d", "flat_field", "srctype", "cal", "x1d"]
)
def test_nircam_wfss_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline on a
    NIRCam WFSS exposure."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = f"jw01076103001_02104_00001_nrcblong_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nircam_wfss_spec2/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.fixture
def run_pipeline_select_sources(rtdata_module):
    """
    Run calwebb_spec2 pipeline enabling source selection in extract_2d step.

    Skip background step for runtime.
    Include wfss_contam step since it occurs downstream of extract_2d.
    """
    rtdata = rtdata_module

    # Get the input data; load individual data files first, load ASN file last
    rtdata.get_data("nircam/wfss/jw01076-o101_t002_nircam_clear-f356w_cat.ecsv")
    from astropy.table import QTable

    cat = QTable.read(rtdata.input, format="ascii.ecsv")
    centroid1 = cat[423]["sky_centroid"]  # source id 424
    centroid2 = cat[1359]["sky_centroid"]  # source id 1360
    bad_centroid = (0.0, 0.0)  # not on detector

    rtdata.get_data("nircam/wfss/jw01076-o101_t002_nircam_clear-f356w_segm.fits")
    rtdata.get_data("nircam/wfss/jw01076-o101_t002_nircam_clear-f356w_i2d.fits")
    rtdata.get_data("nircam/wfss/jw01076101001_02101_00003_nrcalong_rate.fits")
    rtdata.get_data("nircam/wfss/jw01076-o101_20220403t120233_spec2_002_asn.json")

    # note wfss_contam mag limit does not affect source ID selection.
    # just makes the step model fewer sources and thereby run faster
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--output_file=four_sources_cal.fits",  # avoid clobbering the test above
        "--steps.extract_1d.output_file=four_sources_x1d.fits",  # avoid clobbering the test above
        "--steps.bkg_subtract.skip=True",
        "--steps.extract_2d.source_ids=202,2157",  # some source ids that are actually on detector
        f"--steps.extract_2d.source_ra={centroid1.ra.value},{centroid2.ra.value},{bad_centroid[0]}",
        f"--steps.extract_2d.source_dec={centroid1.dec.value},{centroid2.dec.value},{bad_centroid[1]}",
        "--steps.extract_2d.source_max_sep=1.0",
        "--steps.wfss_contam.skip=False",
        "--steps.wfss_contam.magnitude_limit=18.0",
    ]

    Step.from_cmdline(args)
    return rtdata


def test_nircam_wfss_spec2_select_sources(run_pipeline_select_sources):
    """
    Test calwebb_spec2 pipeline enabling source selection in extract_2d step.

    Ensure this source selection doesn't cause issues downstream, and that we end up
    with those sources in the output x1d file.
    """

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline_select_sources
    output = "four_sources_x1d.fits"
    rtdata.output = output
    with dm.open(rtdata.output) as model:
        assert len(model.spec[0].spec_table) == 4
        output_ids = model.spec[0].spec_table["SOURCE_ID"]
        assert 202 in output_ids
        assert 2157 in output_ids
        assert 424 in output_ids
        assert 1360 in output_ids
