import pathlib

import numpy as np
import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from stdatamodels.jwst.datamodels import SossWaveGridModel

from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_tso_spec2(rtdata_module):
    """Run stage 2 pipeline on NIRISS SOSS data."""
    rtdata = rtdata_module

    # Run tso-spec2 pipeline on the first _rateints file, saving intermediate products
    rtdata.get_data("niriss/soss/jwst_niriss_soss_bkg_sub256.fits")
    rtdata.get_data("niriss/soss/jw01091002001_03101_00001-seg001_nis_short_rateints.fits")
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.bkg_subtract.save_results=True",
        "--steps.flat_field.save_results=True",
        "--steps.srctype.save_results=True",
        "--steps.extract_1d.soss_atoca=False",
    ]
    Step.from_cmdline(args)

    # Run tso-spec2 pipeline on the second _rateints file, without saving or
    # checking any results (simply create a fresh input for level-3 test)
    rtdata.get_data("niriss/soss/jw01091002001_03101_00001-seg002_nis_short_rateints.fits")
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.extract_1d.soss_atoca=False",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso_spec3(rtdata_module, run_tso_spec2, resource_tracker):
    """Run stage 3 pipeline on NIRISS SOSS data."""
    rtdata = rtdata_module
    # Get the level3 association json file (though not its members) and run
    # the tso3 pipeline on all _calints files listed in association
    rtdata.get_data("niriss/soss/jw01091-o002_20220714t155100_tso3_001_asn.json")
    args = [
        "calwebb_tso3",
        rtdata.input,
        "--steps.extract_1d.soss_rtol=1.e-3",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_atoca_extras(rtdata_module, resource_tracker):
    """Run stage 2 pipeline on NIRISS SOSS data using enhanced modes via parameter settings."""
    rtdata = rtdata_module

    # Run spec2 pipeline on the second _rateints file, using wavegrid generated from first segment.
    rtdata.get_data("niriss/soss/jw01091002001_03101_00001-seg001_wavegrid.fits")
    rtdata.get_data("niriss/soss/jw01091002001_03101_00001-seg002_nis_short_rateints.fits")
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--output_file=atoca_extras",
        "--steps.extract_1d.soss_modelname=atoca_extras",
        "--steps.extract_1d.soss_wave_grid_in=jw01091002001_03101_00001-seg001_wavegrid.fits",
        "--steps.extract_1d.soss_bad_pix=model",
        "--steps.extract_1d.soss_rtol=1.e-3",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_spec2(log_tracked_resources, run_atoca_extras):
    log_tracked_resources()


def test_log_tracked_resources_spec3(log_tracked_resources, run_tso_spec3):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["calints", "flat_field", "srctype", "x1dints"])
def test_niriss_soss_stage2(rtdata_module, run_tso_spec2, fitsdiff_default_kwargs, suffix):
    """Regression test of tso-spec2 pipeline performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = f"jw01091002001_03101_00001-seg001_nis_short_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    # Ignore the custom bkg reference file because it contains a full path.
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_niriss_soss_stage3_crfints(rtdata_module, run_tso_spec3, fitsdiff_default_kwargs):
    """Regression test of tso3pipeline outlier_detection results performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = "jw01091002001_03101_00001-seg001_nis_short_o002_crfints.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    # Ignore the custom bkg reference file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_BKG")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_niriss_soss_stage3_x1dints(run_tso_spec3, rtdata_module, fitsdiff_default_kwargs):
    """Regression test of tso-spec3 pipeline extract_1d results performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = "jw01091-o002_t001_niriss_clear-gr700xd-substrip256_x1dints.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    # Ignore the custom bkg reference file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_BKG")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_niriss_soss_stage3_whtlt(run_tso_spec3, rtdata_module, diff_astropy_tables):
    """Regression test of tso-spec3 pipeline white_light results performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = "jw01091-o002_t001_niriss_clear-gr700xd-substrip256_whtlt.ecsv"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)


@pytest.mark.parametrize("suffix", ["calints", "x1dints", "AtocaSpectra", "SossExtractModel"])
def test_niriss_soss_extras(rtdata_module, run_atoca_extras, fitsdiff_default_kwargs, suffix):
    """Regression test of ATOCA enhanced algorithm performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = f"atoca_extras_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    if suffix == "AtocaSpectra":
        # Supplemental output from atoca may have system dependent diffs
        # that can't be reasonably compared. Just check for existence.
        assert pathlib.Path(rtdata.output).exists()
    else:
        diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
        assert diff.identical, diff.report()


@pytest.fixture(scope="module")
def run_extract1d_null_order2(rtdata_module):
    """
    Test coverage for fix to error thrown when all of the pixels
    in order 2 are flagged as bad. Ensure graceful handling of the
    MaskOverlapError exception raise.
    Pin tikfac and transform for faster runtime
    """
    rtdata = rtdata_module
    rtdata.get_data("niriss/soss/jw01201008001_04101_00001-seg003_nis_int72.fits")
    args = [
        "extract_1d",
        rtdata.input,
        "--soss_tikfac=4.290665733550672e-17",
    ]
    Step.from_cmdline(args)


def test_extract1d_null_order2(rtdata_module, run_extract1d_null_order2, fitsdiff_default_kwargs):
    rtdata = rtdata_module

    output = "jw01201008001_04101_00001-seg003_nis_int72_extract1dstep.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.fixture(scope="module")
def run_spec2_substrip96(rtdata_module):
    """
    Run stage 2 pipeline on substrip96 data.

    Solving for the optimal Tikhonov factor is time-consuming, and the code to do
    so is identical between substrip96 and substrip256 data. Therefore just set
    it to a reasonable value here.

    Similarly, computing the wave_grid is already tested for other modes and is also
    time-consuming. Therefore, set it to be just 1000 evenly spaced grid points here.
    """
    rtdata = rtdata_module
    rtdata.get_data("niriss/soss/jw03596001001_03102_00001-seg001_nis_ints0-2_rateints.fits")

    # further reduce runtime by setting the input wave_grid instead of calculating it
    # this also serves to test that input parameter
    wave_fname = "jw03596001001_wavegrid.fits"
    wave_grid = np.linspace(0.55, 2.8, 1000)
    wave_grid = SossWaveGridModel(wavegrid=wave_grid)
    wave_grid.save(wave_fname)

    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.extract_1d.soss_tikfac=1.0e-16",
        "--steps.extract_1d.soss_wave_grid_in=jw03596001001_wavegrid.fits",
    ]
    Step.from_cmdline(args)


@pytest.mark.parametrize("suffix", ["calints", "x1dints"])
def test_spec2_substrip96(rtdata_module, run_spec2_substrip96, fitsdiff_default_kwargs, suffix):
    """Regression test of tso-spec2 pipeline performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = f"jw03596001001_03102_00001-seg001_nis_ints0-2_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
