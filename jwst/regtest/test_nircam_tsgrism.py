import pytest

from jwst.lib.set_telescope_pointing import add_wcs
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_spec2_pipeline(rtdata_module, resource_tracker):
    """Run stage 2 pipeline on NIRCAM TSO grism data."""
    rtdata = rtdata_module

    # Run spec2 pipeline on the _rateints file, saving intermediate products
    rtdata.get_data("nircam/tsgrism/jw01366002001_04103_00001-seg001_nrcalong_rateints.fits")
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.flat_field.save_results=True",
        "--steps.extract_2d.save_results=True",
        "--steps.srctype.save_results=True",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


@pytest.fixture(scope="module")
def run_tso3_pipeline(rtdata_module, run_spec2_pipeline, resource_tracker):
    """Run stage 3 pipeline on NIRCAM TSO grism data."""
    rtdata = rtdata_module

    # Get the level3 association json file (though not its members) and run
    # the tso3 pipeline on all _calints files listed in association
    rtdata.get_data("nircam/tsgrism/jw01366-o002_20230107t004627_tso3_00001_asn.json")
    args = ["calwebb_tso3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


@pytest.fixture(
    scope="module",
    params=[
        "jw01185015001_03104_00001-seg004_nrcalong_rate.fits",
        "jw01185013001_04103_00001-seg003_nrcalong_rate.fits",
    ],
)
def run_pipeline_offsetSR(request, rtdata_module):
    rtdata = rtdata_module
    rtdata.get_data("nircam/tsgrism/" + request.param)
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.extract_1d.save_results=True",
    ]
    Step.from_cmdline(args)
    return rtdata


def test_log_tracked_resources_spec2(log_tracked_resources, run_spec2_pipeline):
    log_tracked_resources()


def test_log_tracked_resources_tso3(log_tracked_resources, run_tso3_pipeline):
    log_tracked_resources()


def test_nircam_tsgrism_stage2_offsetSR(run_pipeline_offsetSR, fitsdiff_default_kwargs):
    """
    Test coverage for offset special requirement specifying nonzero offset in X.

    Test data are two observations of Gliese 436, one with offset specified and one without.
    Quantitatively we just ensure that the outputs are identical to the inputs, but qualitatively
    can check that the spectral lines fall at the same wavelengths in both cases,
    which is why the zero-offset case is also included here."""
    rtdata = run_pipeline_offsetSR
    rtdata.output = rtdata.input.replace("rate", "x1d")
    rtdata.get_truth("truth/test_nircam_tsgrism_stages/" + rtdata.output.split("/")[-1])

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "suffix", ["calints", "extract_2d", "flat_field", "o002_crfints", "srctype", "x1dints"]
)
def test_nircam_tsgrism_stage2(run_spec2_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of tso-spec2 pipeline performed on NIRCam TSO grism data."""
    rtdata = run_spec2_pipeline
    rtdata.input = "jw01366002001_04103_00001-seg001_nrcalong_rateints.fits"
    output = "jw01366002001_04103_00001-seg001_nrcalong_" + suffix + ".fits"
    rtdata.output = output

    rtdata.get_truth("truth/test_nircam_tsgrism_stages/" + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nircam_tsgrism_stage3_x1dints(run_tso3_pipeline, fitsdiff_default_kwargs):
    rtdata = run_tso3_pipeline
    rtdata.input = "jw01366-o002_20230107t004627_tso3_00001_asn.json"
    rtdata.output = "jw01366-o002_t001_nircam_f322w2-grismr-subgrism256_x1dints.fits"
    rtdata.get_truth(
        "truth/test_nircam_tsgrism_stages/jw01366-o002_t001_nircam_f322w2-grismr-subgrism256_x1dints.fits"
    )

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nircam_tsgrism_stage3_whtlt(run_tso3_pipeline, diff_astropy_tables):
    rtdata = run_tso3_pipeline
    rtdata.input = "jw01366-o002_20230107t004627_tso3_00001_asn.json"
    rtdata.output = "jw01366-o002_t001_nircam_f322w2-grismr-subgrism256_whtlt.ecsv"
    rtdata.get_truth(
        "truth/test_nircam_tsgrism_stages/jw01366-o002_t001_nircam_f322w2-grismr-subgrism256_whtlt.ecsv"
    )

    assert diff_astropy_tables(rtdata.output, rtdata.truth)


def test_nircam_setpointing_tsgrism(rtdata, fitsdiff_default_kwargs):
    """
    Regression test of the set_telescope_pointing script on a level-1b NIRCam file.
    """
    rtdata.get_data("nircam/tsgrism/jw02459001001_03103_00001-seg001_nrcalong_uncal.fits")
    # The add_wcs function overwrites its input
    rtdata.output = rtdata.input

    # Call the WCS routine
    add_wcs(rtdata.input)

    rtdata.get_truth(
        "truth/test_nircam_setpointing/jw02459001001_03103_00001-seg001_nrcalong_uncal.fits"
    )

    fitsdiff_default_kwargs["rtol"] = 1e-6
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
