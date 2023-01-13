from astropy.io.fits.diff import FITSDiff
import pytest

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """
    Run the calwebb_spec3 pipeline on NIRSpec Fixed-Slit exposures.
    """
    rtdata = rtdata_module

    # Get the ASN file and input exposures
    rtdata.get_asn('nirspec/fs/jw01309-o022_spec3_regtest_asn.json')

    # Run the calwebb_spec3 pipeline; save results from intermediate steps
    args = ["calwebb_spec3", rtdata.input,
            "--steps.outlier_detection.save_results=true",
            "--steps.resample_spec.save_results=true",
            "--steps.extract_1d.save_results=true"]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
@pytest.mark.parametrize("source_id", ["s00001", "s00002", "s00003", "s00004", "s00005"])
def test_nirspec_fs_spec3(run_pipeline, rtdata_module, fitsdiff_default_kwargs, suffix, source_id):
    """Test spec3 pipeline on a set of NIRSpec FS exposures."""
    rtdata = rtdata_module

    output = f"jw01309-o022_{source_id}_nirspec_f290lp-g395h-s200a2-allslits_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_fs_spec3/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "s2d":
        fitsdiff_default_kwargs["rtol"] = 1e-2
        fitsdiff_default_kwargs["atol"] = 2e-4

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
