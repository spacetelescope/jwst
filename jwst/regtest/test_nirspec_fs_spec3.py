from astropy.io.fits.diff import FITSDiff
import pytest

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """
    Run the calwebb_spec3 pipeline on NIRSpec Fixed-Slit exposures.
    """
    rtdata = rtdata_module

    # Get the cfg files
    collect_pipeline_cfgs("config")

    # Get the ASN file and input exposures
    rtdata.get_asn('nirspec/fs/jw93045-o010_20180725t035735_spec3_001_asn.json')

    # Run the calwebb_spec3 pipeline; save results from intermediate steps
    args = ["config/calwebb_spec3.cfg", rtdata.input,
            "--steps.outlier_detection.save_results=true",
            "--steps.resample_spec.save_results=true",
            "--steps.extract_1d.save_results=true"]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
def test_nirspec_fs_spec3(run_pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test spec3 pipeline on a set of NIRSpec FS exposures."""
    rtdata = rtdata_module

    output = f"jw93045-o010_s00003_nirspec_f290lp-g395h-subs400a1_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(f"truth/test_nirspec_fs_spec3/{output}")

    # Set looser tolerance for resampled s2d comparison
    if suffix == "s2d":
        fitsdiff_default_kwargs["atol"] = 1e-5

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
