import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run the calwebb_spec2 pipeline on a single NIRSpec MOS/FS exposure."""

    rtdata = rtdata_module

    # Get the MSA metadata file referenced in the input exposure
    rtdata.get_data("nirspec/mos/jw02674004001_01_msa.fits")

    # Get the input ASN file and exposures
    rtdata.get_data("nirspec/mos/jw02674004001_03101_00001_nrs1_rate.fits")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = ["calwebb_spec2", rtdata.input,
            "--steps.assign_wcs.save_results=true",
            "--steps.msa_flagging.save_results=true",
            "--steps.master_background_mos.save_results=true",
            "--steps.extract_2d.save_results=true",
            "--steps.srctype.save_results=true",
            "--steps.wavecorr.save_results=true",
            "--steps.flat_field.save_results=true",
            "--steps.pathloss.save_results=true",
            "--steps.barshadow.save_results=true"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", [
    "assign_wcs", "msa_flagging", "extract_2d", "srctype",
    "master_background_mos", "wavecorr", "flat_field", "pathloss", "barshadow",
    "wavecorr_fs", "flat_field_fs", "pathloss_fs",
    "cal", "s2d", "x1d"])
def test_nirspec_mos_fs_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test for calwebb_spec2 on a NIRSpec MOS/FS exposure."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = f"jw02674004001_03101_00001_nrs1_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nirspec_mos_fs_spec2/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
