import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run the calwebb_spec2 pipeline on a single NIRSpec MOS exposure."""

    rtdata = rtdata_module

    # Get the MSA metadata file referenced in the input exposure
    rtdata.get_data("nirspec/mos/jw01180025001_01_msa.fits")

    # Get the input ASN file and exposures
    rtdata.get_asn("nirspec/mos/jw01180-o025_20221129t204108_spec2_00037_asn.json")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = ["calwebb_spec2", rtdata.input,
            "--steps.assign_wcs.save_results=true",
            "--steps.msa_flagging.save_results=true",
            "--steps.nsclean.skip=False",
            "--steps.nsclean.save_results=true",
            "--steps.bkg_subtract.save_results=true",
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
    "assign_wcs", "msa_flagging", "nsclean", "bsub", "extract_2d", "wavecorr", "flat_field",
    "srctype", "pathloss", "barshadow", "cal", "s2d", "x1d"])
def test_nirspec_mos_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline on a
       NIRSpec MOS exposure."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = f"jw01180025001_05101_00001_nrs2_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nirspec_mos_spec2/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
