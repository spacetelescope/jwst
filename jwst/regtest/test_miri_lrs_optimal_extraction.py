from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
import pytest

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_spec2_optimal(rtdata_module):
    """Run the calwebb_spec2 pipeline on MIRI LRS fixedslit with optimal extraction."""
    rtdata = rtdata_module

    # Get the spec2 ASN and its members
    rtdata.get_asn("miri/lrs/jw01530-o005_20221202t204827_spec2_00001_asn.json")

    # Run the calwebb_spec2 pipeline with optimal extraction and recommended
    # parameters, saving intermediate files
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--output_file=jw01530005001_03103_00001_mirimage_opt",
        "--steps.resample_spec.skip=true",
        "--steps.pixel_replace.skip=true",
        "--steps.extract_1d.extraction_type=optimal",
        "--steps.extract_1d.use_source_posn=true",
        "--steps.extract_1d.model_nod_pair=true",
        "--steps.extract_1d.optimize_psf_location=true",
        "--steps.extract_1d.save_profile=true",
        "--steps.extract_1d.save_scene_model=true",
        "--steps.extract_1d.save_residual_image=true",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_spec3_optimal(rtdata_module):
    """Run the calwebb_spec3 pipeline on MIRI LRS fixedslit with optimal extraction."""
    rtdata = rtdata_module

    # Get the spec3 ASN and its members
    rtdata.get_asn("miri/lrs/jw01530-o005_20221202t204827_spec3_00001_asn.json")

    # Run the calwebb_spec3 pipeline with optimal extraction and recommended
    # parameters, saving intermediate files
    args = [
        "calwebb_spec3",
        rtdata.input,
        "--steps.resample_spec.skip=true",
        "--steps.pixel_replace.skip=true",
        "--steps.extract_1d.extraction_type=optimal",
        "--steps.extract_1d.use_source_posn=true",
        "--steps.extract_1d.model_nod_pair=true",
        "--steps.extract_1d.optimize_psf_location=true",
        "--steps.extract_1d.save_profile=true",
        "--steps.extract_1d.save_scene_model=true",
        "--steps.extract_1d.save_residual_image=true",
    ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["x1d", "profile", "scene_model", "residual"])
def test_miri_lrs_slit_spec2_optimal(
    run_spec2_optimal, fitsdiff_default_kwargs, rtdata_module, suffix
):
    """Regression test for MIRI LRS FS optimal extraction in spec2."""
    rtdata = rtdata_module
    output = f"jw01530005001_03103_00001_mirimage_opt_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(f"truth/test_miri_lrs_optimal_extraction/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "suffix",
    [
        "0_x1d",
        "1_x1d",
        "0_profile",
        "1_profile",
        "0_scene_model",
        "1_scene_model",
        "0_residual",
        "1_residual",
        "c1d",
    ],
)
def test_miri_lrs_slit_spec3_optimal(
    run_spec3_optimal, fitsdiff_default_kwargs, rtdata_module, suffix
):
    """Regression test for MIRI LRS FS optimal extraction in spec3."""
    rtdata = rtdata_module
    output = f"jw01530-o005_t004_miri_p750l_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(f"truth/test_miri_lrs_optimal_extraction/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
