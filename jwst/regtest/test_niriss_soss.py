import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_tso_spec2(jail, rtdata_module):
    """Run stage 2 pipeline on NIRISS SOSS data."""
    rtdata = rtdata_module

    # Run tso-spec2 pipeline on the first _rateints file, saving intermediate products
    rtdata.get_data("niriss/soss/jw01091002001_03101_00001-seg001_nis_short_rateints.fits")
    args = ["calwebb_spec2", rtdata.input,
            "--steps.flat_field.save_results=True",
            "--steps.srctype.save_results=True",
            "--steps.extract_1d.soss_atoca=False",
            ]
    Step.from_cmdline(args)

    # Run tso-spec2 pipeline on the second _rateints file, without saving or
    # checking any results (simply create a fresh input for level-3 test)
    rtdata.get_data("niriss/soss/jw01091002001_03101_00001-seg002_nis_short_rateints.fits")
    args = ["calwebb_spec2", rtdata.input,
            "--steps.extract_1d.soss_atoca=False",
            ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso_spec3(jail, rtdata_module, run_tso_spec2):
    """Run stage 3 pipeline on NIRISS SOSS data."""
    rtdata = rtdata_module
    # Get the level3 association json file (though not its members) and run
    # the tso3 pipeline on all _calints files listed in association
    rtdata.get_data("niriss/soss/jw01091-o002_20220714t155100_tso3_001_asn.json")
    args = ["calwebb_tso3", rtdata.input,
            "--steps.extract_1d.soss_rtol=1.e-3",
            ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_atoca_extras(jail, rtdata_module):
    """Run stage 2 pipeline on NIRISS SOSS data using enhanced modes via parameter settings."""
    rtdata = rtdata_module

    # Run spec2 pipeline on the second _rateints file, using wavegrid generated from first segment.
    rtdata.get_data("niriss/soss/seg001_wavegrid.fits")
    rtdata.get_data("niriss/soss/atoca_extras_rateints.fits")
    args = ["calwebb_spec2", rtdata.input,
            "--steps.extract_1d.soss_modelname=atoca_extras",
            "--steps.extract_1d.soss_wave_grid_in=seg001_wavegrid.fits",
            "--steps.extract_1d.soss_bad_pix=model",
            "--steps.extract_1d.soss_rtol=1.e-3",
            ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["calints", "flat_field", "srctype", "x1dints"])
def test_niriss_soss_stage2(rtdata_module, run_tso_spec2, fitsdiff_default_kwargs, suffix):
    """Regression test of tso-spec2 pipeline performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = f"jw01091002001_03101_00001-seg001_nis_short_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_niriss_soss_stage3_crfints(rtdata_module, run_tso_spec3, fitsdiff_default_kwargs):
    """Regression test of tso-spec3 pipeline outlier_detection results performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = "jw01091002001_03101_00001-seg001_nis_short_o002_crfints.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_niriss_soss_stage3_x1dints(run_tso_spec3, rtdata_module, fitsdiff_default_kwargs):
    """Regression test of tso-spec3 pipeline extract_1d results performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = "jw01091-o002_t001_niriss_clear-gr700xd-substrip256_x1dints.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_niriss_soss_stage3_whtlt(run_tso_spec3, rtdata_module, diff_astropy_tables):
    """Regression test of tso-spec3 pipeline white_light results performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = "jw01091-o002_t001_niriss_clear-gr700xd-substrip256_whtlt.ecsv"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["calints", "x1dints", "AtocaSpectra", "SossExtractModel"])
def test_niriss_soss_extras(rtdata_module, run_atoca_extras, fitsdiff_default_kwargs, suffix):
    """Regression test of ATOCA enhanced algorithm performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = f"atoca_extras_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
