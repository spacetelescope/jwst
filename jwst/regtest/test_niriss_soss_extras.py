import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_atoca_extras(rtdata_module):
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
@pytest.mark.parametrize("suffix", ["calints", "x1dints", "AtocaSpectra", "SossExtractModel"])
def test_niriss_soss_extras(rtdata_module, run_atoca_extras, fitsdiff_default_kwargs, suffix):
    """Regression test of ATOCA enhanced algorithm performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = f"atoca_extras_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
