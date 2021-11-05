import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step

from jwst.ami import AmiAnalyzeStep


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_ami3 on NIRISS AMI data."""
    rtdata = rtdata_module
    rtdata.get_asn("niriss/ami/jw00793-c1014_20191210t203450_ami3_001_asn.json")

    # Run the calwebb_ami3 pipeline on the association
    args = ["calwebb_ami3", rtdata.input,
            "--steps.ami_analyze.rotation=1.49",
            ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["c1014_ami"])
@pytest.mark.parametrize("exposure", ["022", "025"])
def test_niriss_ami3_exp(run_pipeline, suffix, exposure, fitsdiff_default_kwargs):
    """Check exposure-level results of calwebb_ami3"""
    rtdata = run_pipeline

    output = "jw00793" + exposure + "001_03102_00001_nis_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_niriss_ami3/" + output)

    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["amiavg", "psf-amiavg", "aminorm"])
def test_niriss_ami3_product(run_pipeline, suffix, fitsdiff_default_kwargs):
    """Check final products of calwebb_ami3"""
    rtdata = run_pipeline

    output = "jw00793-c1014_t005_niriss_f480m-nrm-sub80_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_niriss_ami3/" + output)

    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_ami_analyze_with_nans(rtdata, fitsdiff_default_kwargs):
    """Test the AmiAnalyzeStep using an input image with NaNs"""
    data = rtdata.get_data('niriss/ami/jw00042004001_01101_00005_nis_withNAN_cal.fits')

    AmiAnalyzeStep.call(data, save_results=True)
    rtdata.output = 'jw00042004001_01101_00005_nis_withNAN_amianalyzestep.fits'

    rtdata.get_truth('truth/test_niriss_ami3/jw00042004001_01101_00005_nis_withNAN_amianalyzestep.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
