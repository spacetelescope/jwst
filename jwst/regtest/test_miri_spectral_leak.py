"""Test SpectralLeakStep on MIRI MRS"""
import pytest
from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step

@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    ['test_spectral_leak_asn_0_spectralleakstep.fits', 'test_spectral_leak_asn_1_spectralleakstep.fits']
)
def test_miri_spectral_leak(output, _jail, rtdata, fitsdiff_default_kwargs):
    """Run cube_build on single file using coord system = ifu_align"""

    rtdata.get_asn("miri/mrs/test_spectral_leak_asn.json")
    args = ["jwst.spectral_leak.SpectralLeakStep", rtdata.input]
    Step.from_cmdline(args)

    # The spectral_leak correction is part of calwebb_spec3. The top program, calwebb_spec3.py,
    # controls writing the file names according to the assocation rules.
    # For this test we just want to isolate the spectral_leak correction and run the extracted spectrum
    # through the spectral_leak correction. This step does not have the smarts to know how to write
    # the data. The filenames are written as:
    # 1. test_spectral_leak_asn_0_spectralleakstep.fits which is the CH1 short medium extracted spectrum (this should
    # not be changed from the input file).
    # 2. test_spectral_leak_asn_1_spectralleakstep.fits which is the corrected CH3 short extracted spectrum

    rtdata.output = f"{output}"
    # Get the truth files
    rtdata.get_truth(f"truth/test_miri_spectral_leak/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
