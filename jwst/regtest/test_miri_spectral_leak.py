"""Test SpectralLeakStep on MIRI MRS"""

import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "output",
    [
        "test_spectral_leak_asn_0_spectralleakstep.fits",
        "test_spectral_leak_asn_1_spectralleakstep.fits",
    ],
)
def test_miri_spectral_leak(output, rtdata, fitsdiff_default_kwargs):
    """Run cube_build on single file using coord system = ifu_align"""

    # this data does not have the residual fringe corrected columns in the x1d files
    rtdata.get_asn("miri/mrs/test_spectral_leak_asn.json")
    args = ["jwst.spectral_leak.SpectralLeakStep", rtdata.input]
    Step.from_cmdline(args)

    # The spectral_leak correction is part of calwebb_spec3. The top program, calwebb_spec3.py,
    # controls writing the file names according to the association rules.
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


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "output",
    [
        "test_spectral_leak_rf_corrected_asn_0_spectralleakstep.fits",
        "test_spectral_leak_rf_corrected_asn_1_spectralleakstep.fits",
    ],
)
def test_miri_spectral_leak_rf(output, rtdata, fitsdiff_default_kwargs):
    """Run cube_build on single file using coord system = ifu_align"""

    # this data has the residual fringe corrected columns in the x1d files
    rtdata.get_asn("miri/mrs/test_spectral_leak_rf_corrected_asn.json")
    args = ["jwst.spectral_leak.SpectralLeakStep", rtdata.input]
    Step.from_cmdline(args)

    # The spectral_leak correction is part of calwebb_spec3. The top program, calwebb_spec3.py,
    # controls writing the file names according to the association rules.
    # For this test we just want to isolate the spectral_leak correction and run the extracted spectrum
    # through the spectral_leak correction. This step does not have the smarts to know how to write
    # the data. The filenames are written as:
    # 1. test_spectral_leak_rf_corrected_asn_0_spectralleakstep.fits which is the CH1 short medium extracted spectrum (this should
    # not be changed from the input file).
    # 2. test_spectral_leak_rf_corrected_asn_1_spectralleakstep.fits which is the corrected CH3 short extracted spectrum

    rtdata.output = f"{output}"
    # Get the truth files
    rtdata.get_truth(f"truth/test_miri_spectral_leak/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
