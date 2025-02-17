"""Test operations in Combine1dStep."""

import numpy as np
import pytest

from jwst import datamodels
from jwst.combine_1d import Combine1dStep
from jwst.combine_1d.combine1d import InputSpectrumModel, check_exptime


@pytest.fixture
def two_spectra():
    spec1 = create_spec_model(flux=1e-9)
    spec2 = create_spec_model(flux=1e-9)
    ms = datamodels.MultiSpecModel()
    ms.meta.exposure.exposure_time = 1.0
    ms.meta.exposure.integration_time = 2.0
    ms.spec.append(spec1)
    ms.spec.append(spec2)
    yield ms
    ms.close()


@pytest.fixture
def three_spectra():
    spec1 = create_spec_model(flux=1e-9)
    spec2 = create_spec_model(flux=1e-9)
    spec3 = create_spec_model(flux=1e-9)
    ms = datamodels.MultiSpecModel()
    ms.meta.exposure.exposure_time = 1
    ms.spec.append(spec1)
    ms.spec.append(spec2)
    ms.spec.append(spec3)
    yield ms
    ms.close()


def test_dq(two_spectra):
    """Test that DQ exclusion works."""
    ms = two_spectra
    spec1 = ms.spec[0]
    spec2 = ms.spec[1]

    # Make one pixel bad by being large.
    # Result should be large
    bad_pix = 5
    spec1.spec_table["FLUX"][bad_pix] = 1.0
    result = Combine1dStep.call(ms)
    assert np.isclose(result.spec[0].spec_table["FLUX"][bad_pix], 0.5)

    # Now mark that pixel bad.
    # Result should not just contain the second spectrum value.
    spec1.spec_table["DQ"][bad_pix] = datamodels.dqflags.pixel["DO_NOT_USE"]
    result_dq = Combine1dStep.call(ms)
    assert np.isclose(
        result_dq.spec[0].spec_table["FLUX"][bad_pix], spec2.spec_table["FLUX"][bad_pix]
    )


def test_err():
    """Test error propagation."""
    spec1 = create_spec_model(flux=1.0, error=0.1)
    spec2 = create_spec_model(flux=1.0, error=0.2)
    ms = datamodels.MultiSpecModel()
    ms.meta.exposure.exposure_time = 1
    ms.spec.append(spec1)
    ms.spec.append(spec2)

    result = Combine1dStep.call(ms)

    expected_error = np.sqrt(0.1**2 + 0.2**2) / 2
    assert np.allclose(result.spec[0].spec_table["FLUX"], 1.0)
    assert np.allclose(result.spec[0].spec_table["ERROR"], expected_error)
    assert np.allclose(result.spec[0].spec_table["SURF_BRIGHT"], 1.0)
    assert np.allclose(result.spec[0].spec_table["SB_ERROR"], expected_error)


def test_nan(two_spectra):
    """Test that nan exclusion works."""
    ms = two_spectra
    spec1 = ms.spec[0]
    spec2 = ms.spec[1]

    # Make one pixel bad by being nan, without being marked as DO_NOT_USE.
    # Result should just contain the second spectrum value
    bad_pix = 5
    spec1.spec_table["FLUX"][bad_pix] = np.nan
    result = Combine1dStep.call(ms)
    assert np.isclose(result.spec[0].spec_table["FLUX"][bad_pix], spec2.spec_table["FLUX"][bad_pix])


def test_sigmaclip(three_spectra):
    """Test that sigma clipping works."""
    ms = three_spectra
    spec1 = ms.spec[0]
    spec2 = ms.spec[1]

    # Make one pixel bad by being large.
    # Result should be large
    bad_pix = 5
    spec1.spec_table["FLUX"][bad_pix] = 1.0
    result = Combine1dStep.call(ms)
    assert np.isclose(result.spec[0].spec_table["FLUX"][bad_pix], 1 / 3)
    result.close()

    # Now process with sigma_clipping turned on
    # Result should not just contain the second and third spectra values.
    result_sc = Combine1dStep.call(ms, sigma_clip=3)
    assert np.isclose(
        result_sc.spec[0].spec_table["FLUX"][bad_pix], spec2.spec_table["FLUX"][bad_pix]
    )
    result_sc.close()


@pytest.mark.parametrize("casing", ["upper", "lower"])
@pytest.mark.parametrize(
    "exptime",
    ["exposure_time", "integration_time", "unit_weight", "effexptm", "effinttm", "unit weight"],
)
def test_exptime_keys(exptime, casing):
    if casing == "upper":
        exptime = exptime.upper()

    spec1 = create_spec_model(flux=1.0, error=0.1)
    spec2 = create_spec_model(flux=2.0, error=0.2)

    ms1 = datamodels.MultiSpecModel()
    ms1.meta.exposure.exposure_time = 1.0
    ms1.meta.exposure.integration_time = 2.0
    ms1.spec.append(spec1)

    ms2 = datamodels.MultiSpecModel()
    ms2.meta.exposure.exposure_time = 2.0
    ms2.meta.exposure.integration_time = 1.0
    ms2.spec.append(spec2)

    result = Combine1dStep.call(datamodels.ModelContainer([ms1, ms2]), exptime_key=exptime)
    if "exp" in exptime.lower():
        # closer to 2
        assert np.allclose(result.spec[0].spec_table["FLUX"], 1 + 2 / 3)
    elif "int" in exptime.lower():
        # closer to 1
        assert np.allclose(result.spec[0].spec_table["FLUX"], 1 + 1 / 3)
    else:
        # equidistant
        assert np.allclose(result.spec[0].spec_table["FLUX"], 1.5)


def test_bad_exptime(two_spectra):
    # Runtime error if bad key is passed to input model
    with pytest.raises(RuntimeError):
        InputSpectrumModel(two_spectra, two_spectra.spec[0], "bad")

    # Bad key is translated to unit_weight if checked
    assert check_exptime("bad") == "unit_weight"


def create_spec_model(npoints=10, flux=1e-9, error=1e-10, wave_range=(11, 13)):
    """Create a SpecModel."""
    wavelength = np.arange(*wave_range, step=(wave_range[1] - wave_range[0]) / npoints)
    flux = np.full(npoints, flux)
    error = np.full(npoints, error)

    surf_bright = np.full(npoints, flux)
    sb_error = np.full(npoints, error)

    var = np.zeros(npoints)
    dq = np.zeros(npoints)
    background = np.zeros(npoints)
    berror = np.zeros(npoints)
    npixels = np.zeros(npoints)

    # This data type is used for creating an output table.
    spec_dtype = datamodels.SpecModel().spec_table.dtype

    otab = np.array(
        list(
            zip(
                wavelength,
                flux,
                error,
                var,
                var,
                var,
                surf_bright,
                sb_error,
                var,
                var,
                var,
                dq,
                background,
                berror,
                var,
                var,
                var,
                npixels, strict=False,
            ),
        ),
        dtype=spec_dtype,
    )

    spec_model = datamodels.SpecModel(spec_table=otab)

    return spec_model
