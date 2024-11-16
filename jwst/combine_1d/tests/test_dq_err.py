"""
Test DQ operations in Combine1dStep
"""
import numpy as np

from stdatamodels.jwst import datamodels

from jwst.combine_1d import Combine1dStep


def test_dq():
    """Test that DQ exclusion works"""
    spec1 = create_spec_model(flux=1e-9)
    spec2 = create_spec_model(flux=1e-9)
    ms = datamodels.MultiSpecModel()
    ms.meta.exposure.exposure_time = 1
    ms.spec.append(spec1)
    ms.spec.append(spec2)

    # Make one pixel bad by being large.
    # Result should be large
    BAD_PIX = 5
    spec1.spec_table['FLUX'][BAD_PIX] = 1.0
    result = Combine1dStep.call(ms)
    assert np.isclose(result.spec[0].spec_table['FLUX'][BAD_PIX], 0.5)

    # Now mark that pixel bad.
    # Result should not just contain the second spectrum value.
    spec1.spec_table['DQ'][BAD_PIX] = datamodels.dqflags.pixel['DO_NOT_USE']
    result_dq = Combine1dStep.call(ms)
    assert np.isclose(result_dq.spec[0].spec_table['FLUX'][BAD_PIX], spec2.spec_table['FLUX'][BAD_PIX])


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
    assert np.allclose(result.spec[0].spec_table['FLUX'], 1.0)
    assert np.allclose(result.spec[0].spec_table['ERROR'], expected_error)
    assert np.allclose(result.spec[0].spec_table['SURF_BRIGHT'], 1.0)
    assert np.allclose(result.spec[0].spec_table['SB_ERROR'], expected_error)


def create_spec_model(npoints=10, flux=1e-9, error=1e-10, wave_range=(11, 13)):
    """Create a SpecModel"""

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

    spec_dtype = datamodels.SpecModel().spec_table.dtype  # This data type is used for creating an output table.

    otab = np.array(
        list(
            zip(
                wavelength, flux, error, var, var, var,
                surf_bright, sb_error, var, var, var,
                dq, background, berror, var, var, var,
                npixels
            ),
        ), dtype=spec_dtype
    )

    spec_model = datamodels.SpecModel(spec_table=otab)

    return spec_model
