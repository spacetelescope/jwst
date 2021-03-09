"""
Test DQ operations in Combine1dStep
"""
import numpy as np

from jwst import datamodels
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


def create_spec_model(npoints=10, flux=1e-9, wave_range=(11, 13)):
    """Create a SpecModel"""

    flux = np.full(npoints, flux)
    wavelength = np.arange(*wave_range, step=(wave_range[1] - wave_range[0]) / npoints)
    error = np.zeros(npoints)
    surf_bright = np.zeros(npoints)
    sb_error = np.zeros(npoints)
    dq = np.zeros(npoints)
    background = np.zeros(npoints)
    berror = np.zeros(npoints)
    npixels = np.zeros(npoints)

    data = [(i,j,k,l,m,n,o,p,q) for i,j,k,l,m,n,o,p,q in zip(wavelength, flux, error,
                                                             surf_bright, sb_error, dq,
                                                             background, berror, npixels)]

    spec_table = np.array(data, dtype=[('WAVELENGTH', 'f8'),('FLUX', 'f8'),
                                        ('ERROR', 'f8'), ('SURF_BRIGHT', 'f8'),
                                        ('SB_ERROR', 'f8'),('DQ', 'u4'), ('BACKGROUND', 'f8'),
                                        ('BERROR', 'f8'), ('NPIXELS', 'f8')])
    spec_model = datamodels.SpecModel()
    spec_model.spec_table = spec_table

    return spec_model
