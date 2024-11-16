import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.assign_wcs.util import wcs_bbox_from_shape
from jwst.extract_1d.extract_1d_step import Extract1dStep


@pytest.fixture
def simple_wcs():
    shape = (50, 50)
    xcenter = shape[1] // 2.0

    def simple_wcs_function(x, y):
        """ Simple WCS for testing """
        crpix1 = xcenter
        crpix3 = 1.0
        cdelt1 = 0.1
        cdelt2 = 0.1
        cdelt3 = 0.01

        crval1 = 45.0
        crval2 = 45.0
        crval3 = 7.5

        wave = (x + 1 - crpix3) * cdelt3 + crval3
        ra = (x + 1 - crpix1) * cdelt1 + crval1
        dec = np.full_like(ra, crval2 + 1 * cdelt2)

        return ra, dec, wave

    simple_wcs_function.bounding_box = wcs_bbox_from_shape(shape)

    return simple_wcs_function


@pytest.fixture()
def mock_nirspec_fs_one_slit(simple_wcs):
    model = dm.SlitModel()
    model.meta.instrument.name = 'NIRSPEC'
    model.meta.instrument.detector = 'NRS1'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.instrument.fixed_slit = 'S200A1'
    model.meta.exposure.type = 'NRS_FIXEDSLIT'
    model.meta.subarray.name = 'ALLSLITS'

    model.source_type = 'EXTENDED'
    model.meta.wcsinfo.dispersion_direction = 1
    model.meta.wcs = simple_wcs

    model.data = np.arange(50 * 50, dtype=float).reshape((50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05
    yield model
    model.close()


@pytest.fixture()
def mock_nirspec_mos(mock_nirspec_fs_one_slit):
    model = dm.MultiSlitModel()
    model.meta.instrument.name = 'NIRSPEC'
    model.meta.instrument.detector = 'NRS1'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.exposure.type = 'NRS_MSASPEC'

    nslit = 3
    for i in range(nslit):
        slit = mock_nirspec_fs_one_slit.copy()
        slit.name = str(i + 1)
        model.slits.append(slit)

    yield model
    model.close()


@pytest.mark.parametrize('slit_name', [None, 'S200A1', 'S1600A1'])
def test_extract_nirspec_fs_slit(mock_nirspec_fs_one_slit, simple_wcs, slit_name):
    mock_nirspec_fs_one_slit.name = slit_name
    result = Extract1dStep.call(mock_nirspec_fs_one_slit)
    assert result.meta.cal_step.extract_1d == 'COMPLETE'

    if slit_name is None:
        assert (result.spec[0].name
                == mock_nirspec_fs_one_slit.meta.instrument.fixed_slit)
    else:
        assert result.spec[0].name == slit_name

    # output wavelength is the same as input
    _, _, expected_wave = simple_wcs(np.arange(50), np.arange(50))
    assert np.allclose(result.spec[0].spec_table['WAVELENGTH'], expected_wave)

    # output flux and errors are non-zero, exact values will depend
    # on extraction parameters
    assert np.all(result.spec[0].spec_table['FLUX'] > 0)
    assert np.all(result.spec[0].spec_table['FLUX_ERROR'] > 0)
    result.close()


def test_extract_nirspec_mos_multi_slit(mock_nirspec_mos, simple_wcs):
    result = Extract1dStep.call(mock_nirspec_mos)
    assert result.meta.cal_step.extract_1d == 'COMPLETE'

    for i, spec in enumerate(result.spec):
        assert spec.name == str(i + 1)

        # output wavelength is the same as input
        _, _, expected_wave = simple_wcs(np.arange(50), np.arange(50))
        assert np.allclose(spec.spec_table['WAVELENGTH'], expected_wave)

        # output flux and errors are non-zero, exact values will depend
        # on extraction parameters
        assert np.all(spec.spec_table['FLUX'] > 0)
        assert np.all(spec.spec_table['FLUX_ERROR'] > 0)

    result.close()
