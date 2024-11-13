import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.assign_wcs.util import wcs_bbox_from_shape
from jwst.datamodels import ModelContainer
from jwst.exp_to_source import multislit_to_container
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


@pytest.fixture
def simple_wcs_ifu():
    shape = (10, 50, 50)
    xcenter = shape[1] // 2.0

    def simple_wcs_function(x, y, z):
        """ Simple WCS for testing """
        crpix1 = xcenter
        crpix3 = 1.0
        cdelt1 = 0.1
        cdelt2 = 0.1
        cdelt3 = 0.01

        crval1 = 45.0
        crval2 = 45.0
        crval3 = 7.5

        wave = (z + 1 - crpix3) * cdelt3 + crval3
        ra = (z + 1 - crpix1) * cdelt1 + crval1
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


@pytest.fixture()
def mock_nirspec_bots(simple_wcs):
    model = dm.CubeModel()
    model.meta.instrument.name = 'NIRSPEC'
    model.meta.instrument.detector = 'NRS1'
    model.meta.instrument.filter = 'F290LP'
    model.meta.instrument.grating = 'G395H'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.instrument.fixed_slit = 'S1600A1'
    model.meta.exposure.type = 'NRS_BRIGHTOBJ'
    model.meta.subarray.name = 'SUB2048'
    model.meta.exposure.nints = 10

    model.name = 'S1600A1'
    model.meta.target.source_type = 'POINT'
    model.meta.wcsinfo.dispersion_direction = 1
    model.meta.wcs = simple_wcs

    model.data = np.arange(10 * 50 * 50, dtype=float).reshape((10, 50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05
    yield model
    model.close()


@pytest.fixture()
def mock_miri_ifu(simple_wcs_ifu):
    model = dm.IFUCubeModel()
    model.meta.instrument.name = 'MIRI'
    model.meta.instrument.detector = 'MIRIFULONG'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.exposure.type = 'MIR_MRS'

    model.meta.wcsinfo.dispersion_direction = 2
    model.meta.photometry.pixelarea_steradians = 1.0
    model.meta.wcs = simple_wcs_ifu

    model.data = np.arange(10 * 50 * 50, dtype=float).reshape((10, 50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05
    model.weightmap = np.full_like(model.data, 1.0)
    yield model
    model.close()


@pytest.fixture()
def mock_niriss_wfss_l3(mock_nirspec_fs_one_slit):
    model = dm.MultiSlitModel()
    model.meta.instrument.name = 'NIRISS'
    model.meta.instrument.detector = 'NIS'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.exposure.type = 'NIS_WFSS'

    nslit = 3
    for i in range(nslit):
        slit = mock_nirspec_fs_one_slit.copy()
        slit.name = str(i + 1)
        slit.meta.exposure.type = 'NIS_WFSS'
        model.slits.append(slit)

    container = multislit_to_container([model])['0']

    yield container
    container.close()


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


def test_extract_nirspec_bots(mock_nirspec_bots, simple_wcs):
    result = Extract1dStep.call(mock_nirspec_bots, apply_apcorr=False)
    assert result.meta.cal_step.extract_1d == 'COMPLETE'
    assert (result.spec[0].name == 'S1600A1')

    # output wavelength is the same as input
    _, _, expected_wave = simple_wcs(np.arange(50), np.arange(50))
    assert np.allclose(result.spec[0].spec_table['WAVELENGTH'], expected_wave)

    # output flux and errors are non-zero, exact values will depend
    # on extraction parameters
    assert np.all(result.spec[0].spec_table['FLUX'] > 0)
    assert np.all(result.spec[0].spec_table['FLUX_ERROR'] > 0)
    result.close()


@pytest.mark.parametrize('ifu_set_srctype', [None, 'EXTENDED'])
def test_extract_miri_ifu(mock_miri_ifu, simple_wcs_ifu, ifu_set_srctype):
    # Source type defaults to extended, results should be the
    # same with and without override
    result = Extract1dStep.call(mock_miri_ifu, ifu_covar_scale=1.0,
                                ifu_set_srctype=ifu_set_srctype)
    assert result.meta.cal_step.extract_1d == 'COMPLETE'

    # output wavelength is the same as input
    _, _, expected_wave = simple_wcs_ifu(np.arange(50), np.arange(50), np.arange(10))
    assert np.allclose(result.spec[0].spec_table['WAVELENGTH'], expected_wave)

    # output flux for extended data is a simple sum over all data
    # with a conversion factor for Jy
    assert np.allclose(result.spec[0].spec_table['FLUX'],
                       1e6 * np.sum(mock_miri_ifu.data, axis=(1, 2)))

    # output error comes from the sum of the variance components
    variance_sum = (np.sum(mock_miri_ifu.var_rnoise, axis=(1, 2))
                    + np.sum(mock_miri_ifu.var_poisson, axis=(1, 2))
                    + np.sum(mock_miri_ifu.var_flat, axis=(1, 2)))
    assert np.allclose(result.spec[0].spec_table['FLUX_ERROR'],
                       1e6 * np.sqrt(variance_sum))
    result.close()


def test_extract_container(mock_nirspec_mos, mock_nirspec_fs_one_slit, simple_wcs):
    # if not WFSS, the container is looped over, so contents need not match
    container = ModelContainer([mock_nirspec_mos, mock_nirspec_fs_one_slit])
    result = Extract1dStep.call(container)
    assert isinstance(result, ModelContainer)

    for model in result:
        for spec in model.spec:
            # output wavelength is the same as input
            _, _, expected_wave = simple_wcs(np.arange(50), np.arange(50))
            assert np.allclose(spec.spec_table['WAVELENGTH'], expected_wave)

            # output flux and errors are non-zero, exact values will depend
            # on extraction parameters
            assert np.all(spec.spec_table['FLUX'] > 0)
            assert np.all(spec.spec_table['FLUX_ERROR'] > 0)

    result.close()


def test_extract_niriss_wfss(mock_niriss_wfss_l3, simple_wcs):
    # input is a SourceModelContainer
    result = Extract1dStep.call(mock_niriss_wfss_l3)

    # output is a single spectral model (not a container)
    assert isinstance(result, dm.MultiSpecModel)
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
