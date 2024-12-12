import os

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.datamodels import ModelContainer
from jwst.extract_1d.extract_1d_step import Extract1dStep
from jwst.extract_1d.soss_extract import soss_extract


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


@pytest.mark.slow
def test_extract_niriss_soss_256(tmp_path, mock_niriss_soss_256):
    result = Extract1dStep.call(mock_niriss_soss_256, soss_rtol=0.1,
                                soss_modelname='soss_model.fits',
                                output_dir=str(tmp_path))
    assert result.meta.cal_step.extract_1d == 'COMPLETE'

    # output flux and errors are non-zero, exact values will depend
    # on extraction parameters
    assert np.all(result.spec[0].spec_table['FLUX'] > 0)
    assert np.all(result.spec[0].spec_table['FLUX_ERROR'] > 0)
    result.close()

    # soss output files are saved
    assert os.path.isfile(tmp_path / 'soss_model_SossExtractModel.fits')
    assert os.path.isfile(tmp_path / 'soss_model_AtocaSpectra.fits')


@pytest.mark.slow
def test_extract_niriss_soss_96(tmp_path, mock_niriss_soss_96):
    result = Extract1dStep.call(mock_niriss_soss_96, soss_rtol=0.1,
                                soss_modelname='soss_model.fits',
                                output_dir=str(tmp_path))
    assert result.meta.cal_step.extract_1d == 'COMPLETE'

    # output flux and errors are non-zero, exact values will depend
    # on extraction parameters
    assert np.all(result.spec[0].spec_table['FLUX'] > 0)
    assert np.all(result.spec[0].spec_table['FLUX_ERROR'] > 0)
    result.close()

    # soss output files are saved
    assert os.path.isfile(tmp_path / 'soss_model_SossExtractModel.fits')
    assert os.path.isfile(tmp_path / 'soss_model_AtocaSpectra.fits')


def test_extract_niriss_soss_fail(tmp_path, monkeypatch, mock_niriss_soss_96):
    # Mock an error condition in the soss extraction
    def mock(*args, **kwargs):
        return None, None, None
    monkeypatch.setattr(soss_extract, 'run_extract1d', mock)

    # None is returned
    result = Extract1dStep.call(mock_niriss_soss_96)
    assert result is None


def test_save_output_single(tmp_path, mock_nirspec_fs_one_slit):
    mock_nirspec_fs_one_slit.meta.filename = 'test_s2d.fits'
    result = Extract1dStep.call(mock_nirspec_fs_one_slit,
                                save_results=True, save_profile=True,
                                save_scene_model=True, output_dir=str(tmp_path),
                                suffix='x1d')

    output_path = str(tmp_path / 'test_x1d.fits')

    assert os.path.isfile(output_path)
    assert os.path.isfile(output_path.replace('x1d', 'profile'))
    assert os.path.isfile(output_path.replace('x1d', 'scene_model'))

    result.close()


def test_save_output_multislit(tmp_path, mock_nirspec_mos):
    mock_nirspec_mos.meta.filename = 'test_s2d.fits'
    result = Extract1dStep.call(mock_nirspec_mos,
                                save_results=True, save_profile=True,
                                save_scene_model=True, output_dir=str(tmp_path),
                                suffix='x1d')

    output_path = str(tmp_path / 'test_x1d.fits')

    assert os.path.isfile(output_path)

    # intermediate files for multislit data contain the slit name
    for slit in mock_nirspec_mos.slits:
        assert os.path.isfile(output_path.replace('x1d', f'{slit.name}_profile'))
        assert os.path.isfile(output_path.replace('x1d', f'{slit.name}_scene_model'))

    result.close()
