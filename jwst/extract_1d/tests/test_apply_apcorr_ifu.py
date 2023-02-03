import pytest
import numpy as np
from asdf import AsdfFile

from stdatamodels.jwst.datamodels import IFUCubeModel, NrsIfuApcorrModel, MirMrsApcorrModel

from jwst.extract_1d.apply_apcorr import ApCorrRadial, select_apcorr


@pytest.fixture(scope='module')
def dummy_nirspec_ref(tmpdir_factory):
    """ Generate a dummy apcorr ref file """
    filename = tmpdir_factory.mktemp('dummy_apcorr')
    filename = str(filename.join('dummy_nirspec_apcorr.asdf'))

    refap = {}
    refap['meta'] = {}
    refap['apcorr_table'] = {}
    refap['meta']['telescope'] = 'JWST'
    refap['meta']['reftype'] = 'APCORR'
    refap['meta']['exp_type'] = 'P_EXP_TY'
    refap['meta']['detector'] = 'N/A'
    refap['meta']['datamodel'] = 'NrsIfuApcorrModel'
    refap['meta']['version'] = '1.0'
    refap['meta']['name'] = 'NIRSPEC'
    refap['meta']['origin'] = 'STScI'

    dummy_wave = np.zeros(100) + 0.5
    dummy_radius = np.zeros((3, 100)) + 0.5
    dummy_apcorr = np.ones((3, 100))
    dummy_apcorr_error = np.zeros((3, 100))

    refap['apcorr_table'] = {}
    refap['apcorr_table']['wavelength'] = dummy_wave.copy()
    refap['apcorr_table']['radius'] = dummy_radius.copy()
    refap['apcorr_table']['apcorr'] = dummy_apcorr.copy()
    refap['apcorr_table']['apcorr_err'] = dummy_apcorr_error.copy()
    refap['apcorr_table']['wavelength_units'] = 'microns'
    refap['apcorr_table']['radius_units'] = 'arcseconds'
    refap['apcorr_table']['filter'] = 'ANY'
    refap['apcorr_table']['grating'] = 'ANY'

    ff = AsdfFile(refap)
    ff.set_array_storage(refap['apcorr_table']['wavelength'], 'inline')
    ff.set_array_storage(refap['apcorr_table']['radius'], 'inline')
    ff.set_array_storage(refap['apcorr_table']['apcorr'], 'inline')
    ff.set_array_storage(refap['apcorr_table']['apcorr_err'], 'inline')

    ff.write_to(filename)
    return filename


@pytest.fixture(scope='module')
def dummy_miri_ref(tmpdir_factory):
    """ Generate a dummy apcorr ref file """
    filename = tmpdir_factory.mktemp('dummy_apcorr')
    filename = str(filename.join('dummy_miri_apcorr.asdf'))

    refap = {}
    refap['meta'] = {}
    refap['apcorr_table'] = {}
    refap['meta']['telescope'] = 'JWST'
    refap['meta']['reftype'] = 'APCORR'
    refap['meta']['exp_type'] = 'MIRI_MRS'
    refap['meta']['detector'] = 'N/A'
    refap['meta']['datamodel'] = 'MirMrsApcorrModel'
    refap['meta']['version'] = '1.0'
    refap['meta']['name'] = 'MIRI'
    refap['meta']['origin'] = 'STScI'

    dummy_wave = np.zeros(100) + 0.5
    dummy_radius = np.zeros((3, 100)) + 0.5
    dummy_apcorr = np.ones((3, 100))
    dummy_apcorr_error = np.zeros((3, 100))

    refap['apcorr_table'] = {}
    refap['apcorr_table']['wavelength'] = dummy_wave.copy()
    refap['apcorr_table']['radius'] = dummy_radius.copy()
    refap['apcorr_table']['apcorr'] = dummy_apcorr.copy()
    refap['apcorr_table']['apcorr_err'] = dummy_apcorr_error.copy()
    refap['apcorr_table']['wavelength_units'] = 'microns'
    refap['apcorr_table']['radius_units'] = 'arcseconds'
    refap['apcorr_table']['channel'] = 'ANY'
    refap['apcorr_table']['band'] = 'ANY'

    ff = AsdfFile(refap)
    ff.set_array_storage(refap['apcorr_table']['wavelength'], 'inline')
    ff.set_array_storage(refap['apcorr_table']['radius'], 'inline')
    ff.set_array_storage(refap['apcorr_table']['apcorr'], 'inline')
    ff.set_array_storage(refap['apcorr_table']['apcorr_err'], 'inline')

    ff.write_to(filename)
    return filename


@pytest.fixture
def miri_cube():
    model = IFUCubeModel((15, 10, 10))
    instrument = 'MIRI'
    exptype = 'MIR_MRS'
    model.meta.instrument.name = instrument
    model.meta.exposure.type = exptype
    model.meta.instrument.channel = '1'
    model.meta.instrument.band = 'SHORT'
    return model


@pytest.fixture
def nirspec_cube():
    model = IFUCubeModel((15, 10, 10))
    instrument = 'NIRSPEC'
    exptype = 'NRS_IFU'
    model.meta.instrument.name = instrument
    model.meta.exposure.type = exptype
    model.meta.instrument.filter = 'CLEAR'
    model.meta.instrument.prism = 'PRISM'
    return model


def test_select_apcorr_miri(miri_cube):

    apcorr_cls = select_apcorr(miri_cube)
    assert apcorr_cls == ApCorrRadial


def test_select_apcorr_nirspec(nirspec_cube):

    apcorr_cls = select_apcorr(nirspec_cube)
    assert apcorr_cls == ApCorrRadial


def test_table_type_miri(_jail, dummy_miri_ref, miri_cube):
    dummy_wave = np.zeros(100) + 0.5
    with MirMrsApcorrModel(dummy_miri_ref) as apcorr_model:
        table = apcorr_model.apcorr_table

        assert table.channel == 'ANY'
        assert table.band == 'ANY'
        assert np.all(table.wavelength == dummy_wave)


def test_table_type_nirspec(_jail, dummy_nirspec_ref, nirspec_cube):
    dummy_wave = np.zeros(100) + 0.5
    with NrsIfuApcorrModel(dummy_nirspec_ref) as apcorr_model:
        table = apcorr_model.apcorr_table

        assert table.filter == 'ANY'
        assert table.grating == 'ANY'
        assert np.all(table.wavelength == dummy_wave)
