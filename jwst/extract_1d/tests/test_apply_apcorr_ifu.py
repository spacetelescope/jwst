import pytest
import numpy as np
import os

from jwst.datamodels import IFUImageModel, NrsIfuApcorrModel, MirMrsApcorrModel
from jwst.extract_1d.apply_apcorr import ApCorrRadial, select_apcorr

NIR_TEST_FILES = 'data/jwst_nirspec_apcorr_ifu_dummy.asdf'

MIRI_TEST_FILES =  'data/jwst_miri_apcorr_ifu_dummy.asdf'

@pytest.fixture(
# define the tests/mode we want to run
    scope='module',
    params=[
        ('MIRI', 'MIR_MRS'),
        ('NIRSPEC', 'NRS_IFU'),
    ]
)

def inputs(request):
# set up datamodel and apcorr table
    table = None

    instrument, exptype = request.param

    dm = IFUImageModel()
    dm.meta.instrument.name = instrument
    dm.meta.exposure.type = exptype

    if instrument == 'MIRI':
        apcorr_model = MirMrsApcorrModel(MIRI_TEST_FILES)
        table = apcorr_model.apcorr_table
        apcorr_model.close()

    if instrument == 'NIRSPEC':
        apcorr_model = NrsIfuApcorrModel(NIR_TEST_FILES)
        # return the clear_prism table
        table = apcorr_model.apcorr_table_clear_prism
        apcorr_model.close()

    return dm, table



def test_select_apcorr(inputs):
    dm,table = inputs
    apcorr_cls = select_apcorr(dm)

    if 'MRS' in dm.meta.exposure.type:
        assert apcorr_cls == ApCorrRadial


    if dm.meta.instrument.name == 'NIRSPEC':
        assert apcorr_cls == ApCorrRadial


def test_table_type(inputs):
    dm, table = inputs
    dummy_wave = np.zeros(100) + 0.5

    if dm.meta.instrument.name == 'MIRI':
        assert table.channel == 'ANY'
        assert table.band == 'ANY'
        assert np.all(table.wavelength == dummy_wave)

    if dm.meta.instrument.name == 'NIRSPEC':
        assert table.filter == 'clear'
        assert table.grating == 'prism'
        assert np.all(table.wavelength == dummy_wave)



    instrument = 'NIRSPEC'
    exptype = 'NRS_IFU'
    dm = IFUImageModel()
    dm.meta.instrument.name = instrument
    dm.meta.exposure.type = exptype
    dm.meta.instrument.filter = 'CLEAR'
    dm.meta.instrument.grating = 'PRISM'
    apcorr_model = NrsIfuApcorrModel(NIR_TEST_FILES)
    apcorr_model.close()
    assert apcorr_model.apcorr_table_clear_prism.filter == 'clear'
    assert apcorr_model.apcorr_table_clear_prism.grating == 'prism'
