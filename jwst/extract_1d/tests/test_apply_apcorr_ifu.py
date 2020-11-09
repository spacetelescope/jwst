import pytest
import numpy as np
import os

from astropy.table import Table
from jwst.datamodels import IFUImageModel, NrsIfuApcorrModel, MirMrsApcorrModel
from jwst.extract_1d.apply_apcorr import ApCorr, ApCorrRadial, select_apcorr

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
NIR_TEST_FILES = {
    'IFU': os.path.join(data_dir, 'jwst_nirspec_apcorr_ifu_dummy.asdf')
}

MIRI_TEST_FILES = {
    'MRS': os.path.join(data_dir, 'jwst_miri_apcorr_ifu_dummy.asdf')
}

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
    dm.meta.instrument.channel = '1'
    dm.meta.instrument.band = 'SHORT'

    if instrument == 'MIRI':
        apcorr_model = MirMrsApcorrModel(MIRI_TEST_FILES)
        table = apcorr_model.apcorr_table

    if instrument == 'NIRSPEC':
        dm.meta.instrument.filter = 'CLEAR'
        dm.meta.instrument.grating = 'PRISM'
        apcorr_model = NrsIfuApcorrModel(NIR_TEST_FILES)
        table = apcorr_model.apcorr_table_clear_prism

    return dm, table

def inputs_miri():
# set up datamodel and apcorr table
    table = None

    instrument = 'MIRI'
    exptype = 'MIR_MRS'

    dm = IFUImageModel()
    dm.meta.instrument.name = instrument
    dm.meta.exposure.type = exptype
    dm.meta.instrument.channel = '1'
    dm.meta.instrument.band = 'SHORT'

    table_all = MirMrsApcorrModel(MIRI_TEST_FILES)
    table = table_all.apcorr_table
    return dm, table


def inputs_nirspec():
# set up datamodel and apcorr table
    table = None

    instrument = 'NIRSPEC'
    exptype = 'NRS_IFU'
    dm = IFUImageModel()
    dm.meta.instrument.name = instrument
    dm.meta.exposure.type = exptype
    dm.meta.instrument.filter = 'CLEAR'
    dm.meta.instrument.grating = 'PRISM'
    table = NrsIfuApcorrModel(NIR_TEST_FILES)
    return dm, table



def test_select_apcorr():
    
    dm,table = inputs_miri()
    apcorr_cls = select_apcorr(dm)

    if 'MRS' in dm.meta.exposure.type:
        assert apcorr_cls == ApCorrRadial

    dm,table = inputs_nirspec()
    apcorr_cls = select_apcorr(dm)
    if dm.meta.instrument.name == 'NIRSPEC':
        if dm.meta.exposure.type.upper() == 'NRS_IFU' :
            assert apcorr_cls == ApCorrRadial



def test_table_type():
    dm, table = inputs_miri()
    if dm.meta.instrument.name == 'MIRI':
        print(table.channel)
        assert table.channel == 'ANY'
        assert table.band == 'ANY'

    dm, table = inputs_nirspec()
    if dm.meta.instrument.name == 'NIRSPEC':
        assert table.filter == 'CLEAR'
        assert table.grating == 'PRISM'

#def test_table_wave(inputs):
#    dm, table = inputs
#    dummy_wave = np.zeros(100) + 0.5

#    if dm.meta.instrument.name == 'MIRI':
#        assert table.wavelength == dummy_wave

#    if dm.meta.instrument.name == 'NIRSPEC':
#        assert table.wavelength == dummy_wave

