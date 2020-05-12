import pytest
import numpy as np

from astropy.io import fits
from astropy.table import Table

from ...datamodels import DataModel

from ..apply_apcorr import ApCorr, ApCorrRadial, ApCorrPhase, select_apcorr


@pytest.fixture(
    scope='module',
    params=[
        ('MIRI', 'MIR_LRS-FIXEDSLIT'),
        ('MIRI', 'MIR_LRS-SLITLESS'),
        ('MIRI', 'MIR_MRS'),
        ('NIRCAM', 'NRC_GRISM'),
        ('NIRCAM', 'NRC_WFSS'),
        ('NIRISS', 'NIS_WFSS'),
        ('NIRSPEC', 'NRS_BRIGHTOBJ'),
        ('NIRSPEC', 'NRS_FIXEDSLIT'),
        ('NIRSPEC', 'NRS_IFU'),
        ('NIRSPEC', 'NRS_MSASPEC')
    ]
)
def inputs(request):
    table = None

    instrument, exptype = request.param

    dm = DataModel()
    dm.meta.instrument.name = instrument
    dm.meta.exposure.type = exptype

    if instrument == 'MIRI':
        if 'LRS' in exptype:
            dm.meta.subarray.name = 'FULL'
            table = Table(
                {
                    'subarray': ['FULL', 'sub02'],
                    'nelem_wl': [3, 3],
                    'nelem_size': [3, 3],
                    'wavelength': [[1, 2, 3, 0], [1, 2, 3, 0]],
                    'size': [[1, 2, 3, 0], [1, 2, 3, 0]],
                    'apcorr': np.full((2, 4, 4), 0.5)
                }
            )
        if 'MRS' in exptype:
            table = Table(
                {
                    'wavelength': [np.arange(12), ],
                    'radius': [np.arange(12), ],
                    'nelem_wl': [10, ],
                    'apcorr': np.full((1, 12), 0.5)
                }
            )

    if instrument == 'NIRSPEC':
        dm.meta.instrument.filter = 'CLEAR'
        dm.meta.instrument.grating = 'PRISM'
        table = Table(
            {
                'wavelength': [[1, 2, 3, 0], [1, 2, 3, 0]],
                'pixphase': [[0.5], [1.0]],
                'size': [[[1, 2, 3, 0], [1, 2, 3, 0]], [[1, 2, 3, 0], [1, 2, 3, 0]]],
                'nelem_size': [3, 3],
                'nelem_wl': [3, 3],
                'apcorr': np.full((2, 2, 4, 4), 0.5),
                'filter': ['CLEAR', 'filter02'],
                'grating': ['PRISM', 'grating02']
            }

        )
        if 'FIXEDSLIT' in dm.meta.exposure.type:
            dm.meta.instrument.fixed_slit = 'S200A1'
            table.add_column(['S200A1', 'slit2'], name='slit')

    if instrument == 'NIRCAM':
        dm.meta.instrument.filter = 'F322W2'
        dm.meta.instrument.pupil = 'GRISMR'

        table = Table(
            {
                'filter': ['F322W2', 'filter02'],
                'pupil': ['GRISMR', 'pupil02'],
                'nelem_wl': [3, 3],
                'nelem_size': [3, 3],
                'wavelength': [[1, 2, 3, 0], [1, 2, 3, 0]],
                'size': [[1, 2, 3, 0], [1, 2, 3, 0]],
                'apcorr': np.full((2, 4, 4), 0.5)
            }
        )

    if instrument == 'NIRISS':
        dm.meta.instrument.filter = 'GR150R'
        dm.meta.instrument.pupil = 'F090W'

        table = Table(
            {
                'filter': ['GR150R', 'filter02'],
                'pupil': ['F090W', 'pupil02'],
                'nelem_wl': [3, 3],
                'nelem_size': [3, 3],
                'wavelength': [[1, 2, 3, 0], [1, 2, 3, 0]],
                'size': [[1, 2, 3, 0], [1, 2, 3, 0]],
                'apcorr': np.full((2, 4, 4), 0.5)
            }
        )

    return dm, fits.table_to_hdu(table).data


@pytest.fixture(scope='module')
def apcorr_instance(inputs):
    dm, table = inputs
    apcorr = select_apcorr(dm)

    if dm.meta.instrument.name == 'NIRSPEC' and 'FIXEDSLIT' in dm.meta.exposure.type:
        return apcorr(dm, table, apcorr_row_units='pixels', slit='S200A1')

    return apcorr(dm, table, apcorr_row_units='pixels')


def test_select_apcorr(inputs):
    dm, _ = inputs

    apcorr_cls = select_apcorr(dm)

    if dm.meta.instrument.name == 'NIRSPEC':
        assert apcorr_cls == ApCorrPhase
    elif 'MRS' in dm.meta.exposure.type:
        assert apcorr_cls == ApCorrRadial
    else:
        assert apcorr_cls == ApCorr


class TestApCorr:

    def test_match_parameters(self, apcorr_instance):
        if apcorr_instance.model.meta.instrument.name == 'MIRI':
            if 'MRS' in apcorr_instance.model.meta.exposure.type:
                assert apcorr_instance.match_pars == {}
            else:
                assert apcorr_instance.match_pars == {'subarray': 'FULL'}

        if apcorr_instance.model.meta.instrument.name == 'NIRSPEC':
            if 'FIXEDSLIT' in apcorr_instance.model.meta.exposure.type:
                assert apcorr_instance.match_pars == {'filter': 'CLEAR', 'grating': 'PRISM', 'slit': 'S200A1'}
            else:
                assert apcorr_instance.match_pars == {'filter': 'CLEAR', 'grating': 'PRISM'}

        if apcorr_instance.model.meta.instrument.name == 'NIRCAM':
            assert apcorr_instance.match_pars == {'filter': 'F322W2', 'pupil': 'GRISMR'}

        if apcorr_instance.model.meta.instrument.name == 'NIRISS':
            assert apcorr_instance.match_pars == {'filter': 'GR150R', 'pupil': 'F090W'}

    def test_approximate(self, apcorr_instance):
        assert np.isclose(apcorr_instance.apcorr_func(0.5, 0.5)[0], 0.5)
