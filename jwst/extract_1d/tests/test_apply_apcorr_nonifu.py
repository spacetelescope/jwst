import pytest
import numpy as np
import os

from astropy.io import fits
from astropy.table import Table

from stdatamodels.jwst.datamodels import JwstDataModel

from jwst.extract_1d.apply_apcorr import ApCorr, ApCorrPhase, select_apcorr

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
NIR_TEST_FILES = {
    'MSASPEC': os.path.join(data_dir, 'jwst_nirspec_apcorr_msa_dummy.fits'),
    'FIXEDSLIT/BRIGHTOBJ': os.path.join(data_dir, 'jwst_nirspec_apcorr_fs_dummy.fits')
}


@pytest.fixture(
    scope='module',
    params=[
        ('MIRI', 'MIR_LRS-FIXEDSLIT'),
        ('MIRI', 'MIR_LRS-SLITLESS'),
        ('NIRCAM', 'NRC_GRISM'),
        ('NIRCAM', 'NRC_WFSS'),
        ('NIRISS', 'NIS_WFSS'),
        ('NIRSPEC', 'NRS_BRIGHTOBJ'),
        ('NIRSPEC', 'NRS_FIXEDSLIT'),
        ('NIRSPEC', 'NRS_MSASPEC')
    ]
)
def inputs(request):
    table = None

    instrument, exptype = request.param

    dm = JwstDataModel()
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

    if instrument == 'NIRSPEC':  # Too complicated to come up with silly example data; using "dummy" ref file
        dm.meta.instrument.filter = 'CLEAR'
        dm.meta.instrument.grating = 'PRISM'

        if 'FIXEDSLIT' in dm.meta.exposure.type:
            dm.meta.instrument.fixed_slit = 'S200A1'
            table = Table.read(NIR_TEST_FILES['FIXEDSLIT/BRIGHTOBJ'], format='fits')[:2]
        else:
            table = Table.read(NIR_TEST_FILES['MSASPEC'], format='fits')[:2]

        table['APCORR'] = table['APCORR'].reshape((len(table), 3, 2048, 3))  # Reshape test data to expected shape

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
        return apcorr(dm, table, 'pixels', slit='S200A1')

    return apcorr(dm, table, 'pixels')


def test_select_apcorr(inputs):
    dm, _ = inputs
    apcorr_cls = select_apcorr(dm)

    if dm.meta.instrument.name == 'NIRSPEC':
        assert apcorr_cls == ApCorrPhase
    else:
        assert apcorr_cls == ApCorr


class TestApCorr:

    def test_match_parameters(self, apcorr_instance):
        if apcorr_instance.model.meta.instrument.name == 'MIRI':
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
        if isinstance(apcorr_instance, ApCorrPhase):
            assert np.isclose(
                apcorr_instance.apcorr_func(
                    apcorr_instance.reference['wavelength'][0],
                    apcorr_instance.reference['size'][0],
                    0.5
                )[0],
                1
            )
        else:
            assert np.isclose(apcorr_instance.apcorr_func(0.5, 0.5)[0], 0.5)
