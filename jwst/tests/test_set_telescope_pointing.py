"""
Test suite for set_telescope_pointing

Notes
-----
This file has been named specifically so it is not
automatically found by py.test. This is because, to test,
a connection to the internal engineering service is needed,
which is generally not available.
"""

from __future__ import absolute_import

from astropy.io import fits
from astropy.time import Time
import copy
import numpy as np
import os
import pytest
import requests_mock
import sys
from tempfile import TemporaryDirectory

from jwst.lib import engdb_tools

sys.path.insert(
    0,
    os.path.join(os.path.dirname(__file__), '../../scripts')
)

import set_telescope_pointing as stp


# Setup mock engineering service
GOOD_MNEMONIC = 'INRSI_GWA_Y_TILT_AVGED'
STARTTIME = Time('2014-01-01')
ENDTIME = Time('2014-01-04')

# DB population
POINTING_MNEMONICS = {
    'SA_ZATTEST1':  -3.691528614E-01,
    'SA_ZATTEST2':   3.376328232E-01,
    'SA_ZATTEST3':   5.758532642E-02,
    'SA_ZATTEST4':   8.639526444E-01,
    'SA_ZRFGS2J11': -1.004440003E-03,
    'SA_ZRFGS2J21':  3.381458357E-03,
    'SA_ZRFGS2J31':  9.999937784E-01,
    'SA_ZRFGS2J12':  9.999994955E-01,
    'SA_ZRFGS2J22': -3.900000000E-14,
    'SA_ZRFGS2J32':  1.004445746E-03,
    'SA_ZRFGS2J13':  3.396491461E-06,
    'SA_ZRFGS2J23':  9.999942829E-01,
    'SA_ZRFGS2J33': -3.381456651E-03,
    'SA_ZADUCMDX':   0.000000000E+00,
    'SA_ZADUCMDY':   0.000000000E+00,
}

# Header defaults
TARG_RA = 345.0
TARG_DEC = -87.0
V2_REF = 200.0
V3_REF = -350.0
V3I_YANG = 42.0
VPARITY = -1


def register_responses(mocker, mnemonics, starttime, endtime):
    request_url = ''.join([
        engdb_tools.ENGDB_BASE_URL,
        'Data/',
        '{mnemonic}',
        '?sTime={starttime}',
        '&eTime={endtime}'
    ])

    starttime_mil = int(starttime.unix * 1000)
    endtime_mil = int(endtime.unix * 1000)
    response_generic = {
        'AllPoints': 1,
        'Count': 2,
        'ReqSTime': '/Date({:013d}+0000)/'.format(starttime_mil),
        'ReqETime': '/Date({:013d}+0000)/'.format(endtime_mil),
        'TlmMnemonic': None,
        'Data': [
            {
                'EUValue': 0.1968553,
                'ObsTime': '/Date({:013d}+0000)/'.format(starttime_mil)
            },
            {
                'EUValue': 0.1968553,
                'ObsTime': '/Date({:013d}+0000)/'.format(endtime_mil)
            }
        ],
    }

    responses = {}
    for mnemonic, value in mnemonics.items():
        response = copy.deepcopy(response_generic)
        response['TlmMnemonic'] = mnemonic
        response['Data'][0]['EUValue'] = value
        response['Data'][1]['EUValue'] = value
        mocker.get(
            request_url.format(
                mnemonic=mnemonic,
                starttime=starttime.iso,
                endtime=endtime.iso
            ),
            json=response
        )
        responses[mnemonic] = response

    return responses


@pytest.fixture
def eng_db():
    with requests_mock.Mocker() as rm:

        # Define response for aliveness
        url = ''.join([
            engdb_tools.ENGDB_BASE_URL,
            engdb_tools.ENGDB_METADATA
        ])
        rm.get(url, text='Success')

        # Define required DB responses.
        responses = register_responses(
            rm,
            POINTING_MNEMONICS,
            STARTTIME,
            ENDTIME
        )
        edb = engdb_tools.ENGDB_Service()
        for mnemonic in POINTING_MNEMONICS:
            r = edb.get_records(mnemonic, STARTTIME, ENDTIME)
            assert r == responses[mnemonic]
        yield edb


@pytest.fixture
def fits_file():
    hdu = fits.PrimaryHDU()
    header = hdu.header
    header['EXPSTART'] = STARTTIME.mjd
    header['EXPEND'] = ENDTIME.mjd
    header['TARG_RA'] = TARG_RA
    header['TARG_DEC'] = TARG_DEC
    header['V2_REF'] = V2_REF
    header['V3_REF'] = V3_REF
    header['V3I_YANG'] = V3I_YANG
    header['VPARITY'] = VPARITY
    hdul = fits.HDUList([hdu])

    with TemporaryDirectory() as path:
        file_path = os.path.join(path, 'fits.fits')
        hdul.writeto(file_path)
        yield file_path


def test_get_pointing_fail():
    with pytest.raises(Exception):
        q, j2fgs_matrix, fmscorr = stp.get_pointing(47892.0, 48256.0)


def test_get_pointing(eng_db):
        assert stp.get_pointing(STARTTIME, ENDTIME)


def test_add_wcs_default(fits_file):
    try:
        stp.add_wcs(fits_file)
    except:
        pytest.skip('Live ENGDB service is not accessible.')

    hdul = fits.open(fits_file)
    header = hdul[0].header
    assert header['RA_V1'] == TARG_RA
    assert header['DEC_V1'] == TARG_DEC
    assert header['PA_V3'] == 0.
    assert header['CRVAL1'] == TARG_RA
    assert header['CRVAL2'] == TARG_DEC
    assert header['PC1_1'] == 1.0
    assert header['PC1_2'] == 0.0
    assert header['PC2_1'] == 0.0
    assert header['PC2_2'] == 1.0
    assert header['RA_REF'] == TARG_RA
    assert header['DEC_REF'] == TARG_DEC
    assert np.isclose(header['ROLL_REF'], 0.07993869)
    assert header['WCSAXES'] == 0.


def test_add_wcs_with_db(eng_db, fits_file):
        stp.add_wcs(fits_file)

        hdul = fits.open(fits_file)
        header = hdul[0].header
        assert np.isclose(header['RA_V1'], 348.9278669)
        assert np.isclose(header['DEC_V1'], -38.749239)
        assert np.isclose(header['PA_V3'], 50.1767077)
        assert np.isclose(header['CRVAL1'], 348.8776709)
        assert np.isclose(header['CRVAL2'], -38.854159)
        assert np.isclose(header['PC1_1'], -0.0385309)
        assert np.isclose(header['PC1_2'], -0.9992574)
        assert np.isclose(header['PC2_1'], 0.9992574)
        assert np.isclose(header['PC2_2'], -0.0385309)
        assert np.isclose(header['RA_REF'], 348.8776709)
        assert np.isclose(header['DEC_REF'], -38.854159)
        assert np.isclose(header['ROLL_REF'], 354.76818)
        assert header['WCSAXES'] == 0.
