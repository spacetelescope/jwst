"""Test calwebb_spec2 for NIRSpec MSA"""

from os import path
from shutil import copy as file_copy

import pytest

from .helpers import (
    SCRIPT_DATA_PATH,
    abspath,
)

from ...associations import load_asn
from ...stpipe.step import Step

DATAPATH = abspath(
    '$TEST_BIGDATA/nirspec/test_datasets/msa/simulated-3nod'
)


@pytest.mark.bigdata
def test_msa_missing(mk_tmp_dirs, caplog):
    """Test MSA missing failure"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    input_file = 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod.fits'
    file_copy(
        path.join(
            DATAPATH,
            'level2a_twoslit',
            input_file
        ),
        tmp_data_path
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec2_basic.cfg'),
        path.join(tmp_data_path, input_file)
    ]

    with pytest.raises(Exception):
        Step.from_cmdline(args)

    assert 'Unable to open MSA FITS file (MSAMETFL)' in caplog.text


@pytest.mark.bigdata
def test_msa_missing_nofail(mk_tmp_dirs, caplog):
    """Test MSA missing failure"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    input_file = 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod.fits'
    file_copy(
        path.join(
            DATAPATH,
            'level2a_twoslit',
            input_file
        ),
        tmp_data_path
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec2_basic.cfg'),
        path.join(tmp_data_path, input_file),
        '--fail_on_exception=false'
    ]

    Step.from_cmdline(args)

    assert 'Unable to open MSA FITS file (MSAMETFL)' in caplog.text


@pytest.mark.bigdata
def test_msa_missing_skip(mk_tmp_dirs, caplog):
    """Test MSA missing failure"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    input_file = 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod.fits'
    file_copy(
        path.join(
            DATAPATH,
            'level2a_twoslit',
            input_file
        ),
        tmp_data_path
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec2_basic.cfg'),
        path.join(tmp_data_path, input_file),
        '--steps.assign_wcs.skip=true'
    ]

    Step.from_cmdline(args)

    assert 'Aborting remaining processing for this exposure.' in caplog.text


@pytest.mark.slow
@pytest.mark.bigdata
def test_run_msaflagging(mk_tmp_dirs, caplog):
    """Test msa flagging operation"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'jw95065006001_0_msa_twoslit.fits'),
        tmp_data_path
    )
    file_copy(
        path.join(DATAPATH, 'mos_udf_g235m_twoslit_spec2_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'mos_udf_g235m_twoslit_spec2_asn.json')
    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2a_twoslit', member['expname']),
                tmp_data_path
            )

    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec2_basic.cfg'),
        asn_path,
        '--steps.msa_flagging.skip=false'
    ]

    Step.from_cmdline(args)

    assert 'Step msa_flagging running with args' in caplog.text
    assert 'Step msa_flagging done' in caplog.text

    for product in asn['products']:
        prod_name = product['name'] + '_cal.fits'
        assert path.isfile(prod_name)
