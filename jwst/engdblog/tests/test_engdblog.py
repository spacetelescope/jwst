"""
Test the engdblog step
"""
import os
import pytest
from tempfile import TemporaryDirectory

from jwst.engdblog import EngDBLogStep
from jwst.lib.tests.engdb_mock import EngDB_Mocker


def test_engdblogstep(caplog, engdb):
    mnemonic = 'INRSI_GWA_Y_TILT_AVGED'
    result = EngDBLogStep.call([mnemonic])
    assert isinstance(result, dict)
    assert mnemonic in result
    assert 'EngDBLogStep instance created' in caplog.text
    assert mnemonic in caplog.text
    assert "Step EngDBLogStep running with args (['{}'],)".format(mnemonic) in caplog.text
    assert '{}[2016-01-01:2016-01-31] = '.format(mnemonic) in caplog.text
    assert 'Step EngDBLogStep done' in caplog.text


@pytest.mark.xfail(
    reason='See Jira JP-1108',
    run=False
)
def test_barestring(caplog, engdb):
    mnemonic = 'INRSI_GWA_Y_TILT_AVGED'
    result = EngDBLogStep.call(mnemonic)
    assert isinstance(result, dict)
    assert mnemonic in result
    assert 'EngDBLogStep instance created' in caplog.text
    assert mnemonic in caplog.text
    assert f"Step EngDBLogStep running with args ('{mnemonic}')." in caplog.text
    assert '{}[2016-01-01:2016-01-31] = '.format(mnemonic) in caplog.text
    assert 'Step EngDBLogStep done' in caplog.text


def test_badmnemonic(caplog, engdb):
    mnemonic = 'NOSUCHMNEMONIC'
    result = EngDBLogStep.call([mnemonic])
    assert isinstance(result, dict)
    assert len(result) == 0
    assert 'Cannot retrieve info for {}'.format(mnemonic) in caplog.text


def test_novalues(caplog, engdb):
    mnemonic = 'INRSI_GWA_Y_TILT_AVGED'
    result = EngDBLogStep.call([mnemonic], etime='2016-01-02')
    assert isinstance(result, dict)
    assert len(result) == 0
    assert '{} has no entries in time range'.format(mnemonic) in caplog.text


def test_all(caplog, engdb):
    mnemonic = 'INRSI_GWA_Y_TILT_AVGED'
    cfg_file_name = 'engdblogste.cfg'
    with TemporaryDirectory() as cfg_dir:
        cfg_path = os.path.join(cfg_dir, cfg_file_name)
        with open(cfg_path, 'w') as cfg:
            cfg.write('verbosity = "all"\n')
        result = EngDBLogStep.call([mnemonic], config_file=cfg_path)
        assert len(result[mnemonic]) > 1


def test_multi_mnemonics(caplog, engdb):
    mnemonics = ['INRSI_GWA_Y_TILT_AVGED', 'SA_ZATTEST1']
    result = EngDBLogStep.call(mnemonics)
    assert len(result) == 2
    for mnemonic in mnemonics:
        assert '{}[2016-01-01:2016-01-31] = '.format(mnemonic) in caplog.text



# #####################
# Utilities for testing
# #####################
@pytest.fixture
def engdb():
    with EngDB_Mocker() as mocker: # noqa: F841
        yield
