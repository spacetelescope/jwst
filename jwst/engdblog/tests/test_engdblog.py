"""
Test the engdblog step
"""

import pytest

from jwst.lib import engdb_mast
from jwst.engdblog import EngDBLogStep


def test_engdblogstep(caplog, engdb):
    mnemonic = "INRSI_GWA_Y_TILT_AVGED"
    result = EngDBLogStep.call([mnemonic])
    assert isinstance(result, dict)
    assert mnemonic in result
    assert "EngDBLogStep instance created" in caplog.text
    assert mnemonic in caplog.text
    assert "Step EngDBLogStep running with args (['{}'],)".format(mnemonic) in caplog.text
    assert "{}[2022-01-25 02:00:00:2022-01-26 02:10:00] = ".format(mnemonic) in caplog.text
    assert "Step EngDBLogStep done" in caplog.text


def test_barestring(caplog, engdb):
    mnemonic = "INRSI_GWA_Y_TILT_AVGED"
    result = EngDBLogStep.call(mnemonic)
    assert isinstance(result, dict)
    assert mnemonic in result
    assert "EngDBLogStep instance created" in caplog.text
    assert mnemonic in caplog.text
    assert f"Step EngDBLogStep running with args ('{mnemonic}',)." in caplog.text
    assert "{}[2022-01-25 02:00:00:2022-01-26 02:10:00] = ".format(mnemonic) in caplog.text
    assert "Step EngDBLogStep done" in caplog.text


def test_badmnemonic(caplog, engdb):
    mnemonic = "NOSUCHMNEMONIC"
    result = EngDBLogStep.call([mnemonic])
    assert isinstance(result, dict)
    assert len(result) == 0
    assert "{} has no entries in time range".format(mnemonic) in caplog.text


def test_novalues(caplog, engdb):
    mnemonic = "INRSI_GWA_Y_TILT_AVGED"
    result = EngDBLogStep.call([mnemonic], etime="2016-01-02")
    assert isinstance(result, dict)
    assert len(result) == 0
    assert "{} has no entries in time range".format(mnemonic) in caplog.text


def test_all(caplog, engdb, tmp_path):
    mnemonic = "INRSI_GWA_Y_TILT_AVGED"
    cfg_file_name = "engdblogstep.cfg"
    cfg_path = str(tmp_path / cfg_file_name)
    with open(cfg_path, "w") as cfg:
        cfg.write('verbosity = "all"\n')
    result = EngDBLogStep.call([mnemonic], config_file=cfg_path)
    assert len(result[mnemonic]) > 1


def test_multi_mnemonics(caplog, engdb):
    mnemonics = ["INRSI_GWA_Y_TILT_AVGED", "SA_ZATTEST1"]
    result = EngDBLogStep.call(mnemonics)
    assert len(result) == 2
    for mnemonic in mnemonics:
        assert "{}[2022-01-25 02:00:00:2022-01-26 02:10:00] = ".format(mnemonic) in caplog.text


# #####################
# Utilities for testing
# #####################
@pytest.fixture
def engdb():
    yield engdb_mast.EngdbMast()
