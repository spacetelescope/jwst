"""
Test the engdblog step
"""

import numpy as np

from jwst.engdblog import EngDBLogStep
from jwst.lib.engdb_mast import EngdbMast


def _mock_init(*args, **kwargs):
    pass


def test_engdblogstep(caplog, tmp_path, monkeypatch):
    n_val = 10

    def _mock_values(*args, **kwargs):
        return [np.float64(0.0)] * n_val

    monkeypatch.setattr(EngdbMast, "__init__", _mock_init)
    monkeypatch.setattr(EngdbMast, "get_values", _mock_values)

    mnemonic = "INRSI_GWA_Y_TILT_AVGED"

    # Test normal input and bare string
    for arg in ([mnemonic], mnemonic):
        result = EngDBLogStep.call(arg)
        assert result == {mnemonic: np.float64(0.0)}
        assert "EngDBLogStep instance created" in caplog.text
        assert f"Step EngDBLogStep running with args (['{mnemonic}'],)" in caplog.text
        assert f"{mnemonic}[2022-01-25 02:00:00:2022-01-26 02:10:00] = " in caplog.text
        assert "Step EngDBLogStep done" in caplog.text

    # Test config file
    cfg_path = str(tmp_path / "engdblogstep.cfg")
    with open(cfg_path, "w") as cfg:
        cfg.write('verbosity = "all"\n')
    result = EngDBLogStep.call([mnemonic], config_file=cfg_path)
    assert len(result[mnemonic]) == n_val


def test_badmnemonic_and_novalues(caplog, monkeypatch):
    def _mock_values(*args, **kwargs):
        return []

    monkeypatch.setattr(EngdbMast, "__init__", _mock_init)
    monkeypatch.setattr(EngdbMast, "get_values", _mock_values)

    # Bad mnemonic
    mnemonic = "NOSUCHMNEMONIC"
    result = EngDBLogStep.call([mnemonic])
    assert result == dict()
    assert f"{mnemonic} has no entries in time range" in caplog.text

    # No values
    mnemonic = "INRSI_GWA_Y_TILT_AVGED"
    result = EngDBLogStep.call([mnemonic], etime="2016-01-02")
    assert result == dict()
    assert f"{mnemonic} has no entries in time range" in caplog.text


def test_multi_mnemonics(caplog, monkeypatch):
    # Actual DB values below but they are not checked, so no point to mock
    # carefully
    #
    # INRSI_GWA_Y_TILT_AVGED --> np.float64(0.0)
    # SA_ZATTEST1 --> np.float64(-0.4829944968)
    def _mock_values(*args, **kwargs):
        return [1.0]

    monkeypatch.setattr(EngdbMast, "__init__", _mock_init)
    monkeypatch.setattr(EngdbMast, "get_values", _mock_values)

    mnemonics = ["INRSI_GWA_Y_TILT_AVGED", "SA_ZATTEST1"]
    result = EngDBLogStep.call(mnemonics)
    assert len(result) == 2
    for mnemonic in mnemonics:
        assert f"{mnemonic}[2022-01-25 02:00:00:2022-01-26 02:10:00] = " in caplog.text
