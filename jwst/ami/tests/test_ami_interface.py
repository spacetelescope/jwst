import pytest

from jwst import datamodels
from jwst.ami import AmiAnalyzeStep
import jwst


def test_ami_analyze_calints_fail(tmpdir):
    """Make sure ami_analyze fails if input is CubeModel (_calints)"""
    input_file = str(tmpdir.join("ami_analyze_input.fits"))
    model = datamodels.CubeModel((25, 19, 19))
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.filter = "F277W"
    model.meta.observation.date = "2019-01-01"
    model.meta.observation.time = "00:00:00"
    model.save(input_file)
    with pytest.raises(RuntimeError):
        AmiAnalyzeStep.call(input_file)


def test_ami_analyze_cube_fail():
    """Make sure ami_analyze fails if input is CubeModel (_calints)"""
    model = datamodels.ImageModel((25, 19, 19))
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.filter = "F277W"
    model.meta.observation.date = "2019-01-01"
    model.meta.observation.time = "00:00:00"
    with pytest.raises(RuntimeError):
        AmiAnalyzeStep.call(model)


def test_ami_analyze_no_reffile_fail(monkeypatch):
    """Make sure that ami_analyze fails if no throughput reffile is available"""
    model = datamodels.ImageModel((19, 19))
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.filter = "F277W"
    model.meta.observation.date = "2019-01-01"
    model.meta.observation.time = "00:00:00"

    def mockreturn(input_model, reftype, observatory=None):
        return("N/A")
    monkeypatch.setattr(jwst.stpipe.crds_client, 'get_reference_file', mockreturn)

    with pytest.raises(RuntimeError):
        AmiAnalyzeStep.call(model)
