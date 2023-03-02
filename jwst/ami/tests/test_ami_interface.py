import pytest

import stpipe

from stdatamodels.jwst import datamodels

from jwst.ami import AmiAnalyzeStep


def test_ami_analyze_cube_fail():
    """Make sure ami_analyze fails if input is CubeModel (_calints)"""
    model = datamodels.CubeModel((25, 19, 19))
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

    def mockreturn(input_model, reftype, observatory=None, asn_exptypes=None):
        return "N/A"
    monkeypatch.setattr(stpipe.crds_client, 'get_reference_file', mockreturn)

    with pytest.raises(RuntimeError):
        AmiAnalyzeStep.call(model)
