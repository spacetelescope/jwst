"""Unit tests for ami_analyze module and step."""
import pytest
import numpy as np
import math
import stpipe
from stdatamodels.jwst import datamodels

from jwst.ami import AmiAnalyzeStep


@pytest.fixture()
def example_model():
    model = datamodels.CubeModel((2, 80, 80))
    # some non-zero data is required as this step will center
    # the image and find the centroid (both fail with all zeros)
    model.data[:, 24, 24] = 1
    model.data[:, 28, 28] = 1
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.filter = "F277W"
    model.meta.subarray.name = "SUB80"
    model.meta.observation.date = "2021-12-26"
    model.meta.observation.time = "00:00:00"
    model.meta.target.proposer_name = ""
    model.meta.program.pi_name = "someone"
    model.meta.target.catalog_name = ""
    model.meta.visit.start_time = "2022-06-05 12:15:41.5020000"
    model.meta.wcsinfo.roll_ref = 171.8779402866089
    model.meta.wcsinfo.v3yangle = 0.56126717
    model.meta.filename = "test_calints.fits"
    model.meta.instrument.pupil = "NRM"
    model.meta.exposure.type = "NIS_AMI"
    return model


@pytest.mark.parametrize("oversample", [2, 4])
def test_ami_analyze_even_oversample_fail(example_model, oversample):
    """Make sure ami_analyze fails if oversample is even"""
    with pytest.raises(ValueError, match="Oversample value must be an odd integer."):
        AmiAnalyzeStep.call(example_model, oversample=oversample)


def test_ami_analyze_no_reffile_fail(monkeypatch, example_model):
    """Make sure that ami_analyze fails if no throughput reffile is available"""

    def mockreturn(input_model, reftype, observatory=None, asn_exptypes=None):
        return "N/A"

    monkeypatch.setattr(stpipe.crds_client, "get_reference_file", mockreturn)

    with pytest.raises(RuntimeError, match="No THROUGHPUT reference file found."):
        AmiAnalyzeStep.call(example_model)


def test_ami_analyze_step(example_model):
    AmiAnalyzeStep.call(example_model)


def create_throughput(nelem):
    """Create a symmetric dummy throughput function that has values near
    0 on the wings and near 1 at the center.
    """
    ctr = int(nelem / 2.0)

    lower_half = [2.0 / (1.0 + math.e ** (-5.0 * i / ctr)) - 1.0 for i in range(ctr)]

    throughput = np.zeros(nelem, dtype=np.float32)
    throughput[:ctr] = lower_half
    throughput[ctr:] = lower_half[::-1]  # mirror image for upper half

    return throughput
