"""Unit tests for ami_analyze module and step."""

import pytest
import numpy as np
import math
import stpipe
from astropy.io import fits

from jwst.ami import AmiAnalyzeStep


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


def test_ami_analyze_step(example_model, tmp_cwd):
    model, _, _ = AmiAnalyzeStep.call(example_model)

    # ensure the primary header contains expected wcs info, but no SCI extension
    model.save("ami_analyze_step_output.fits")
    model.close()
    with fits.open("ami_analyze_step_output.fits") as hdul:
        extensions = [ext.name for ext in hdul]
        assert "SCI" not in extensions
        for kw in ["VPARITY", "V3I_YANG", "ROLL_REF"]:
            assert kw in list(hdul[0].header.keys())


def test_ami_analyze_step_no_affine(example_model):
    AmiAnalyzeStep.call(example_model, affine2d=None)


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
