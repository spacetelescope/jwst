import pytest

from jwst.associations.asn_from_list import asn_from_list
from jwst.extract_1d.tests.helpers import mock_niriss_soss_96_func
from jwst.pipeline.calwebb_tso3 import Tso3Pipeline


def niriss_soss_tso(subarray=None):
    """
    Mock a NIRISS SOSS TSO calints model.

    Returns
    -------
    CubeModel
        An open model with just enough metadata to run through tso3.
    """
    input_model = mock_niriss_soss_96_func()
    if subarray is not None:
        input_model.meta.subarray.name = subarray
    input_model.meta.wcs = None
    input_model.meta.visit.tsovisit = True
    input_model.int_times = input_model.int_times
    return input_model


def tso3_asn(tmp_path, input_model):
    input_name = "test_calints.fits"
    input_model.meta.filename = input_name
    input_model.save(str(tmp_path / input_name))
    input_model.close()

    asn_path = tmp_path / "test_tso3_asn.json"
    asn = asn_from_list([input_name], product_name="test_tso3")

    with asn_path.open("w") as outfile:
        name, serialized = asn.dump(format="json")
        outfile.write(serialized)

    return str(asn_path)


@pytest.mark.slow
def test_niriss_soss(tmp_path):
    """Smoke test for tso3 for a valid NIRISS SOSS TSO mode."""
    asn = tso3_asn(tmp_path, niriss_soss_tso())

    # Reduce runtime for soss extraction
    steps = {"extract_1d": {"soss_rtol": 0.1}}

    # No errors
    Tso3Pipeline.call(asn, output_dir=str(tmp_path), steps=steps)

    # Check for expected output files
    expected = ["test_a3001_crfints.fits", "test_tso3_x1dints.fits", "test_tso3_whtlt.ecsv"]
    for filename in expected:
        assert (tmp_path / filename).exists()


def test_niriss_soss_full(tmp_path):
    """Smoke test for tso3 for an invalid NIRISS SOSS TSO mode."""
    asn = tso3_asn(tmp_path, niriss_soss_tso(subarray="FULL"))

    # No errors
    Tso3Pipeline.call(asn, output_dir=str(tmp_path))

    # Check for expected output files
    expected = ["test_a3001_crfints.fits"]
    not_expected = ["test_tso3_x1dints.fits", "test_tso3_whtlt.ecsv"]
    for filename in expected:
        assert (tmp_path / filename).exists()
    for filename in not_expected:
        assert not (tmp_path / filename).exists()
