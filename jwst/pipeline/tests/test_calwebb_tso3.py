import stdatamodels.jwst.datamodels as dm

from jwst.associations.asn_from_list import asn_from_list
from jwst.datamodels import ModelContainer
from jwst.extract_1d.tests.helpers import mock_niriss_soss_96_func, mock_niriss_soss_full_func
from jwst.pipeline.calwebb_tso3 import Tso3Pipeline


def niriss_soss_tso(subarray="SUBSTRIP96"):
    """
    Mock a NIRISS SOSS TSO calints model.

    Parameters
    ----------
    subarray : str
        May be either "SUBSTRIP96" (expected to complete processing)
        or "FULL" (expected to fail in extract_1d).

    Returns
    -------
    CubeModel
        An open model with just enough metadata to run through tso3.
    """
    if subarray == "FULL":
        input_model = mock_niriss_soss_full_func()
    else:
        input_model = mock_niriss_soss_96_func()
    input_model.meta.wcs = None
    input_model.meta.wcsinfo.s_region = "POLYGON ICRS 0 0 0 1 1 1 1 0"
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


def test_niriss_soss(tmp_path):
    """Smoke test for tso3 for a valid NIRISS SOSS TSO mode."""
    asn = tso3_asn(tmp_path, niriss_soss_tso())

    # Reduce runtime for soss extraction
    steps = {"extract_1d": {"soss_rtol": 0.1, "soss_tikfac": 2.434559775e-13}}

    # No errors
    Tso3Pipeline.call(asn, output_dir=str(tmp_path), steps=steps)

    # Check for expected output files
    expected = ["test_a3001_crfints.fits", "test_tso3_x1dints.fits", "test_tso3_whtlt.ecsv"]
    for filename in expected:
        assert (tmp_path / filename).exists()

    with dm.open(tmp_path / "test_tso3_x1dints.fits") as x1d:
        assert x1d.spec[0].s_region == "POLYGON ICRS 0 0 0 1 1 1 1 0"


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


def test_populate_tso_spectral_sregion(log_watcher):
    model = dm.TSOMultiSpecModel()
    model.spec.append(dm.SpecModel())
    model.spec.append(dm.SpecModel())
    cal_model = dm.CubeModel((10, 10, 10))
    cal_model_list = ModelContainer([cal_model, cal_model.copy()])

    # no s_region attributes present
    watcher = log_watcher(
        "jwst.pipeline.calwebb_tso3",
        message="No input model(s) have an `s_region` attribute; output S_REGION will not be set.",
        level="warning",
    )
    Tso3Pipeline()._populate_tso_spectral_sregion(model, cal_model_list)
    watcher.assert_seen()
    assert not model.spec[0].hasattr("s_region")

    # s_regions are all the same, should run without issues
    for m in cal_model_list:
        m.meta.wcsinfo.s_region = "POLYGON ICRS 0 0 0 1 1 1 1 0"
    Tso3Pipeline()._populate_tso_spectral_sregion(model, cal_model_list)
    assert model.spec[0].s_region == "POLYGON ICRS 0 0 0 1 1 1 1 0"
    assert not model.spec[1].hasattr("s_region")

    # s_regions differ, should warn and set to first
    cal_model_list[1].meta.wcsinfo.s_region = "POLYGON ICRS 1 1 1 2 2 2 2 1"
    watcher = log_watcher(
        "jwst.pipeline.calwebb_tso3",
        message="Input models have different S_REGION values;",
        level="warning",
    )
    Tso3Pipeline()._populate_tso_spectral_sregion(model, cal_model_list)
    watcher.assert_seen()
    assert model.spec[0].s_region == "POLYGON ICRS 0 0 0 1 1 1 1 0"
    assert not model.spec[1].hasattr("s_region")

    # only one model has s_region, should warn and set to that one
    del cal_model_list[0].meta.wcsinfo.s_region
    watcher = log_watcher(
        "jwst.pipeline.calwebb_tso3",
        message="One or more input model(s) are missing an `s_region` attribute;",
        level="warning",
    )
    Tso3Pipeline()._populate_tso_spectral_sregion(model, cal_model_list)
    watcher.assert_seen()
    assert model.spec[0].s_region == "POLYGON ICRS 1 1 1 2 2 2 2 1"
    assert not model.spec[1].hasattr("s_region")


def test_populate_tso_spectral_sregion_empty_container(log_watcher):
    """Test with an empty ModelContainer."""
    model = dm.TSOMultiSpecModel()
    model.spec.append(dm.SpecModel())
    cal_model_list = ModelContainer()

    # Should handle empty container gracefully
    watcher = log_watcher(
        "jwst.pipeline.calwebb_tso3", message="No input or output models provided;", level="warning"
    )
    Tso3Pipeline()._populate_tso_spectral_sregion(model, cal_model_list)
    watcher.assert_seen()
    assert not model.spec[0].hasattr("s_region")


def test_populate_tso_spectral_sregion_no_spec(log_watcher):
    """Test with a model that has no spec entries."""
    model = dm.TSOMultiSpecModel()
    cal_model = dm.CubeModel((10, 10, 10))
    cal_model.meta.wcsinfo.s_region = "POLYGON ICRS 0 0 0 1 1 1 1 0"
    cal_model_list = ModelContainer([cal_model])

    # Should handle gracefully when model.spec is empty
    watcher = log_watcher(
        "jwst.pipeline.calwebb_tso3", message="No input or output models provided;", level="warning"
    )
    Tso3Pipeline()._populate_tso_spectral_sregion(model, cal_model_list)
    watcher.assert_seen()
