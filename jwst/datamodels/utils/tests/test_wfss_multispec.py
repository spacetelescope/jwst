import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.io import fits
from numpy.testing import assert_allclose

from jwst.datamodels.utils.tests.wfss_helpers import (
    N_EXPOSURES,
    N_ROWS,
    N_SOURCES,
    wfss_comb,
    wfss_multi,
    wfss_spec2_multi,
    wfss_spec3_multi,
)
from jwst.datamodels.utils.wfss_multispec import (
    make_wfss_multiexposure,
    wfss_multiexposure_to_multispec,
)


@pytest.fixture
def wfss_spec2_multispec():
    return wfss_spec2_multi()


@pytest.fixture
def wfss_spec3_multispec():
    return wfss_spec3_multi()


@pytest.fixture
def wfss_multiexposure():
    return wfss_multi()


@pytest.fixture
def multi_combined():
    return wfss_comb()


@pytest.mark.parametrize("input_model_maker", ["wfss_spec2_multispec", "wfss_spec3_multispec"])
def test_make_wfss_multiexposure(input_model_maker, request):
    """
    Test reorganization of x1d data to flat file format.

    The two fixtures wfss_spec2_multispec and wfss_spec3_multispec have spectra that are
    identical, except that in the first case, all spectra are from the same exposure,
    and in the second case, all spectra are from the same source.
    make_wfss_multiexposure should be able to handle both cases.
    """
    input_model = request.getfixturevalue(input_model_maker)
    output_model = make_wfss_multiexposure(input_model)

    assert isinstance(output_model, dm.WFSSMultiSpecModel)
    if input_model_maker == "wfss_spec2_multispec":
        assert len(output_model.spec) == 1
        assert output_model.spec[0].spec_table.shape == (N_SOURCES,)
    elif input_model_maker == "wfss_spec3_multispec":
        assert len(output_model.spec) == 4
        assert output_model.spec[0].spec_table.shape == (1,)

    # check the required metadata attributes
    assert getattr(output_model.meta, "wcs", None) is None
    for i, exposure in enumerate(output_model.spec):
        assert exposure.group_id == str(i + 1)

    # check that units are present
    # test one vector-like column, which should come from the input specmodels
    # and one meta column, which should be copied from the schema by set_schema_units
    to_check = ["WAVELENGTH", "SOURCE_RA"]
    expected_units = ["um", "deg"]
    for exposure in output_model.spec:
        for col in to_check:
            assert col in exposure.spec_table.columns.names
            assert exposure.spec_table.columns[col].unit == expected_units[to_check.index(col)]


@pytest.mark.parametrize("input_model_maker", ["wfss_spec2_multispec", "wfss_spec3_multispec"])
def test_detector_meta_preserved(input_model_maker, request, tmp_path):
    """Preserve per-exposure detectors through WFSS conversion and FITS serialization."""
    input_model = request.getfixturevalue(input_model_maker)
    input_model.meta.instrument.detector = "NRCALONG"
    detectors = {}
    for spec in input_model.spec:
        # All sources in an exposure have the same detector, but exposures may differ.
        group_id = spec.meta.group_id
        detector = "NRCALONG" if int(group_id) % 2 else "NRCBLONG"
        spec.detector = detector
        detectors[group_id] = detector

    output_path = tmp_path / "wfss_x1d.fits"
    with make_wfss_multiexposure(input_model) as output_model:
        for exposure in output_model.spec:
            assert exposure.detector == detectors[exposure.group_id]
        output_model.save(output_path)

    with fits.open(output_path) as hdul:
        for hdu in hdul:
            if hdu.name == "EXTRACT1D":
                assert hdu.header["DETECTOR"] == detectors[hdu.header["EXPGRPID"]]

    with dm.open(output_path) as saved_model:
        for exposure in saved_model.spec:
            assert exposure.detector == detectors[exposure.group_id]
        expanded_models = wfss_multiexposure_to_multispec(saved_model)
        try:
            for model in expanded_models:
                for spec in model.spec:
                    assert spec.detector == detectors[spec.meta.group_id]
        finally:
            for model in expanded_models:
                model.close()


@pytest.mark.parametrize("flat_input", [False, True])
@pytest.mark.parametrize("reverse_inputs", [False, True])
@pytest.mark.parametrize("invalid_blong", [False, True])
def test_detectors_are_separated(wfss_spec3_multispec, flat_input, reverse_inputs, invalid_blong):
    """Keep both detectors for a source, even if one has no valid flux values."""
    inputs = []
    expected = {}
    for detector in ["NRCALONG", "NRCBLONG"]:
        model = wfss_spec3_multispec.copy()
        # The parent detector must not override the individual spectrum's detector.
        model.meta.instrument.detector = "NRCALONG"
        for spec in model.spec:
            spec.detector = detector
            spec.meta.filename = f"exposure_{spec.meta.group_id}_{detector}.fits"
            spec.spectral_order = 1
            if detector == "NRCBLONG":
                spec.spec_table["FLUX"] = np.nan if invalid_blong else 2.0
            key = (spec.meta.group_id, detector)
            expected[key] = (spec.meta.filename, spec.spec_table.copy())
        inputs.append(model)

    if flat_input:
        flat_models = [make_wfss_multiexposure(model) for model in inputs]
        for model in inputs:
            model.close()
        inputs = flat_models
    if reverse_inputs:
        inputs.reverse()

    def check_output(output):
        assert len(output.spec) == N_EXPOSURES * 2
        assert {(spec.group_id, spec.detector) for spec in output.spec} == set(expected)
        for spec in output.spec:
            filename, table = expected[(spec.group_id, spec.detector)]
            assert spec.filename == filename
            assert spec.spectral_order == 1
            assert len(spec.spec_table) == 1
            row = spec.spec_table[0]
            assert row["SOURCE_ID"] == wfss_spec3_multispec.spec[0].source_id
            assert row["N_ALONGDISP"] == len(table)
            assert len(row["FLUX"]) == len(table)
            for column in ["WAVELENGTH", "FLUX", "FLUX_ERROR", "DQ", "CONTAM", "CONTAMCORR_FLUX"]:
                assert_allclose(row[column], table[column])

    try:
        with make_wfss_multiexposure(inputs) as output:
            check_output(output)
            # Regrouping a mixed-detector product must also preserve both detectors.
            with make_wfss_multiexposure(output) as regrouped:
                check_output(regrouped)
    finally:
        for model in inputs:
            model.close()


def test_orders_are_separated(wfss_spec3_multispec):
    """Ensure that if there are multiple spectral orders, they end up in separate extensions."""

    multi_list = [wfss_spec3_multispec.copy() for _ in range(2)]
    # two input MultiSpecModel objects with the same source_id but different spectral orders
    for i, multi in enumerate(multi_list):
        for spec in multi.spec:
            spec.spectral_order = i + 1
            spec.meta.group_id = spec.meta.group_id + str(i)  # make group_id unique for each order

    output_model = make_wfss_multiexposure(multi_list)
    assert isinstance(output_model, dm.WFSSMultiSpecModel)
    assert len(output_model.spec) == 8  # 4 different exposures, 2 spectral orders each

    # ensure spectral order is in the metadata
    for i, exposure in enumerate(output_model.spec):
        assert exposure.spectral_order == (i // 4) + 1  # first 4 are order 1, next 4 are order 2
        assert exposure.group_id == str(i % 4 + 1) + str(exposure.spectral_order - 1)


def test_wfss_flat_to_multispec(wfss_multiexposure):
    """Test conversion of a WFSSMultiSpecModel back to a list of MultiSpecModel objects."""
    # first test that the fixture is giving a model with the correct dimensions
    assert isinstance(wfss_multiexposure, dm.WFSSMultiSpecModel)
    assert len(wfss_multiexposure.spec) == N_EXPOSURES
    assert wfss_multiexposure.spec[0].spec_table.shape == (N_SOURCES,)
    assert wfss_multiexposure.spec[0].dispersion_direction == 3

    # convert back to a list of MultiSpecModel objects
    multispec_list = wfss_multiexposure_to_multispec(wfss_multiexposure)

    # now test that the conversion to MultiSpecModel returns us to what we had before
    assert len(multispec_list) == N_SOURCES
    for i, multispec in enumerate(multispec_list):
        assert isinstance(multispec, dm.MultiSpecModel)
        assert len(multispec.spec) == N_EXPOSURES
        for j, spec in enumerate(multispec.spec):
            assert spec.source_id == i + 1  # they will now be sorted

            # check that the data is the same as the original
            assert_allclose(spec.spec_table["WAVELENGTH"], np.linspace(1.0, 10.0, N_ROWS))
            assert_allclose(spec.spec_table["FLUX"], np.ones(N_ROWS))

            # check meta that were explicitly input to spec
            assert spec.meta.filename == f"exposure_{j}.fits"
            assert spec.meta.group_id == str(j + 1)
            assert spec.dispersion_direction == 3

            # test that the rest of the metadata exist and are default values
            assert spec.source_type == "POINT"
            for name in [
                "source_ra",
                "source_dec",
                "extract2d_xstart",
                "extract2d_ystart",
                "extract2d_xstop",
                "extract2d_ystop",
            ]:
                assert getattr(spec, name) == 0.0


def test_wfss_multi_from_wfss_multi(wfss_multiexposure):
    """
    Test that a WFSSMultiSpecModel can be made from a list of WFSSMultiSpecModel objects.

    This situation is encountered when running the Spec3Pipeline,
    looping over extract_1d multiple times and compiling the results as a list.
    """
    # make a list of WFSSMultiSpecModel objects
    # and give them different sources
    inputs_list = []
    for i in range(2):
        this_source = wfss_multiexposure.copy()
        for exp in this_source.spec:
            exp.spec_table["SOURCE_ID"] = exp.spec_table["SOURCE_ID"] + N_SOURCES * i
        inputs_list.append(this_source)

    # combine the list into a single WFSSMultiSpecModel
    output_model = make_wfss_multiexposure(inputs_list)

    # check that the output model has the correct dimensions
    assert isinstance(output_model, dm.WFSSMultiSpecModel)
    assert len(output_model.spec) == N_EXPOSURES
    assert output_model.spec[0].spec_table.shape == (N_SOURCES * 2,)

    # test that the data has all the appropriate data and metadata
    for i, exposure in enumerate(output_model.spec):
        assert exposure.group_id == str(i + 1)
        assert exposure.spec_table.shape == (N_SOURCES * 2,)


@pytest.fixture
def comb1d_list(multi_combined):
    """
    Make a list of MultiCombinedSpecModel objects with N_SOURCES sources in list.

    Each MultiCombinedSpecModel object has only one spec, but the source_id is different.
    This looks like the output of calwebb_spec3 after calling combine_1d a bunch of times.
    """
    results_list = []
    for i in range(N_SOURCES):
        multi = multi_combined.copy()
        for spec in multi.spec:
            spec.source_id = N_SOURCES - i
        results_list.append(multi)
    return results_list
