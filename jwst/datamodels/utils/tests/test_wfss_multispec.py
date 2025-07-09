import pytest
import numpy as np
from numpy.testing import assert_allclose

import stdatamodels.jwst.datamodels as dm
from jwst.datamodels.utils.wfss_multispec import (
    make_wfss_multiexposure,
    make_wfss_multicombined,
    wfss_multiexposure_to_multispec,
)

from jwst.datamodels.utils.tests.wfss_helpers import (
    wfss_spec2_multi,
    wfss_spec3_multi,
    wfss_multi,
    wfss_comb,
    N_EXPOSURES,
    N_SOURCES,
    N_ROWS,
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
    assert not hasattr(output_model.meta, "wcs")
    for i, exposure in enumerate(output_model.spec):
        assert exposure.group_id == str(i + 1)

    # check that units are present
    # test one vector-like column, which should come from the input specmodels
    # and one meta column, which should be copied from the schema by set_schema_units
    to_check = ["WAVELENGTH", "SOURCE_RA"]
    expected_units = ["um", "degrees"]
    for exposure in output_model.spec:
        for col in to_check:
            assert col in exposure.spec_table.columns.names
            assert exposure.spec_table.columns[col].unit == expected_units[to_check.index(col)]


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


def test_make_wfss_combined(comb1d_list):
    """Test restructuring of a list of combine_1d output models into a WFSSMultiCombinedSpecModel."""
    output_model = make_wfss_multicombined(comb1d_list)
    assert isinstance(output_model, dm.WFSSMultiCombinedSpecModel)
    assert len(output_model.spec) == 2  # 2 spectral orders

    for i, spec in enumerate(output_model.spec):
        assert spec.spec_table.shape == (N_SOURCES,)

        # test units
        to_check = ["WAVELENGTH", "SOURCE_RA"]
        expected_units = ["um", "degrees"]
        for col in to_check:
            assert col in spec.spec_table.columns.names
            assert spec.spec_table.columns[col].unit == expected_units[to_check.index(col)]

        # check metadata
        assert spec.dispersion_direction == 3
        assert spec.spectral_order == i + 1
