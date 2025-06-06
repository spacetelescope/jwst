import pytest
import numpy as np
from numpy.testing import assert_allclose
from astropy.io import fits
import stdatamodels.jwst.datamodels as dm
from jwst.datamodels.utils.wfss_multispec import make_wfss_multiexposure, make_wfss_multicombined, wfss_multiexposure_to_multispec


N_SOURCES = 5
N_EXPOSURES = 4
N_ROWS = 3

@pytest.fixture
def example_spec():

    def mock_wcs(*args, **kwargs):
        return 0.0, 0.0, 0.0
    spec = dm.SpecModel()
    spectable_dtype = spec.schema["properties"]["spec_table"]["datatype"]
    recarray_dtype = [(d["name"], d["datatype"]) for d in spectable_dtype]
    spec.meta.wcs = mock_wcs
    spec_table = np.recarray((N_ROWS,), dtype=recarray_dtype)
    spec_table["WAVELENGTH"] = np.linspace(1.0, 10.0, N_ROWS)
    spec_table["FLUX"] = np.ones(N_ROWS)
    spec.spec_table = spec_table
    spec.spec_table.columns["wavelength"].unit = "um"
    return spec

def _add_multispec_meta(spec):
    """
    Add specmeta attributes to a spec-like ObjectNode inside a MultiSpecModel.
    
    This only includes attributes that do NOT need to change for each source/exposure
    for the purposes of this test.

    input spec is updated in place.
    """
    spec.source_type = "POINT"
    spec.source_ra = 0.0
    spec.source_dec = 0.0
    spec.extract2d_xstart = 0.0
    spec.extract2d_ystart = 0.0
    spec.extract2d_xstop = 0.0
    spec.extract2d_ystop = 0.0
    spec.extraction_xstart = 0.0
    spec.extraction_ystart = 0.0
    spec.extraction_xstop = 0.0
    spec.extraction_ystop = 0.0


@pytest.fixture
def wfss_spec2_multispec(example_spec):
    """
    Set up a MultiSpecModel object that looks like outputs from extract_1d during calwebb_spec2.
    
    Each of the spectra is from the SAME exposure, but DIFFERENT sources.
    """
    multi = dm.MultiSpecModel()
    for i in range(N_SOURCES):
        # create a new SpecModel for each source
        spec = example_spec.copy()
        spec.meta.filename = f"exposure_1.fits" # all sources in the same exposure
        spec.meta.group_id = "1" # all sources in the same exposure
        spec.source_id = N_SOURCES - i # reverse the order to test sorting
        spec.name = str(spec.source_id)
        _add_multispec_meta(spec)
        multi.spec.append(spec)

    return multi


@pytest.fixture
def wfss_spec3_multispec(example_spec):
    """
    Set up a MultiSpecModel object that looks like outputs from extract_1d during calwebb_spec3.
    
    Each of the spectra is from the SAME source, but DIFFERENT exposures.
    """
    multi = dm.MultiSpecModel()
    multi.meta.exposure.exposure_time = 7.0
    multi.meta.exposure.integration_time = 7.0
    for j in range(N_EXPOSURES):
        # create a new SpecModel for each exposure
        spec = example_spec.copy()
        spec.meta.filename = f"exposure_{j}.fits"
        spec.meta.group_id = str(j + 1)
        spec.source_id = 999 # all sources are the same
        spec.dispersion_direction = 3
        spec.name = str(spec.source_id)
        _add_multispec_meta(spec)
        multi.spec.append(spec)

    return multi


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

    assert isinstance(output_model, dm.WFSSMultiExposureSpecModel)
    if input_model_maker == "wfss_spec2_multispec":
        assert len(output_model.exposures) == 1
        assert output_model.exposures[0].spec_table.shape == (N_SOURCES,)
    elif input_model_maker == "wfss_spec3_multispec":
        assert len(output_model.exposures) == 4
        assert output_model.exposures[0].spec_table.shape == (1,)
    
    # check the required metadata attributes
    assert not hasattr(output_model.meta, "wcs")
    for i, exposure in enumerate(output_model.exposures):
        assert exposure.group_id == str(i + 1)
    
    # check that units are present
    # test one vector-like column, which should come from the input specmodels
    # and one meta column, which should be copied from the schema by set_schema_units
    to_check = ["WAVELENGTH", "SOURCE_RA"]
    expected_units = ["um", "degrees"]
    for exposure in output_model.exposures:
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
    assert isinstance(output_model, dm.WFSSMultiExposureSpecModel)
    assert len(output_model.exposures) == 8  # 4 different exposures, 2 spectral orders each

    # ensure spectral order is in the metadata
    for i, exposure in enumerate(output_model.exposures):
        assert exposure.spectral_order == (i // 4) + 1  # first 4 are order 1, next 4 are order 2
        assert exposure.group_id == str(i%4 + 1) + str(exposure.spectral_order - 1)


@pytest.fixture
def wfss_multiexposure(wfss_spec3_multispec):
    """Make a MultiExposureSpecModel object with N_EXPOSURES exposures and N_SOURCES sources."""
    inputs_list = []
    for i in range(N_SOURCES):
        this_source = wfss_spec3_multispec.copy()
        for spec in this_source.spec:
            spec.source_id = N_SOURCES - i
        inputs_list.append(this_source)
    output_model = make_wfss_multiexposure(inputs_list)
    return output_model


def test_wfss_flat_to_multispec(wfss_multiexposure):
    """Test conversion of a WFSSMultiExposureSpecModel back to a list of MultiSpecModel objects."""
    # first test that the fixture is giving a model with the correct dimensions
    assert isinstance(wfss_multiexposure, dm.WFSSMultiExposureSpecModel)
    assert len(wfss_multiexposure.exposures) == N_EXPOSURES
    assert wfss_multiexposure.exposures[0].spec_table.shape == (N_SOURCES,)
    assert wfss_multiexposure.exposures[0].dispersion_direction == 3

    # convert back to a list of MultiSpecModel objects
    multispec_list = wfss_multiexposure_to_multispec(wfss_multiexposure)

    # now test that the conversion to MultiSpecModel returns us to what we had before
    assert len(multispec_list) == N_SOURCES
    for i, multispec in enumerate(multispec_list):
        assert isinstance(multispec, dm.MultiSpecModel)
        assert len(multispec.spec) == N_EXPOSURES
        for j, spec in enumerate(multispec.spec):
            assert spec.source_id == i + 1 # they will now be sorted

            # check that the data is the same as the original
            assert_allclose(spec.spec_table["WAVELENGTH"], np.linspace(1.0, 10.0, N_ROWS))
            assert_allclose(spec.spec_table["FLUX"], np.ones(N_ROWS))

            # check meta that were explicitly input to spec
            assert spec.meta.filename == f"exposure_{j}.fits"
            assert spec.meta.group_id == str(j + 1)
            assert spec.dispersion_direction == 3

            # test that the rest of the metadata exist and are default values
            assert spec.source_type == "POINT"
            for name in ["source_ra", "source_dec", "extract2d_xstart", "extract2d_ystart",
                         "extract2d_xstop", "extract2d_ystop"]:
                assert getattr(spec, name) == 0.0


def test_wfss_multi_from_wfss_multi(wfss_multiexposure):
    """
    Test that a WFSSMultiExposureSpecModel can be made from a list of WFSSMultiExposureSpecModel objects.

    This situation is encountered when running the Spec3Pipeline,
    looping over extract_1d multiple times and compiling the results as a list.
    """
    # make a list of WFSSMultiExposureSpecModel objects
    # and give them different sources
    inputs_list = []
    for i in range(2):
        this_source = wfss_multiexposure.copy()
        for exp in this_source.exposures:
            exp.spec_table["SOURCE_ID"] = exp.spec_table["SOURCE_ID"] + N_SOURCES * i
        inputs_list.append(this_source)

    # combine the list into a single WFSSMultiExposureSpecModel
    output_model = make_wfss_multiexposure(inputs_list)

    # check that the output model has the correct dimensions
    assert isinstance(output_model, dm.WFSSMultiExposureSpecModel)
    assert len(output_model.exposures) == N_EXPOSURES
    assert output_model.exposures[0].spec_table.shape == (N_SOURCES*2,)

    # test that the data has all the appropriate data and metadata
    for i, exposure in enumerate(output_model.exposures):
        assert exposure.group_id == str(i + 1)
        assert exposure.spec_table.shape == (N_SOURCES*2,)
    

@pytest.fixture
def multi_combined():
    """
    Make a MultiCombinedSpecModel object with N_SOURCES sources, N_ROWS rows, and 2 orders.
    
    This looks like the output of combine_1d.
    """
    multi = dm.MultiCombinedSpecModel()
    spec = dm.CombinedSpecModel()
    _add_multispec_meta(spec)
    spectable_dtype = spec.schema["properties"]["spec_table"]["datatype"]
    recarray_dtype = [(d["name"], d["datatype"]) for d in spectable_dtype]
    spec_table = np.recarray((N_ROWS,), dtype=recarray_dtype)
    spec_table["WAVELENGTH"] = np.linspace(1.0, 10.0, N_ROWS)
    spec_table["FLUX"] = np.ones(N_ROWS)
    spec.spec_table = spec_table
    spec.spec_table.columns["wavelength"].unit = "um"
    spec.dispersion_direction = 3

    spec.spectral_order = 1
    spec2 = spec.copy()
    spec2.spectral_order = 2
    multi.spec.append(spec)
    multi.spec.append(spec2)
    return multi


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
