"""Test utilities"""
import copy
import pytest

import numpy as np

from stdatamodels.jwst import datamodels
from stdatamodels.properties import ObjectNode

from jwst.lib import pipe_utils
from jwst.associations.lib import dms_base
from jwst.lib.pipe_utils import is_tso

all_exp_types = dms_base.IMAGE2_NONSCIENCE_EXP_TYPES + \
    dms_base.IMAGE2_SCIENCE_EXP_TYPES + \
    dms_base.SPEC2_SCIENCE_EXP_TYPES

exp_types = [
    (exp_type, False)
    for exp_type in all_exp_types
    if exp_type not in dms_base.TSO_EXP_TYPES
]
exp_types.extend([
    (exp_type, True)
    for exp_type in dms_base.TSO_EXP_TYPES
])


@pytest.mark.parametrize(
    'exp_type, expected',
    exp_types
)
def test_is_tso_from_exptype(exp_type, expected):
    """Test is_tso integrity based on exp_type"""
    model = datamodels.JwstDataModel()
    model.meta.exposure.type = exp_type.upper()
    assert pipe_utils.is_tso(model) is expected


@pytest.mark.parametrize(
    'tsovisit, expected',
    [
        (True, True),
        (False, False),
    ]
)
def test_is_tso_from_tsoflag(tsovisit, expected):
    """Test is_tso integrity based on the TSO flag"""
    model = datamodels.JwstDataModel()
    model.meta.visit.tsovisit = tsovisit
    assert pipe_utils.is_tso(model) is expected


def test_is_tso_nrcgrism_nints1():
    """Test is_tso with NRC_TSGRISM and NINTS=1"""
    model = datamodels.JwstDataModel()
    model.meta.exposure.type = "NRC_TSGRISM"
    model.meta.visit.tsovisit = False
    model.meta.exposure.nints = 10

    # with NINTS>1, should be True
    assert pipe_utils.is_tso(model)

    # with NINTS=1, should be False
    model.meta.exposure.nints = 1
    assert not pipe_utils.is_tso(model)

    # with hardwired TSO EXP_TYPE's, should always be True
    assert (is_tso(model) or model.meta.exposure.type.lower() in
            ['nrc_tsimage', 'nrc_tsgrism'])


def test_is_irs2_1():
    """Test is_irs2 using a numpy array"""
    x = np.ones((2048, 2), dtype=np.float32)
    assert not pipe_utils.is_irs2(x)
    x = np.ones((3200, 2), dtype=np.float32)
    assert pipe_utils.is_irs2(x)


def test_is_irs2_2():
    """Test is_irs2 using a jwst data model"""
    x = np.ones((2048, 2), dtype=np.float32)
    model = datamodels.ImageModel(data=x)
    assert not pipe_utils.is_irs2(model)
    x = np.ones((3200, 2), dtype=np.float32)
    model = datamodels.ImageModel(data=x)
    assert pipe_utils.is_irs2(model)


@pytest.mark.parametrize('shape', [(10,), (10, 10), (10, 10, 10)])
@pytest.mark.parametrize('extensions',
                         [('data',), ('data', 'dq'),
                          ('data', 'err', 'var_rnoise'),
                          ('data', 'err', 'var_rnoise', 'var_poisson',
                           'var_flat', 'dq')])
def test_match_nans_and_flags(shape, extensions):
    # Set up a model with matching data shapes, variable extensions
    dnu = datamodels.dqflags.pixel['DO_NOT_USE']
    model = datamodels.SlitModel()
    invalid = np.full(shape, False)
    for i, extname in enumerate(extensions):
        if extname == 'dq':
            data = np.zeros(shape, dtype=np.uint32)
            data.flat[i] = dnu
        else:
            data = np.ones(shape)
            data.flat[i] = np.nan
        setattr(model, extname, data)
        invalid.flat[i] = True

    # Match flags and NaNs across all extensions
    pipe_utils.match_nans_and_flags(model)

    # Check that all extensions have all invalid flags
    for extname in extensions:
        data = getattr(model, extname)
        if extname == 'dq':
            assert np.all(data[invalid] == dnu)
            assert np.all(data[~invalid] == 0)
        else:
            assert np.all(np.isnan(data[invalid]))
            assert np.all(data[~invalid] == 1)

    model.close()


def test_match_nans_and_flags_dq_only():
    # Set up a model with only a DQ array, to test edge conditions
    dnu = datamodels.dqflags.pixel['DO_NOT_USE']
    model = datamodels.MaskModel()
    shape = (10, 10)
    invalid = np.full(shape, False)

    model.dq = np.zeros(shape, dtype=np.uint32)
    model.dq.flat[5] = dnu
    invalid.flat[5] = True

    # Match flags and NaNs across all extensions
    dq_copy = model.dq.copy()
    pipe_utils.match_nans_and_flags(model)

    # Check that DQ extension is unchanged
    assert np.all(model.dq == dq_copy)

    model.close()


def test_match_nans_and_flags_no_data(caplog):
    # Set up a model with no data, to test edge conditions
    model = datamodels.JwstDataModel()
    model_copy = model.copy()
    pipe_utils.match_nans_and_flags(model)

    # nothing happens
    assert model == model_copy
    assert 'WARNING' not in caplog.text
    assert 'ERROR' not in caplog.text

    model.close()
    model_copy.close()


def test_match_nans_and_flags_not_a_model():
    # Make sure the function throws an error if input is not a model
    model = "not_a_model"
    with pytest.raises(TypeError, match="not a datamodel"):
        pipe_utils.match_nans_and_flags(model)

    # empty JwstDataModel is okay
    model = datamodels.JwstDataModel()
    pipe_utils.match_nans_and_flags(model)

    # a datamodel object node is also okay
    multislit = datamodels.MultiSlitModel()
    multislit.slits.append(datamodels.SlitModel())
    assert not isinstance(multislit.slits[0], datamodels.JwstDataModel)
    assert isinstance(multislit.slits[0], ObjectNode)
    pipe_utils.match_nans_and_flags(multislit.slits[0])


def test_match_nans_and_flags_shape_mismatch():
    # Set up a model with mismatched data shapes
    dnu = datamodels.dqflags.pixel['DO_NOT_USE']
    model = datamodels.SlitModel()

    shape = (10, 10)
    bad_shape = (10, 5)
    invalid = np.full(shape, False)
    extensions = ['data', 'err', 'dq', 'var_poisson']
    for i, extname in enumerate(extensions):
        if extname == 'dq':
            data = np.zeros(bad_shape, dtype=np.uint32)
            data.flat[i] = dnu
        elif extname == 'var_poisson':
            data = np.ones(bad_shape, dtype=np.uint32)
            data.flat[i] = dnu
        else:
            data = np.ones(shape)
            data.flat[i] = np.nan
            invalid.flat[i] = True
        setattr(model, extname, data)

    # Match flags and NaNs across all extensions:
    # will warn for data mismatch
    model_copy = model.copy()
    pipe_utils.match_nans_and_flags(model)

    # Only data and err are updated
    for extname in extensions:
        data = getattr(model, extname)
        if extname == 'data' or extname == 'err':
            assert np.all(np.isnan(data[invalid]))
            assert np.all(data[~invalid] == 1)
        else:
            assert np.allclose(data, getattr(model_copy, extname))

    model.close()
    model_copy.close()


def test_determine_vector_and_meta_columns():

    input_schema = datamodels.SpecModel().schema
    output_schema = datamodels.WFSSMultiSpecModel().schema

    vector_columns, meta_columns = pipe_utils.determine_vector_and_meta_columns(
        input_schema, output_schema
    )

    # Check that the names ended up in the right place
    input_names = [s["name"] for s in output_schema["properties"]["spec_table"]["datatype"]]
    output_names = [s["name"] for s in output_schema["properties"]["spec_table"]["datatype"]]
    vector_names = [s[0].upper() for s in vector_columns]
    meta_names = [s[0].upper() for s in meta_columns]

    assert np.all(np.isin(vector_names, input_names))
    assert set(vector_names + meta_names) == set(output_names)
    assert set(vector_names).isdisjoint(set(meta_names))
    
    # test that the datatypes are all numpy types
    for _name, dtype in vector_columns:
        assert np.issubdtype(dtype, np.number)
    for _name, dtype in meta_columns:
        assert (np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.bytes_))


def test_make_empty_recarray():
    n_rows = 10
    n_sources = 5
    vector_columns = [('FLUX', np.float32), ('WAVELENGTH', np.float64)]
    meta_columns = [('NAME', np.str_)]
    
    recarray = pipe_utils.make_empty_recarray(n_rows, n_sources, vector_columns, meta_columns)
    
    # Check the shapes
    assert recarray.shape == (n_sources,)
    assert recarray["FLUX"].shape == (n_sources, n_rows)
    assert recarray["WAVELENGTH"].shape == (n_sources, n_rows)
    assert recarray["NAME"].shape == (n_sources,)
    
    # Check the data types of the columns
    for name, dtype in vector_columns + meta_columns:
        assert recarray[name].dtype == dtype


@pytest.fixture
def empty_recarray():
    n_rows = 10
    n_sources = 5
    vector_columns = [('FLUX', np.float32), ('WAVELENGTH', np.float64)]
    meta_columns = [('NAME', "S20"), ('SOURCE_ID', np.int32)]
    return pipe_utils.make_empty_recarray(n_rows, n_sources, vector_columns, meta_columns)


@pytest.mark.parametrize("ignore_columns", [["SOURCE_ID"], ["FLUX"], []])
def test_populate_recarray(empty_recarray, ignore_columns):
    n_rows = 10
    n_sources = 5
    vector_columns = [('FLUX', np.float32), ('WAVELENGTH', np.float64)]
    meta_columns = [('NAME', "S20"), ('SOURCE_ID', np.int32)]
    output_table = copy.deepcopy(empty_recarray)

    for i in range(n_sources):
        input_spec = datamodels.SpecModel()
        input_spec.name = f"Source {i}"
        input_spec.source_id = i

        this_output = output_table[i]

        # hack input model schema to have only the flux and wavelength columns
        input_spec.schema["properties"]["spec_table"]["datatype"] = \
            input_spec.schema["properties"]["spec_table"]["datatype"][:2]

        # give all the input spectra different numbers of rows
        input_spec.spec_table = np.ones((n_rows-i,), dtype=[(name, dtype) for name, dtype in vector_columns])
        pipe_utils.populate_recarray(
            this_output,
            input_spec,
            n_rows,
            vector_columns,
            meta_columns,
            ignore_columns=ignore_columns
        )
        
    # Check that the output table has been populated correctly
    for i in range(n_sources):
        for name, _ in vector_columns:
            if name not in ignore_columns:
                # ensure proper nan-padding of output
                expected = np.ones((n_rows,))*np.nan
                expected[:n_rows-i] = 1.0
                assert np.array_equal(output_table[name][i], expected, equal_nan=True)
            else:
                # These should all be zero because that's the initial value from the schema
                assert np.allclose(output_table[name][i], 0.0, equal_nan=False, atol=1e-10)
                pass
        assert output_table["NAME"][i] == np.bytes_(f"Source {i}")
        if "SOURCE_ID" not in ignore_columns:
            assert output_table["SOURCE_ID"][i] == i
        else:
            # If SOURCE_ID is ignored, it should be set to 0 because that's the schema default value
            assert output_table["SOURCE_ID"][i] == 0
