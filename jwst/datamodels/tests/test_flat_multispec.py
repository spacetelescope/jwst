import pytest
import logging
import copy
import numpy as np

import stdatamodels.jwst.datamodels as dm
from jwst.datamodels.flat_multispec import (determine_vector_and_meta_columns, make_empty_recarray, populate_recarray)
from jwst.tests.helpers import LogWatcher


def test_determine_vector_and_meta_columns():

    input_schema = dm.MRSSpecModel().schema
    in_cols = input_schema["properties"]["spec_table"]["datatype"]
    out_cols = in_cols.copy()
    out_cols.append({"name": "SOURCE_ID", "datatype": "int32"})
    out_cols.append({"name": "NAME", "datatype": ["ascii", 20]})

    vector_columns, meta_columns = determine_vector_and_meta_columns(
        in_cols, out_cols
    )

    # Check that the names ended up in the right place
    input_names = [s["name"] for s in in_cols]
    output_names = [s["name"] for s in out_cols]
    vector_names = [s[0].upper() for s in vector_columns]
    meta_names = [s[0].upper() for s in meta_columns]

    assert np.all(np.isin(vector_names, input_names))
    assert set(vector_names + meta_names) == set(output_names)
    assert set(vector_names).isdisjoint(set(meta_names))
    assert meta_names == ["SOURCE_ID", "NAME"]
    
    # test that the datatypes are all numpy types
    for _name, dtype in vector_columns:
        assert np.issubdtype(dtype, np.number)
    for _name, dtype in meta_columns:
        assert (np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.bytes_))


@pytest.mark.parametrize("defaults", [None, [np.nan, 1.0, "foo"]])
def test_make_empty_recarray(defaults):
    n_rows = 10
    n_sources = 5
    vector_columns = [('FLUX', np.float32), ('WAVELENGTH', np.float64)]
    meta_columns = [('NAME', np.str_)]
    columns = vector_columns + meta_columns
    is_vector = [True] * len(vector_columns) + [False] * len(meta_columns)
    
    recarray = make_empty_recarray(n_rows, n_sources, columns, is_vector)
    
    # Check the shapes
    assert recarray.shape == (n_sources,)
    assert recarray["FLUX"].shape == (n_sources, n_rows)
    assert recarray["WAVELENGTH"].shape == (n_sources, n_rows)
    assert recarray["NAME"].shape == (n_sources,)
    
    # Check the data types of the columns
    for name, dtype in vector_columns + meta_columns:
        assert recarray[name].dtype == dtype
    
    # If no defaults, array should be empty (size zero)
    if defaults is None:
        assert recarray[name].size == 0
    # Otherwise, check the values are equal to defaults
    else:
        for i, (name, dtype) in enumerate(columns):
            default = defaults[i]
            if isinstance(default, str):
                # allclose has trouble with string comparisons
                for elem in dm.table[name]:
                    assert elem.decode("utf-8") == default
            else:
                assert np.allclose(recarray[name], default, equal_nan=True)


@pytest.fixture(params=[False, True])
def empty_recarray(request):
    """
    Make an empty output table.

    The parameter `request.param` is used to determine whether to
    add a metadata column to the output table that is NOT in SpecModel
    in order to test the warning logging for problem columns.
    """
    n_rows = 10
    n_sources = 5
    vector_columns = [('FLUX', np.float32), ('WAVELENGTH', np.float64)]
    meta_columns = [('NAME', "S20"), ('SOURCE_ID', np.int32)]
    if request.param:
        # Add a column that is not in the input model
        meta_columns.append(('EXTRA_COLUMN', np.float32))
    columns = vector_columns + meta_columns
    is_vector = [True] * len(vector_columns) + [False] * len(meta_columns)
    defaults = [0,]*len(columns)
    return make_empty_recarray(n_rows, n_sources, columns, is_vector, defaults)


@pytest.mark.parametrize("ignore_columns", [["SOURCE_ID"], ["FLUX"], None])
def test_populate_recarray(empty_recarray, ignore_columns, monkeypatch):

    # make a log watcher
    watcher = LogWatcher("Metadata could not be determined from input spec_table")
    monkeypatch.setattr(logging.getLogger("jwst.lib.pipe_utils"), "warning", watcher)

    # First re-construct vector and meta columns from the input recarray
    output_table = copy.deepcopy(empty_recarray)
    all_columns = output_table.dtype.names
    all_datatypes = [output_table[name].dtype for name in all_columns]
    vector_columns = [(name, dtype) for name, dtype in zip(all_columns[:2], all_datatypes[:2])]
    meta_columns = [(name, dtype) for name, dtype in zip(all_columns[2:], all_datatypes[2:])]
    (n_sources, n_rows) = output_table["FLUX"].shape

    for i in range(n_sources):
        input_spec = dm.SpecModel()
        input_spec.name = f"Source {i}"
        input_spec.source_id = i

        this_output = output_table[i]

        # hack input model schema to have only the flux and wavelength columns
        input_spec.schema["properties"]["spec_table"]["datatype"] = \
            input_spec.schema["properties"]["spec_table"]["datatype"][:2]

        # give all the input spectra different numbers of rows
        input_spec.spec_table = np.ones((n_rows-i,), dtype=[(name, dtype) for name, dtype in vector_columns])
        populate_recarray(
            this_output,
            input_spec,
            n_rows,
            vector_columns,
            meta_columns,
            ignore_columns=ignore_columns
        )

    if ignore_columns is None:
        ignore_columns = []
        
    # Check that the output table has been populated correctly
    for i in range(n_sources):
        for name, _ in vector_columns:
            if name not in ignore_columns:
                # ensure proper nan-padding of output
                expected = np.ones((n_rows,))*np.nan
                expected[:n_rows-i] = 1.0
                assert np.array_equal(output_table[name][i], expected, equal_nan=True)
            else:
                # These should all be set to default
                assert np.allclose(output_table[name][i], 0.0, equal_nan=False)
        assert output_table["NAME"][i] == np.bytes_(f"Source {i}")
        if "SOURCE_ID" not in ignore_columns:
            assert output_table["SOURCE_ID"][i] == i
        else:
            # If ignored, it should be set to default
            assert output_table["SOURCE_ID"][i] == 0

    # Check for "problems" warning message
    if "EXTRA_COLUMN" in all_columns and "EXTRA_COLUMN" not in ignore_columns:
        watcher.assert_seen()
