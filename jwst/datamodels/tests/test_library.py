import json

import pytest
import stdatamodels.jwst.datamodels
from stdatamodels.jwst.datamodels import ImageModel

from stpipe.library import NoGroupID
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.load_as_asn import load_asn
from jwst.datamodels.library import ModelLibrary
import jwst.datamodels as dm


# for the example association, set 2 different observation numbers
# so the association will have 2 groups (since all other group_id
# determining meta is the same, see `example_asn_path`)
_OBSERVATION_NUMBERS = ["1", "1", "2"]
_N_MODELS = len(_OBSERVATION_NUMBERS)
_PRODUCT_NAME = "foo_out"
_POOL_NAME = "some_pool"


@pytest.fixture
def example_asn_path(tmp_path):
    """
    Fixture that creates a simple association, saves it (and the models)
    to disk, and returns the path of the saved association
    """
    fns = []
    for i in range(_N_MODELS):
        m = ImageModel()
        m.meta.observation.program_number = "0001"
        m.meta.observation.observation_number = _OBSERVATION_NUMBERS[i]
        m.meta.observation.visit_number = "1"
        m.meta.observation.visit_group = "1"
        m.meta.observation.sequence_id = "01"
        m.meta.observation.activity_id = "1"
        m.meta.observation.exposure_number = "1"
        m.meta.instrument.name = "NIRCAM"
        m.meta.instrument.channel = "SHORT"
        base_fn = f"{i}.fits"
        m.meta.filename = base_fn
        m.save(str(tmp_path / base_fn))
        fns.append(base_fn)

    asn = asn_from_list(fns, product_name=_PRODUCT_NAME)
    base_fn, contents = asn.dump(format="json")
    contents_as_dict = json.loads(contents)
    contents_as_dict["asn_pool"] = _POOL_NAME
    contents = json.dumps(contents_as_dict)
    asn_filename = tmp_path / base_fn
    with open(asn_filename, "w") as f:
        f.write(contents)
    return asn_filename


@pytest.fixture
def example_library(request, example_asn_path):
    """
    Fixture that builds off of `example_asn_path` and returns a
    library created from the association with default options
    """
    return ModelLibrary(example_asn_path, on_disk=request.param)


def _set_custom_member_attr(example_asn_path, member_index, attr, value):
    """
    Helper function to modify the association at `example_asn_path`
    by adding an attribute `attr` to the member list (at index
    `member_index`) with value `value`. This is used to modify
    the `group_id` or `exptype` of a certain member for some tests.
    """
    with open(example_asn_path, "r") as f:
        asn_data = load_asn(f)
    asn_data["products"][0]["members"][member_index][attr] = value
    with open(example_asn_path, "w") as f:
        json.dump(asn_data, f)


@pytest.mark.parametrize(
    "example_library",
    [
        False,
    ],
    indirect=True,
)
def test_load_asn(request, example_library, example_asn_path):
    """
    Test that __len__ returns the number of models/members loaded
    from the association (and does not require opening the library).

    Test that the asn_dir and on_disk properties are set correctly.
    """
    assert len(example_library) == _N_MODELS

    expected_asn_dir = str(example_asn_path.parent)
    assert example_library.asn_dir == expected_asn_dir
    assert example_library.on_disk == False

    # test that asn_dir is not settable via public API
    with pytest.raises(AttributeError):
        example_library.asn_dir = "foo"

    # test that on_disk is not settable via public API
    with pytest.raises(AttributeError):
        example_library.on_disk = True


@pytest.mark.parametrize("attr", ["group_names", "group_indices"])
def test_group_with_no_datamodels_open(example_asn_path, attr, monkeypatch):
    """
    Test that the "grouping" methods do not call datamodels.open
    """

    # patch datamodels.open to always raise an exception
    # this will serve as a smoke test to see if any of the attribute
    # accesses (or instance creation) attempts to open models
    def no_open(*args, **kwargs):
        raise Exception()

    monkeypatch.setattr(stdatamodels.jwst.datamodels, "open", no_open)

    # use example_asn_path here to make the instance after we've patched
    # datamodels.open
    library = ModelLibrary(example_asn_path)
    getattr(library, attr)


@pytest.mark.parametrize(
    "asn_group_id, meta_group_id, expected_group_id",
    [
        ("42", None, "42"),
        (None, "42", "42"),
        ("42", "26", "42"),
    ],
)
def test_group_id_override(example_asn_path, asn_group_id, meta_group_id, expected_group_id):
    """
    Test that overriding a models group_id via:
        - the association member entry
        - the model.meta.group_id
    overwrites the automatically calculated group_id (with the asn taking precedence)
    """
    if asn_group_id:
        _set_custom_member_attr(example_asn_path, 0, "group_id", asn_group_id)
    if meta_group_id:
        model_filename = example_asn_path.parent / "0.fits"
        with dm.open(model_filename) as model:
            model.meta.group_id = meta_group_id
            model.save(model_filename)
    library = ModelLibrary(example_asn_path)
    group_names = library.group_names
    assert len(group_names) == 3
    assert expected_group_id in group_names
    with library:
        model = library.borrow(0)
        assert model.meta.group_id == expected_group_id
        library.shelve(model, 0, modify=False)


@pytest.mark.parametrize("example_library", [True, False], indirect=True)
def test_asn_attributes_assignment(example_library):
    expected_table_name = "jwnoprogram-a3001"
    assert example_library.asn["table_name"].startswith(expected_table_name)
    assert example_library.asn["asn_pool"] == _POOL_NAME

    # test that the association attributes are assigned to the models
    with example_library:
        for i in range(_N_MODELS):
            model = example_library.borrow(i)
            assert model.meta.asn.table_name.startswith(expected_table_name)
            assert model.meta.asn.pool_name == _POOL_NAME
            example_library.shelve(model, i, modify=False)


@pytest.mark.parametrize("modify", [True, False])
@pytest.mark.parametrize("example_library", [True, False], indirect=True)
def test_get_crds_parameters(example_library, modify):
    """
    Test the JWST override to get_crds_parameters.

    Ensure that parameters are up-to-date with the models in the library
    if those models have been borrowed before.
    """
    if modify:
        with example_library:
            for i in range(_N_MODELS):
                model = example_library.borrow(i)
                model.meta.instrument.name = "MIRI"
                example_library.shelve(model, i, modify=True)
    params = example_library.get_crds_parameters()
    assert isinstance(params, dict)

    instrument_name = params["meta.instrument.name"]
    channel = params["meta.instrument.channel"]
    assert channel == "SHORT"
    assert instrument_name == "MIRI" if modify else "NIRCAM"


@pytest.mark.parametrize(
    "example_library",
    [
        False,
    ],
    indirect=True,
)
def test_model_to_group_id(example_library):
    """
    Test that the model_to_group_id method returns the correct group_id
    based on the model's meta.observation attributes.
    """
    with example_library:
        model = example_library.borrow(0)
        # test that if a group_id is already set, it is returned unmodified
        assert model.meta.hasattr("group_id")
        group_id_0 = example_library._model_to_group_id(model)

        # test that if the group_id is not set, it is calculated
        # from the meta.observation attributes in the same way as when library was initialized
        del model.meta.group_id
        assert not model.meta.hasattr("group_id")
        group_id = example_library._model_to_group_id(model)
        assert group_id == group_id_0

        # test error raise from missing key
        del model.meta.observation.program_number
        with pytest.raises(NoGroupID):
            example_library._model_to_group_id(model)

        # test error raise from missing meta.observation
        del model.meta.observation
        with pytest.raises(NoGroupID):
            example_library._model_to_group_id(model)

        example_library.shelve(model, 0, modify=False)
