import os
import warnings

import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_filename
from stdatamodels.exceptions import NoTypeWarning
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import JwstDataModel

from jwst.associations import load_as_asn
from jwst.associations.asn_from_list import asn_from_list
from jwst.datamodels import ModelContainer
from jwst.lib.file_utils import pushdir

FITS_FILE = get_pkg_data_filename("data/test.fits", package="jwst.datamodels.tests")
ASN_FILE = get_pkg_data_filename("data/association.json", package="jwst.datamodels.tests")
CUSTOM_GROUP_ID_ASN_FILE = get_pkg_data_filename(
    "data/association_group_id.json", package="jwst.datamodels.tests"
)


@pytest.fixture
def container():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", NoTypeWarning)
        asn_file_path, asn_file_name = os.path.split(ASN_FILE)
        with pushdir(asn_file_path):
            with ModelContainer(asn_file_name) as c:
                for m in c:
                    m.meta.observation.program_number = "0001"
                    m.meta.observation.observation_number = "1"
                    m.meta.observation.visit_number = "1"
                    m.meta.observation.visit_group = "1"
                    m.meta.observation.sequence_id = "01"
                    m.meta.observation.activity_id = "1"
                    m.meta.observation.exposure_number = "1"
                    m.meta.instrument.name = "NIRCAM"
                    m.meta.instrument.channel = "SHORT"
            yield c


def reset_group_id(container):
    """Remove group_id from all models in container"""
    for m in container:
        try:
            del m.meta.group_id
        except AttributeError:
            pass


def test_modelcontainer_iteration(container):
    for model in container:
        assert model.meta.telescope == "JWST"


def test_modelcontainer_indexing(container):
    assert isinstance(container[0], JwstDataModel)


def test_modelcontainer_group1(container):
    for group in container.models_grouped:
        assert len(group) == 2
        for model in group:
            pass


def test_modelcontainer_group2(container):
    container[0].meta.observation.exposure_number = "2"
    for group in container.models_grouped:
        assert len(group) == 1
        for model in group:
            pass
    container[0].meta.observation.exposure_number = "1"


def test_modelcontainer_group_names(container):
    assert len(container.group_names) == 1
    reset_group_id(container)
    container[0].meta.observation.exposure_number = "2"
    assert len(container.group_names) == 2


def test_modelcontainer_error_from_asn(tmp_path):
    asn = asn_from_list(["foo.fits"], product_name="foo_out")
    name, serialized = asn.dump(format="json")
    # The following Path object needs to be stringified because
    # datamodels.open() doesn't deal with pathlib objects when they are
    # .json files
    path = str(tmp_path / name)
    with open(path, "w") as f:
        f.write(serialized)

    # The foo.fits file doesn't exist
    with pytest.raises(FileNotFoundError):
        datamodels.open(path)


def test_model_container_ind_asn_exptype(container):
    ind = container.ind_asn_type("science")
    assert ind == [0, 1]


def test_group_id(tmp_path):
    c = ModelContainer(CUSTOM_GROUP_ID_ASN_FILE)
    groups = list(c.models_grouped)

    assert len(groups) == 5
    assert sorted(map(len, groups)) == [1, 1, 1, 1, 4]

    with open(CUSTOM_GROUP_ID_ASN_FILE) as f:
        asn_data = load_as_asn.load_asn(f)

    asn_group_ids = set()
    for d in asn_data["products"][0]["members"]:
        group_id = d.get("group_id")
        if group_id is None:
            asn_group_ids.add("jw00001001001_02201_00001")
        else:
            asn_group_ids.add(group_id)

    model_droup_ids = set()
    for g in groups:
        for m in g:
            model_droup_ids.add(m.meta.group_id)

    assert asn_group_ids == model_droup_ids


@pytest.mark.parametrize("path", ["foo", "foo.fits"])
def test_save(tmp_cwd, container, path):
    # container pushes us to data/ directory so need to go back to tmp_cwd
    # to avoid polluting the data/ directory
    with pushdir(tmp_cwd):
        # test default just saves things at model meta filename
        container.save()
        expected_fnames = []
        for model in container:
            expected_fnames.append(model.meta.filename)
        for fname in expected_fnames:
            assert os.path.exists(fname)

        # test specifying path saves to custom path with indices
        container.save(path)
        expected_fnames = [
            path.replace(".fits", "") + str(i) + ".fits" for i in range(len(container))
        ]
        for fname in expected_fnames:
            assert os.path.exists(fname)

        # test saving path when the container has length 1: no index appended
        container = ModelContainer([container[0]])
        container.save(path)
        expected_fname = path.replace(".fits", "") + ".fits"
        assert os.path.exists(expected_fname)


def test_open_guess(container):
    """Test that `guess` keyword argument works in ModelContainer."""
    asn_file_path, _asn_file_name = os.path.split(ASN_FILE)
    fnames = [m.meta.filename for m in container]
    with pushdir(asn_file_path):
        # opening it normally works fine
        ModelContainer(fnames, guess=True)
        with pytest.raises(
            TypeError, match="Model type is not specifically defined and guessing has been disabled"
        ):
            # but if you don't allow guessing the model type, it raises TypeError
            ModelContainer(fnames, guess=False)


def test_open_kwargs(container):
    wrong_schema = datamodels.NRMModel()._schema
    asn_file_path, _asn_file_name = os.path.split(ASN_FILE)
    fnames = [m.meta.filename for m in container]
    with pushdir(asn_file_path):
        # opening it normally works fine
        ModelContainer(fnames)
        with pytest.raises(AttributeError):
            # but schema can be passed all the way through to DataModel.__init__ on the
            # individual datamodels, and cause AttributeError
            ModelContainer(fnames, schema=wrong_schema)


def test_copy(container):
    # make a deep copy of a container
    container_copy = container.copy()

    assert container.asn_table is not container_copy.asn_table
    assert container.asn_table == container_copy.asn_table
    assert container.asn_exptypes == container_copy.asn_exptypes
    assert container.asn_n_members == container_copy.asn_n_members
    assert container.asn_table_name == container_copy.asn_table_name
    assert container.asn_pool_name == container_copy.asn_pool_name
    assert container.asn_file_path == container_copy.asn_file_path

    assert len(container._models) == len(container_copy._models)
    for in_exp, out_exp in zip(container._models, container_copy._models, strict=True):
        assert in_exp is not out_exp
        assert in_exp.meta.filename == out_exp.meta.filename
        assert np.all(in_exp.data == out_exp.data)

        # Equal comparison fails when underlying instance is not a shallow copy
        with pytest.raises(ValueError, match="truth value of an array"):
            assert in_exp == out_exp
