import os
import warnings

import pytest

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import JwstDataModel, MultiExposureModel
from stdatamodels.jwst.datamodels import _defined_models as defined_models

from jwst.datamodels import ModelContainer, SourceModelContainer
from jwst.associations.asn_from_list import asn_from_list

from jwst.lib.file_utils import pushdir


ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')
FITS_FILE = os.path.join(ROOT_DIR, 'test.fits')
ASN_FILE = os.path.join(ROOT_DIR, 'association.json')


@pytest.fixture
def container():
    warnings.simplefilter("ignore")
    asn_file_path, asn_file_name = os.path.split(ASN_FILE)
    with pushdir(asn_file_path):
        with ModelContainer(asn_file_name) as c:
            for m in c:
                m.meta.observation.program_number = '0001'
                m.meta.observation.observation_number = '1'
                m.meta.observation.visit_number = '1'
                m.meta.observation.visit_group = '1'
                m.meta.observation.sequence_id = '01'
                m.meta.observation.activity_id = '1'
                m.meta.observation.exposure_number = '1'
                m.meta.instrument.name = 'NIRCAM'
                m.meta.instrument.channel = 'SHORT'
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
        assert model.meta.telescope == 'JWST'


def test_modelcontainer_indexing(container):
    assert isinstance(container[0], JwstDataModel)


def test_modelcontainer_group1(container):
    for group in container.models_grouped:
        assert len(group) == 2
        for model in group:
            pass


def test_modelcontainer_group2(container):
    container[0].meta.observation.exposure_number = '2'
    for group in container.models_grouped:
        assert len(group) == 1
        for model in group:
            pass
    container[0].meta.observation.exposure_number = '1'


def test_modelcontainer_group_names(container):
    assert len(container.group_names) == 1
    reset_group_id(container)
    container[0].meta.observation.exposure_number = '2'
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


@pytest.mark.parametrize("model", [v for v in defined_models.values()])
def test_all_datamodels_init(model):
    """
    Test that all current datamodels can be initialized.
    """
    if model is SourceModelContainer:
        # SourceModelContainer cannot have init=None
        model(MultiExposureModel())
    else:
        model()


def test_model_container_ind_asn_exptype(container):
    ind = container.ind_asn_type('science')
    assert ind == [0, 1]
