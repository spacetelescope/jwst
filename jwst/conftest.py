"""Project default for pytest"""
import os
import tempfile
import pytest
import importlib

from astropy.io import fits
from gwcs import wcs

from .associations import AssociationRegistry, AssociationPool
from .associations.tests.helpers import t_path
from .lib.tests import helpers as lib_helpers
from .lib import s3_utils
from . import datamodels
from .assign_wcs.assign_wcs_step import AssignWcsStep


@pytest.fixture(scope='session')
def full_pool_rules():
    """Setup to use the full example pool and registry"""
    pool_fname = t_path('data/mega_pool.csv')
    pool = AssociationPool.read(pool_fname)
    rules = AssociationRegistry()

    return pool, rules, pool_fname


@pytest.fixture
def mk_tmp_dirs():
    """Create a set of temporary directorys and change to one of them."""
    tmp_current_path = tempfile.mkdtemp()
    tmp_data_path = tempfile.mkdtemp()
    tmp_config_path = tempfile.mkdtemp()

    old_path = os.getcwd()

    try:
        os.chdir(tmp_current_path)
        yield tmp_current_path, tmp_data_path, tmp_config_path

    finally:
        os.chdir(old_path)


@pytest.fixture(autouse=True)
def monkey_patch_s3_client(monkeypatch):
    monkeypatch.setattr(s3_utils, "_CLIENT", lib_helpers.MockS3Client())


@pytest.fixture
def slow(request):
    """Setup slow fixture for tests to identify if --slow
    has been specified
    """
    return request.config.getoption('--slow')


# assign_wcs
class DataBundle:
    """Class that creates and contains several objects required for assign_wcs tests. Builds the requested datamodel \
    instance from a fits HDUList that contains given keywords for the primary HDU and SCI HDU. Collects references and
    from the datamodel instance and references, creates a WCS object.
    """

    # noinspection PyTypeChecker
    def __init__(self, primary_hdu_keys: dict, sci_hdu_keys: dict, model_class_name: str,
                 hdu_class_name: str = 'ImageHDU'):
        self.hdu_list = fits.HDUList()
        self.datamodel = None
        self.references = None
        self.wcs_object = None

        if 'instrume' not in primary_hdu_keys:
            raise KeyError('INSTRUME keyword is required in the primary_hdu_keys to build the correct WCS object.')

        self.create_hdu_list(primary_hdu_keys, sci_hdu_keys, hdu_class_name)
        self.create_datamodel(self.hdu_list, model_class_name)
        self.get_references(self.datamodel)
        self.create_wcs(self.hdu_list[0].header['instrume'], self.datamodel, self.references)

    def create_hdu_list(self, primary_hdu_keys: dict, sci_hdu_keys: dict, hdu_class_name: str) -> fits.HDUList:
        """Create a basic HDUList for use in generating a DataModel instance for testing."""
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header.update(primary_hdu_keys)

        sci_hdu_class = getattr(fits, hdu_class_name)
        sci_hdu = sci_hdu_class()
        sci_hdu.header['EXTNAME'] = 'SCI'
        sci_hdu.header.update(sci_hdu_keys)

        self.hdu_list.extend([primary_hdu, sci_hdu])

        return self.hdu_list

    def create_datamodel(self, hdu_list: fits.HDUList, model_class_name: str) -> datamodels.DataModel:
        """Create a DataModel instance from an HDUList for testing."""
        model_class = getattr(datamodels, model_class_name)
        self.datamodel = model_class(hdu_list)

        return self.datamodel

    def get_references(self, datamodel: datamodels.DataModel) -> dict:
        """Get references required for AssignWcsStep for the given DataModel."""
        step = AssignWcsStep()
        self.references = {
            reftype: step.get_reference_file(datamodel, reftype) for reftype in AssignWcsStep.reference_file_types
        }

        return self.references

    def create_wcs(self, instrument: str, datamodel: datamodels.DataModel, references: dict) -> wcs.WCS:
        """Create a WCS instance for use in testing."""
        module = importlib.import_module(f'.assign_wcs.{instrument.lower()}', package='jwst')
        create_pipeline = getattr(module, 'create_pipeline')

        pipeline = create_pipeline(datamodel, references)
        self.wcs_object = wcs.WCS(pipeline)

        return self.wcs_object


@pytest.fixture
def assign_wcs_objects():
    def _create_data_bundle(*args, **kwargs):
        return DataBundle(*args, **kwargs)

    return _create_data_bundle
