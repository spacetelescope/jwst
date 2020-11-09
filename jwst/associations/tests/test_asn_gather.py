"""Test asn_gather functionality"""
import pytest

from pathlib import Path

from jwst.associations import asn_gather
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.load_as_asn import LoadAsAssociation
from jwst.lib.file_utils import pushdir

# Testing constants
PRIMARY_STEM = 'primary'
PRIMARY_NAME = PRIMARY_STEM + '_asn.json'

@pytest.fixture(scope='module')
def source_folder(tmp_path_factory):
    """Create a set of source associations"""
    primary_members = [
        ('sci_1.txt', 'science'),
        ('sci_2.txt', 'science'),
        ('bkg_1.txt', 'background'),
        ('imprint_1.txt', 'imprint'),
    ]

    source_folder = tmp_path_factory.mktemp('asn_gather_source')
    with pushdir(source_folder):

        # Create all the files
        for expname, exptype in primary_members:
            with open(expname, 'w') as fh:
                fh.write(expname)

        # Create the association
        primary = asn_from_list(primary_members, product_name=PRIMARY_STEM, with_exptype=True)
        _, serialized = primary.dump()
        with open(PRIMARY_NAME, 'w') as fh:
            fh.write(serialized)

    return source_folder


@pytest.fixture(scope='module')
def gather(source_folder, tmp_path_factory):
    """Do the actual gathering"""
    dest_folder = tmp_path_factory.mktemp('asn_gather_dest')
    asn_path = asn_gather.asn_gather(source_folder / PRIMARY_NAME, destination=dest_folder)

    return dest_folder, asn_path, source_folder


def test_ispath(gather):
    """Test that an association is generated"""
    dest_folder, asn_path, source_folder = gather

    assert isinstance(asn_path, Path), f'Return a Path: {asn_path}'

def test_all_members(gather):
    """Test to ensure all members are accounted for"""
    dest_folder, asn_path, source_folder = gather

    source_asn = LoadAsAssociation.load(source_folder / PRIMARY_NAME)
    asn = LoadAsAssociation.load(asn_path)

    assert len(source_asn['products']) == len(asn['products'])
    for source_product, product in zip(source_asn['products'], asn['products']):
        assert len(source_product['members']) == len(product['members'])


def test_copy(gather):
    """Test that members are copied"""
    dest_folder, asn_path, source_folder = gather

    asn = LoadAsAssociation.load(asn_path)

    for product in asn['products']:
        for member in product['members']:
            path = Path(dest_folder / member['expname'])
            assert path.exists(), f'File does not exist {path}'
