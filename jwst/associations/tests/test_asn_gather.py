"""Test asn_gather functionality"""
import pytest

from pathlib import Path

from jwst import associations as asnpkg
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.load_as_asn import LoadAsAssociation
from jwst.lib.file_utils import pushdir

# Testing constants
PRIMARY_NAME = 'primary'
PRIMARY_PATH = PRIMARY_NAME + '_asn.json'

@pytest.fixture(scope='module')
def source_folder(tmp_path_factory):
    """Create a set of source associations"""
    primary_members = [
        'primary_1.txt',
        'primary_2.txt',
        'primary_3.txt',
        'primary_4.txt',
    ]
    primary = asn_from_list(primary_members, product_name=PRIMARY_NAME)
    source_folder = tmp_path_factory.mktemp('asn_gather_source')
    with pushdir(source_folder):
        _, serialized = primary.dump()
        with open(PRIMARY_PATH, 'w') as fh:
            fh.write(serialized)

    return source_folder


@pytest.fixture(scope='module')
def gather(source_folder, tmp_path_factory):
    """Do the actual gathering"""
    dest_folder = tmp_path_factory.mktemp('asn_gather_dest')
    with pushdir(dest_folder):
        asn_path = asnpkg.asn_gather(source_folder / PRIMARY_PATH)

    return dest_folder, asn_path, source_folder


def test_ispath(gather):
    """Test that an association is generated"""
    dest_folder, asn_path, source_folder = gather

    assert isinstance(asn_path, Path), f'Return a Path: {asn_path}'

def test_all_members(gather):
    """Test to ensure all members are accounted for"""
    dest_folder, asn_path, source_folder = gather

    source_asn = LoadAsAssociation.load(source_folder / PRIMARY_PATH)
    asn = LoadAsAssociation.load(dest_folder / PRIMARY_PATH)

    assert len(source_asn['products']) == len(asn['products'])
    for source_product, product in zip(source_asn['products'], asn['products']):
        assert len(source_product['members']) == len(product['members'])


def test_copy(gather):
    """Test that members copied"""
