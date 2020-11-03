"""Test asn_gather functionality"""
import pytest

from jwst import associations as asnpkg
from jwst.associations.asn_from_list import asn_from_list
from jwst.lib.file_utils import pushdir

# Testing constants
PRIMARY_NAME = 'primary'


@pytest.fixture(scope='module')
def source_path(tmp_path_factory):
    """Create a set of source associations"""
    primary_members = [
        'primary_1.txt',
        'primary_2.txt',
        'primary_3.txt',
        'primary_4.txt',
    ]
    primary = asn_from_list(primary_members, product_name=PRIMARY_NAME)
    source_path = tmp_path_factory.mktemp('asn_gather_source')
    with pushdir(source_path):
        primary_path = PRIMARY_NAME + '_asn.json'
        _, serialized = primary.dump()
        with open(primary_path, 'w') as fh:
            fh.write(serialized)

    return source_path


@pytest.fixture(scope='module')
def gather(source_path, tmp_path_factory):
    """Do the actual gathering"""
    dest_path = tmp_path_factory.mktemp('asn_gather_dest')
    with pushdir(dest_path):
        asn = asnpkg.asn_gather(source_path)

    return dest_path, asn, source_path


def test_isasn(gather):
    """Test that an association is generated"""
    asn_path, asn, source_asn_path = gather

    assert isinstance(asn, (asnpkg.Association, dict)), f'Returned association is not an Association: {asn}'
