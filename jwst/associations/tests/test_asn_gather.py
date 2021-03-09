"""Test asn_gather functionality"""
import pytest

from pathlib import Path
import shutil

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


@pytest.fixture(scope='module',
                params = [
                    (None, None),
                    (['science'], None),
                    (None, 'science')
                ])
def gather_simple(source_folder, tmp_path_factory, request):
    """Do the actual gathering"""
    exp_types, excludes = request.param

    dest_folder = tmp_path_factory.mktemp('asn_gather_dest')
    asn_path = asn_gather.asn_gather(source_folder / PRIMARY_NAME, destination=dest_folder,
                                     exp_types=exp_types, exclude_types=excludes)

    return dest_folder, asn_path, source_folder, exp_types, excludes


@pytest.fixture(scope='module',
                params = [
                    (None, None),
                    (['science'], None),
                    (None, 'science')
                ])
def gather_alternate(source_folder, tmp_path_factory, request):
    """Gather but using an alternate source folder"""
    exp_types, excludes = request.param

    dest_folder = tmp_path_factory.mktemp('asn_gather_alternate_dest')
    asn_folder = tmp_path_factory.mktemp('asn_alternate')
    shutil.copy2(source_folder / PRIMARY_NAME, asn_folder / PRIMARY_NAME)

    asn_path = asn_gather.asn_gather(asn_folder / PRIMARY_NAME, destination=dest_folder,
                                     exp_types=exp_types, exclude_types=excludes,
                                     source_folder=source_folder)

    return dest_folder, asn_path, source_folder, exp_types, excludes


@pytest.fixture
def gather(request, gather_simple, gather_alternate):
    """Parametrize the gather fixtures"""
    type = request.param
    if type == 'gather_simple':
        return gather_simple
    elif type == 'gather_alternate':
        return gather_alternate
    else:
        raise ValueError(f'Unknown gather fixture: "{type}"')


@pytest.mark.parametrize('gather', ['gather_simple', 'gather_alternate'], indirect=True)
def test_ispath(gather, request):
    """Test that an association is generated"""
    dest_folder, asn_path, source_folder, exptypes, excludes = gather

    assert isinstance(asn_path, Path), f'Return a Path: {asn_path}'


@pytest.mark.parametrize('gather', ['gather_simple', 'gather_alternate'], indirect=True)
def test_all_members(gather):
    """Test to ensure all members are accounted for"""
    dest_folder, asn_path, source_folder, exptypes, excludes = gather

    source_asn = LoadAsAssociation.load(source_folder / PRIMARY_NAME)
    asn = LoadAsAssociation.load(asn_path)

    excludes = [] if excludes is None else excludes
    if exptypes is None:
        exptypes = {
            member['exptype']
            for source_product in source_asn['products']
            for member in source_product['members']
        }

    assert len(source_asn['products']) == len(asn['products'])
    for source_product, product in zip(source_asn['products'], asn['products']):
        expected = [
            member
            for member in source_product['members']
            if member['exptype'] in exptypes and member['exptype'] not in excludes
        ]
        assert len(expected) == len(product['members'])


@pytest.mark.parametrize('gather', ['gather_simple', 'gather_alternate'], indirect=True)
def test_copy(gather):
    """Test that members are copied"""
    dest_folder, asn_path, source_folder, exptypes, excludes = gather

    asn = LoadAsAssociation.load(asn_path)

    for product in asn['products']:
        for member in product['members']:
            path = Path(dest_folder / member['expname'])
            assert path.exists(), f'File does not exist {path}'
