"""Test for various WFSS modes"""
from os import path

from jwst.associations import AssociationPool
from jwst.associations.main import Main
from jwst.associations.tests.helpers import t_path

REQUIRED_ASN_TYPES = set(['image2', 'spec2', 'image3', 'spec3'])


def test_niriss_wfss():
    """Test association properties for NIRISS WFSS"""

    pool = AssociationPool.read(
        t_path(path.join('data', 'jw87800_20180412T163456_pool.csv'))
    )
    cmd_args = [
        '--dry-run',
        '--D'
    ]
    results = Main.cli(
        cmd_args,
        pool=pool
    )
    asns = results.associations

    # Need 4 associations: image2, spec2, image3, spec3
    assert len(asns) == 12
    asn_types = [
        asn['asn_type']
        for asn in asns
    ]
    assert REQUIRED_ASN_TYPES == set(asn_types)

    # Arrange associations by type
    asn_by_type = {
        asn['asn_type']: asn
        for asn in asns
    }

    # Ensure catalog and segmentation map names are correct in the spec2 associations
    l3name = asn_by_type['image3']['products'][0]['name']
    source_cat = l3name + '_cat.ecsv'
    segmap = l3name + '_segm.fits'
    for product in asn_by_type['spec2']['products']:
        members_by_type = {
            member['exptype']: member
            for member in product['members']
        }
        assert members_by_type['sourcecat']['expname'] == source_cat
        assert members_by_type['segmap']['expname'] == segmap
