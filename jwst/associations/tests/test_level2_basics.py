"""Test basic usage of Level2 associations"""
from __future__ import absolute_import

import pytest
import re

from .helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)

from .. import (
    generate,
    load_asn,
)
from ..main import Main

NONSSCIENCE = ['background']
REGEX_LEVEL2A = '(?P<path>.+)(?P<type>_rate(ints)?)(?P<extension>\..+)'


def from_level2_schema():
    with open(t_path('data/asn_level2.json')) as asn_file:
        asn = load_asn(asn_file)
    return [asn]


def generate_from_pool(pool_path):
    """Generate associations from pools"""
    rules = registry_level2_only()
    pool = combine_pools(t_path(pool_path))
    asns, orphaned = generate(pool, rules)
    return asns


def cmd_from_pool(pool_path, args):
    """Run commandline on pool

    Parameters
    ---------
    pool_path: str
        The pool to run on.

    args: [arg(, ...)]
        Additional command line arguments in the form `sys.argv`
    """
    full_args = [
        '--dry-run',
        '-D',
        '-r',
        t_path('../lib/rules_level2b.py'),
        '--ignore-default'
    ]
    full_args.extend(args)
    result = Main(full_args, pool=pool_path)
    return result


@pytest.mark.parametrize(
    'pool_path, args, n_asns, n_products, n_members',
    [
        (
            t_path('data/pool_010_spec_nirspec_lv2bkg.csv'),
            [],
            1,
            [8],
            [3]*8
        )
    ]
)
def test_from_cmdline(pool_path, args, n_asns, n_products, n_members):
    """Test results from commandline

    Parameters
    ----------
    pool_path: str
        Pool to run on.

    args: [arg(, ...)]
        Other command-line arguments

    n_asns: int
        Number of associations that should be produced

    n_products: [int(, ...)]
        Number of products for each association

    n_members: [int(, ...)]
        Number of members for each product
    """
    pool = combine_pools([pool_path])
    results = cmd_from_pool(pool, args)
    assert len(results.associations) == n_asns
    for asn in results.associations:
        products = asn['products']
        assert len(products) == n_products
        for product in products:
            assert len(product['members']) == n_members


@pytest.mark.parametrize(
    'asns, n_asns, n_products, asn_type, required_rules, required_members',
    [
        (
            from_level2_schema(), 1, 6,
            'image2', ['Asn_Lv2Image'], []
        ),
        (
            generate_from_pool('data/pool_002_image_miri.csv'), 8, 1,
            'image2', ['Asn_Lv2Image'], []
        ),
        (
            generate_from_pool('data/pool_007_spec_miri.csv'), 11, 1,
            'spec2', ['Asn_Lv2Spec'], []
        ),
        (
            generate_from_pool('data/pool_011_spec_miri_lv2bkg_lrs.csv'), 1, 1,
            'spec2', ['Asn_Lv2Spec', 'Asn_Lv2SpecBkg'], ['background']
        ),
        (
            generate_from_pool('data/pool_009_spec_miri_lv2bkg.csv'), 9, 1,
            'spec2', ['Asn_Lv2Spec', 'Asn_Lv2SpecBkg'], ['background']
        ),
        (
            generate_from_pool('data/pool_010_spec_nirspec_lv2bkg.csv'), 8, 1,
            'spec2', ['Asn_Lv2Spec', 'Asn_Lv2SpecBkg'], ['background']
        ),
        (
            generate_from_pool('data/pool_015_spec_nirspec_lv2bkg_reversed.csv'), 8, 1,
            'spec2', ['Asn_Lv2Spec', 'Asn_Lv2SpecBkg'], ['background']
        ),
        (
            generate_from_pool('data/pool_016_spec_nirspec_lv2bkg_double.csv'), 8, 1,
            'spec2', ['Asn_Lv2Spec', 'Asn_Lv2SpecBkg'], ['background']
        ),
    ]
)
def test_level2(
        asns,
        n_asns,
        n_products,
        asn_type,
        required_rules,
        required_members
):
    assert len(asns) == n_asns
    for asn in asns:
        assert asn['asn_type'] == asn_type
        if len(required_rules) > 0:
            assert asn['asn_rule'] in required_rules
        products = asn['products']
        assert len(products) == n_products
        for product in products:
            members = product['members']
            assert len(members) > 0
            science = [
                member
                for member in members
                if member['exptype'] == 'science'
            ]
            assert len(science) == 1
            if len(members) > 1:
                nonscience = [
                    member
                    for member in members
                    if member['exptype'] != 'science'
                ]
                assert len(nonscience) > 0
                exptypes = set(
                    member['exptype']
                    for member in nonscience
                )
                assert exptypes.issubset(NONSSCIENCE)
                if len(required_members) > 0:
                    assert set(required_members).issubset(exptypes)


def test_level2_productname():
    asns = generate_from_pool('data/pool_002_image_miri.csv')
    for asn in asns:
        for product in asn['products']:
            science = [
                member
                for member in product['members']
                if member['exptype'] == 'science'
            ]
            assert len(science) == 1
            match = re.match(REGEX_LEVEL2A, science[0]['expname'])
            assert match.groupdict()['path'] == product['name']
