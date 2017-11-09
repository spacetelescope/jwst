"""Test against Level2 standard associations

Notes
-----
Most of the standard associations which are compared
against are built in the jupyter notebook

./notebooks/make_tests.ipynb
"""
from __future__ import absolute_import
from glob import glob
from os import path
import pytest

from .helpers import (
    combine_pools,
    compare_asns,
    t_path,
)

from .. import load_asn
from ..main import Main

# Main test args
TEST_ARGS = ['--dry-run']

# Produce Level2b only associations
LV2_ONLY_ARGS = [
    '-r',
    t_path('../lib/rules_level2b.py'),
    '--ignore-default'
]

# Produce Level3 only associations
LV3_ONLY_ARGS = [
    '-r',
    t_path('../lib/rules_level3.py'),
    '--ignore-default'
]

# Produce general associations
GENERAL_ARGS = []


@pytest.yield_fixture(
    scope='module',
    params=[
        (
            LV3_ONLY_ARGS,
            'pool_002_image_miri'
        ),
        (
            GENERAL_ARGS,
            'pool_004_wfs'
        ),
        (
            GENERAL_ARGS,
            'pool_005_spec_niriss'
        ),
        (
            GENERAL_ARGS,
            'pool_006_spec_nirspec'
        ),
        (
            GENERAL_ARGS,
            'pool_007_spec_miri'
        ),
        (
            LV2_ONLY_ARGS,
            'pool_009_spec_miri_lv2bkg'
        ),
        (
            LV2_ONLY_ARGS,
            'pool_010_spec_nirspec_lv2bkg'
        ),
        (
            LV2_ONLY_ARGS,
            'pool_011_spec_miri_lv2bkg_lrs'
        ),
        (
            GENERAL_ARGS,
            'pool_013_coron_nircam'
        ),
        (
            GENERAL_ARGS,
            'pool_014_ami_niriss'
        ),
        (
            LV2_ONLY_ARGS,
            'pool_015_spec_nirspec_lv2bkg_reversed'
        ),
        (
            LV2_ONLY_ARGS,
            'pool_016_spec_nirspec_lv2bkg_double'
        ),
        (
            LV2_ONLY_ARGS,
            'pool_017_spec_nirspec_lv2imprint'
        ),
        (
            LV2_ONLY_ARGS,
            'pool_018_all_exptypes'
        ),
        (
            GENERAL_ARGS,
            'pool_019_niriss_wfss'
        ),
        (
            GENERAL_ARGS,
            'pool_021_tso'
        ),
        (
            GENERAL_ARGS,
            'pool_022_tso_noflag'
        ),
    ]
)
def generate_asns(request):
    """Test exp_type inclusion based on standard associations"""
    main_args, pool_root = request.param

    standards_paths = glob(t_path(path.join('data', pool_root + '*_asn.json')))
    standards = []
    for standard_path in standards_paths:
        with open(standard_path) as fp:
            asn = load_asn(fp)
        standards.append(asn)

    pool_path = t_path(path.join('data', pool_root + '.csv'))
    pool = combine_pools([pool_path])
    args = TEST_ARGS + main_args
    results = Main(args, pool=pool)

    asns = results.associations
    assert len(asns) == len(standards)
    yield asns, standards


def test_against_standard(generate_asns):
    """Compare a generated assocaition against a standard
    """
    generated, standards = generate_asns
    for asn in generated:
        for idx, standard in enumerate(standards):
            try:
                compare_asns(asn, standard)
            except AssertionError as e:
                last_err = e
            else:
                del standards[idx]
                break
        else:
            raise last_err
