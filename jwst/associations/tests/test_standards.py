"""Test against Level2 standard associations

Notes
-----
Most of the standard associations which are compared
against are built in the jupyter notebook

./notebooks/make_tests.ipynb
"""
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
    '--ignore-default',
    '--no-merge'
]

# Produce Level3 only associations
LV3_ONLY_ARGS = [
    '-r',
    t_path('../lib/rules_level3.py'),
    '--ignore-default'
]

# Do not merge Level2b
LV2_NOMERGE_ARGS = [
    '--no-merge'
]

# Produce general associations
DEF_ARGS = []


# Define the standards
class MakePars():
    def __init__(
            self,
            pool_root,
            main_args=DEF_ARGS,
            source=None,
            outdir=None,
            execute=True,
            xfail=None
    ):
        self.pool_root = pool_root
        self.main_args = main_args
        self.source = source
        self.outdir = outdir
        self.execute = execute
        self.xfail = xfail


standards = [
    MakePars('pool_002_image_miri', main_args=LV3_ONLY_ARGS),
    MakePars('pool_004_wfs'),
    MakePars('pool_005_spec_niriss'),
    MakePars('pool_006_spec_nirspec', main_args=LV2_NOMERGE_ARGS),
    MakePars('pool_007_spec_miri'),
    MakePars('pool_009_spec_miri_lv2bkg'),
    MakePars('pool_010_spec_nirspec_lv2bkg'),
    MakePars('pool_011_spec_miri_lv2bkg_lrs'),
    MakePars('pool_013_coron_nircam', main_args=LV2_NOMERGE_ARGS),
    MakePars('pool_014_ami_niriss'),
    MakePars('pool_015_spec_nirspec_lv2bkg_reversed', main_args=LV2_ONLY_ARGS),
    MakePars('pool_016_spec_nirspec_lv2bkg_double', main_args=LV2_ONLY_ARGS),
    MakePars('pool_017_spec_nirspec_lv2imprint'),
    MakePars('pool_018_all_exptypes', main_args=LV2_ONLY_ARGS),
    MakePars('pool_019_niriss_wfss'),
    MakePars('pool_020_00009_image_miri'),
    MakePars('pool_021_tso'),
    MakePars('pool_022_tso_noflag'),
    MakePars('pool_023_nirspec_msa_3nod', main_args=LV2_ONLY_ARGS),
    MakePars('pool_024_nirspec_fss_nods'),
    MakePars('pool_025_nirspec_fss_nod_chop'),
    MakePars('pool_026_mir_image_tso'),
    MakePars('pool_027_nirspec_ifu_nods', main_args=LV2_NOMERGE_ARGS),
    MakePars('pool_028_mir_lrsfs_nods', main_args=LV2_NOMERGE_ARGS),
]


@pytest.mark.slow
@pytest.mark.parametrize(
    'standard_pars',
    standards
)
def test_against_standard(standard_pars):
    """Compare a generated assocaition against a standard
    """
    if standard_pars.xfail is not None:
        pytest.xfail(reason=standard_pars.xfail)

    generated, standards = generate_asns(standard_pars)
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
            assert False, '{}'.format(last_err)


def generate_asns(standard):
    """Test exp_type inclusion based on standard associations"""
    standards_paths = glob(t_path(path.join(
        'data',
        'asn_standards',
        standard.pool_root + '*_asn.json'))
    )
    standards = []
    for standard_path in standards_paths:
        with open(standard_path) as fp:
            asn = load_asn(fp)
        standards.append(asn)

    pool_path = t_path(path.join('data', standard.pool_root + '.csv'))
    pool = combine_pools([pool_path])
    args = TEST_ARGS + standard.main_args
    results = Main(args, pool=pool)

    asns = results.associations
    assert len(asns) == len(standards)
    return asns, standards
