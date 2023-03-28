"""Test using SDP-generated pools
"""
import logging
from pathlib import Path
import re

import pytest

from jwst.associations.lib.diff import (
    compare_asn_files,
)
from jwst.associations.main import Main as asn_generate
from jwst.lib.file_utils import pushdir

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Decompose pool name to retrieve proposal and version id.
pool_regex = re.compile(r'(?P<proposal>jw.+?)_(?P<versionid>.+)_pool')


# Mark expected failures. Key is the pool name
# and value is the reason message.
EXPECTED_FAILS = {
}

# Pools that require special handling
SPECIAL_DEFAULT = {
    'args': [],
    'xfail': None,
    'slow': False,
}
SPECIAL_POOLS = {
    'jw00623_20190607t021101_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw00628_20191102t153956_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw00629_20190605t025157_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw00632_20190819t190005_pool': {
        'args': [],
        'xfail': 'JSOCINT-TDB: WFSC ROUTINE VISIT issue',
        'slow': False,
    },
    'jw00663_20221218t111937_pool': {
        'args': ['-i', 'o004', 'c1000'],
        'xfail': None,
        'slow': False,
    },
    'jw00676_20210403t114320_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw00675_20211225t181823_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw00676_20210403t114320_c1007_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw00839_20221220t025418_pool': {
        'args': ['-i', 'o002', 'c1000'],
        'xfail': None,
        'slow': False,
    },
    'jw01194_20230115t113819_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw01257_20221201t192226_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw01290_20230304t140931_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw01288_c1005_mostilno12_pool': {
        'args': ['-i', 'o003', 'c1001', 'c1005'],
        'xfail': None,
        'slow': True,
    },
    'jw01290_20230304t140931_withids_pool': {
        'args': ['-i', 'o012', 'c1018'],
        'xfail': None,
        'slow': False,
    },
    'jw01355_20230109t002554_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw01493_20230307t040130_withids_pool': {
        'args': ['-i', 'o003', 'c1000'],
        'xfail': None,
        'slow': False,
    },
    'jw02064_20230302t112350_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw02064_20230302t112350_withids_pool': {
        'args': ['-i', 'o061', 'c1008', 'c1017'],
        'xfail': None,
        'slow': False,
    },
    'jw80600_20171108T041522_pool': {
        'args': [],
        'xfail': 'PR #3450',
        'slow': False,
    },
    'jw82600_20180921T023255_pool': {
        'args': [],
        'xfail': None,
        'slow': True
    },
    'jw93065_20171108T041402_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw93135_20171108T041617_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw93135_20171108T041617-fixed_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw94015_20171108T041516_pool': {
        'args': [],
        'xfail': None,
        'slow': True,
    },
    'jw98010_20171108T062332_pool': {
        'args': [],
        'xfail': 'PR #3450',
        'slow': False,
    },
}


# #####
# Tests
# #####
@pytest.mark.filterwarnings('error')
def test_against_standard(sdpdata_module, pool_path, slow):
    """Compare a generated association against a standard

    Success is when no other AssertionError occurs.
    """

    # Parse pool name
    pool = Path(pool_path).stem
    proposal, version_id = pool_regex.match(pool).group('proposal', 'versionid')
    special = SPECIAL_POOLS.get(pool, SPECIAL_DEFAULT)

    if special['slow'] and not slow:
        pytest.skip(f'Pool {pool} requires "--slow" option')

    # Setup test path
    cwd = Path(pool)
    cwd.mkdir()
    with pushdir(cwd):

        # Create the generator running arguments
        output_path = Path(pool)
        output_path.mkdir()
        sdpdata_module.output = str(output_path)
        args = special['args'] + [
            '-p', sdpdata_module.output,
            '--version-id', version_id,
            sdpdata_module.get_data(pool_path)
        ]

        # Create the associations
        asn_generate.cli(
            args)

        # Compare to the truth associations.
        truth_paths = sdpdata_module.truth_paths(pool)
        try:
            compare_asn_files(output_path.glob('*.json'), truth_paths)
        except AssertionError:
            if special['xfail']:
                pytest.xfail(special['xfail'])
            else:
                raise
