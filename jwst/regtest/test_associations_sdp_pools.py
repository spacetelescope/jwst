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

from jwst.regtest.sdp_pools_source import SDPPoolsSource

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
class TestSDPPools(SDPPoolsSource):
    """Test creation of association from SDP-created pools"""
    def test_against_standard(self, pool_path, slow):
        """Compare a generated association against a standard

        Success is when no other AssertionError occurs.
        """

        # Parse pool name
        pool = Path(pool_path).stem
        proposal, version_id = pool_regex.match(pool).group('proposal', 'versionid')
        special = SPECIAL_POOLS.get(pool, SPECIAL_DEFAULT)

        if special['slow'] and not slow:
            pytest.skip('Pool {pool} requires "--slow" option')

        # Create the generator running arguments
        generated_path = Path('generate')
        generated_path.mkdir()
        args = special['args'] + [
            '-p', str(generated_path),
            '--version-id', version_id,
            self.get_data(pool_path)
        ]

        # Create the associations
        asn_generate(args)

        # Retrieve the truth files
        asn_regex = re.compile(
            r'.+{proposal}.+{version_id}(_[^_]+?_[^_]+?_asn\.json)$'.format(
                proposal=proposal, version_id=version_id
            ),
            flags=re.IGNORECASE
        )
        truth_paths = [
            self.get_data(truth_path)
            for truth_path in self.truth_paths
            if asn_regex.match(truth_path)
        ]

        # Compare the association sets.
        try:
            compare_asn_files(generated_path.glob('*.json'), truth_paths)
        except AssertionError:
            if special['xfail']:
                pytest.xfail(special['xfail'])
            else:
                raise
