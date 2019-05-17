"""Test using SDP-generated pools
"""
from collections import Counter
import logging
from pathlib import Path
import re

import pytest

from jwst.associations.lib.diff import (
    compare_asn_files,
)
from jwst.associations.main import Main as asn_generate

from jwst.tests_nightly.general.associations.sdp_pools_source import SDPPoolsSource

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Main test args
TEST_ARGS = ['--dry-run', '--no-merge']

# Decompose pool name to retrieve proposal and version id.
pool_regex = re.compile(r'(?P<proposal>jw.+?)_(?P<versionid>.+)_pool')


# Mark expected failures. Key is the pool name
# and value is the reason message.
EXPECTED_FAILS = {
    'jw00624_20190205t031003_pool': 'Issue #655',
    'jw80600_20171108T041522_pool': 'PR #3450',
    'jw87600_20180824T213416_pool': 'Issue #3039',
    'jw98010_20171108T062332_pool': 'PR #3450',
}


# #####
# Tests
# #####
class TestSDPPools(SDPPoolsSource):
    """Test createion of association from SDP-created pools"""
    def test_against_standard(self, pool_path):
        """Compare a generated association against a standard

        Success is when no other AssertionError occurs.
        """

        # Parse pool name
        pool = Path(pool_path).stem
        proposal, version_id = pool_regex.match(pool).group('proposal', 'versionid')

        # Create the associations
        generated_path = Path('generate')
        generated_path.mkdir()
        asn_generate([
            '--no-merge',
            '-p', str(generated_path),
            '--version-id', version_id,
            self.get_data(pool_path)
        ])

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
            if pool in EXPECTED_FAILS:
                pytest.xfail(EXPECTED_FAILS[pool])
            else:
                raise

    def test_dup_product_names(self, pool_path):
        """Check for duplicate product names for a pool"""

        results = asn_generate([
            '--dry-run',
            '--no-merge',
            self.get_data(pool_path)
        ])
        asns = results.associations

        product_names = Counter(
            product['name']
            for asn in asns
            for product in asn['products']
        )

        multiples = [
            product_name
            for product_name, count in product_names.items()
            if count > 1
        ]

        assert not multiples, 'Multiple product names: {}'.format(multiples)

    def test_specified_sdp_pool(self, sdp_pool):
        """Test a command-line specified pool"""
        if sdp_pool:
            pool_path = Path(self.test_dir) / 'pools' / (sdp_pool + '.csv')
            self.test_against_standard(pool_path)
        else:
            pytest.skip('No SDP pool specified using `--sdp-pool` command-line option.')
