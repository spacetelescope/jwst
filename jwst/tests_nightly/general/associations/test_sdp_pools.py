"""Test using SDP-generated pools

Notes
-----
Most of the standard associations which are compared
against are built in the jupyter notebook

./notebooks/make_tests.ipynb
"""
from pathlib import Path
import pytest
import re

from jwst.associations.lib.diff import (
    compare_asn_files,
)
from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations.main import Main

# Main test args
TEST_ARGS = ['--dry-run', '--no-merge']

# Get version id out of pool name
version_id_regex = re.compile(r'(?P<proposal>jw.+?)_(?P<versionid>.+)_pool')


# #############################################################
# Setup a base class and instantiate it in order to provide the
# file lists for the test parametrization.
# #############################################################
class AssociationBase(BaseJWSTTest):
    input_loc = 'associations'
    test_dir = 'sdp'
    ref_loc = [test_dir, 'truth']

    _pool_paths = None
    _truth_paths = None

    @property
    def pool_paths(self):
        if self._pool_paths is None:
            self._pool_paths = self.data_glob(self.test_dir, 'pools', glob='*.csv')
        return self._pool_paths

    @property
    def truth_paths(self):
        if self._truth_paths is None:
            self._truth_paths = self.data_glob(*self.ref_loc, glob='*.json')
        return self._truth_paths

asn_base = AssociationBase()


# #####
# Tests
# #####
class TestAgainstStandards(AssociationBase):
    @pytest.mark.parametrize(
        'pool_path',
        asn_base.pool_paths
    )
    def test_against_standard(self, pool_path):
        """Compare a generated association against a standard

        Success is when no other AssertionError occurs.
        """

        # Parse pool name
        pool = Path(pool_path).stem
        proposal, version_id = version_id_regex.match(pool).group('proposal', 'versionid')

        # Create the associations
        generated_path = Path('generate')
        generated_path.mkdir()
        Main([

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
            for truth_path in asn_base.truth_paths
            if asn_regex.match(truth_path)
        ]

        # Compare the association sets.
        try:
            compare_asn_files(generated_path.glob('*.json'), truth_paths)
        except AssertionError as error:
            if 'Associations do not share a common set of products' in str(error):
                pytest.xfail('Possibly due to Issue #3039 but manually confirm.')
            else:
                raise
