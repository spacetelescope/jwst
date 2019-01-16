"""Test using SDP-generated pools

Notes
-----
Most of the standard associations which are compared
against are built in the jupyter notebook

./notebooks/make_tests.ipynb
"""
from pathlib import Path

import pytest

from jwst.associations.lib.diff import (
    compare_asn_files,
)
from jwst.lib.file_utils import pushdir
from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations.main import Main

# Main test args
TEST_ARGS = ['--dry-run', '--no-merge']


class AssociationBase(BaseJWSTTest):
    input_loc = 'associations'
    test_dir = 'sdp'
    ref_loc = [test_dir, 'truth']

    @property
    def pool_paths(self):
        return self.data_glob(self.test_dir, 'pools', glob='*.csv')


class TestAgainstStandards(AssociationBase):
    @pytest.mark.parametrize(
        'pool_path',
        AssociationBase().pool_paths
    )
    def test_against_standard(self, pool_path):
        """Compare a generated association against a standard
        """

        # Create the associations
        generated_path = Path('generate')
        generated_path.mkdir()
        local_pool_path = self.get_data(pool_path)
        Main([
            '--no-merge',
            '-p', str(generated_path),
            local_pool_path
        ])

        # Retrieve the truth files
        pool = Path(pool_path).stem
        sdp_path = Path('sdp')
        sdp_path.mkdir()
        with pushdir(sdp_path):
            for standard_path in self.data_glob(*self.ref_loc, glob=pool + '*.json'):
                self.get_data(standard_path)

        # Compare the association sets.
        assert compare_asn_files(
            generated_path.glob('*.json'), standard_path.glob('*.json')
        )
