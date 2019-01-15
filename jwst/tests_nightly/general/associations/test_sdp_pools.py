"""Test using SDP-generated pools

Notes
-----
Most of the standard associations which are compared
against are built in the jupyter notebook

./notebooks/make_tests.ipynb
"""
from pathlib import Path

import pytest

from jwst.associations.tests.helpers import (
    compare_asns,
)
from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations import load_asn
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
        local_pool_path = self.get_data(pool_path)
        args = TEST_ARGS + [local_pool_path]
        generated = Main(args).associations

        # Retrieve the truth files
        pool = Path(pool_path).stem
        standards = []
        for standard_path in self.data_glob(*self.ref_loc, glob=pool + '*.json'):
            local_path = self.get_data(standard_path)
            with open(local_path, 'r') as fh:
                standards.append(load_asn(fh))

        # Start assertions
        assert len(generated) == len(standards)
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
