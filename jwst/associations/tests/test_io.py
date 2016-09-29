from __future__ import absolute_import

from glob import glob
import os
import pytest

from astropy.table import Table

from .helpers import (
    TemporaryDirectory,
    full_pool_rules,
)

from ..main import Main
from .. import Association


@pytest.yield_fixture(
    scope='module',
    params=['yaml', 'json']
)
def make_asns(request):
    asn_format = request.param
    pool, rules, pool_fname = full_pool_rules(None)
    with TemporaryDirectory() as path:
        generated = Main([
            pool_fname,
            '-p', path,
            '--save-orphans',
            '--format', asn_format
        ])
        yield generated, path, asn_format


class TestIO(object):

    def test_roundtrip(self, make_asns):
        generated, path, asn_format = make_asns
        asn_files = glob(os.path.join(path, '*.' + asn_format))
        assert len(asn_files) == len(generated.associations)

        for asn_file in asn_files:
            with open(asn_file, 'r') as asn_fp:
                asn = Association.load(asn_fp)
            valid_schemas = generated.rules.validate(asn)
            assert isinstance(valid_schemas, list)

        orphaned_files = glob(os.path.join(path, '*.csv'))
        assert len(orphaned_files) == 1
        orphaned = Table.read(
            orphaned_files[0],
            format='ascii',
            delimiter='|'
        )
        assert len(orphaned) == len(generated.orphaned)
