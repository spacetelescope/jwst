from __future__ import absolute_import

from glob import glob
import os

from .helpers import (
    TemporaryDirectory,
    full_pool_rules,
)

from ..main import Main
from .. import Association


class TestIO(object):

    def test_roundtrip(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules
        with TemporaryDirectory() as path:
            generated = Main([pool_fname, '-p', path])

            asn_files = glob(os.path.join(path, '*.json'))
            assert len(asn_files) == len(generated.associations)

            for asn_file in asn_files:
                with open(asn_file, 'r') as asn_fp:
                    asn = Association.load(asn_fp)
                valid_schemas = generated.rules.validate(asn)
                assert isinstance(valid_schemas, list)
