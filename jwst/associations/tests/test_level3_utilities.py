import pytest

from . import helpers
from .helpers import full_pool_rules

from .. import (AssociationRegistry, AssociationPool, generate)


class TestUtilities():
    @pytest.mark.xfail
    def test_filter_cross_candidates(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules
        (asns, orphaned) = generate(pool, rules)
        assert len(asns) == 11
        assert len(orphaned) == 298
        filtered = rules.Utility.filter_cross_candidates(asns)
        assert len(filtered) == 5
