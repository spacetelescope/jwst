"""Helpers for tests."""
from collections import namedtuple

from jwst.associations.association import AssociationRegistry
from jwst.associations.pool import AssociationPool
from jwst.associations.generate import generate


# Define how to setup initial conditions with pools.
class PoolParams(namedtuple('PoolParams', [
                            'path',
                            'n_asns',
                            'n_orphaned',
                            'n_candidates',
                            'kwargs'])):
    def __new__(cls, path='',
                n_asns=0,
                n_orphaned=0,
                n_candidates=None,
                kwargs=None):
        if not kwargs:
            kwargs = {}
        if n_candidates is None:
            n_candidates = []
        return super(PoolParams, cls).__new__(
            cls,
            path,
            n_asns,
            n_orphaned,
            n_candidates,
            kwargs
        )


# Basic Pool/Rule test class
class BasePoolRule(object):

    # Define the pools and testing parameters related to them.
    # Each entry is a tuple starting with the path of the pool.
    pools = []

    # Define the rules that SHOULD be present.
    # Each entry is the class name of the rule.
    valid_rules = []

    def setup(self):
        """You set me up...."""

    def tearDown(self):
        """You tear me down..."""

    def test_rules_exist(self):
        rules = AssociationRegistry()
        assert len(rules) >= len(self.valid_rules)
        for rule in self.valid_rules:
            yield check_in_list, rule, rules

    def test_run_generate(self):
        rules = AssociationRegistry()
        for ppars in self.pools:
            pool = AssociationPool.read(ppars.path, **ppars.kwargs)
            (asns, orphaned) = generate(pool, rules)
            yield check_equal, len(asns), ppars.n_asns
            yield check_equal, len(orphaned), ppars.n_orphaned
            for asn, n_candidates in zip(asns, ppars.n_candidates):
                yield check_equal, len(asn.candidates), n_candidates


# Basic utilities.
def check_in_list(element, alist):
    assert element in alist


def check_not_in_list(element, alist):
    assert element not in alist


def check_equal(left, right):
    assert left == right


def not_none(value):
    assert value is not None
