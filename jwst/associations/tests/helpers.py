"""Helpers for tests."""
from collections import namedtuple
import os

from astropy.table import (Table, vstack)

from .. import (AssociationRegistry, AssociationPool, generate)
from ..association import is_iterable

# Define how to setup initial conditions with pools.
class PoolParams(namedtuple('PoolParams', [
                            'path',
                            'n_asns',
                            'n_orphaned',
                            'candidates',
                            'kwargs'])):
    def __new__(cls, path='',
                n_asns=0,
                n_orphaned=0,
                candidates=None,
                kwargs=None):
        if not kwargs:
            kwargs = {}
        if candidates is None:
            candidates = []
        return super(PoolParams, cls).__new__(
            cls,
            path,
            n_asns,
            n_orphaned,
            candidates,
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
            pool = combine_pools(ppars.path, **ppars.kwargs)
            (asns, orphaned) = generate(pool, rules)
            yield check_equal, len(asns), ppars.n_asns
            yield check_equal, len(orphaned), ppars.n_orphaned
            for asn, candidates in zip(asns, ppars.candidates):
                yield check_equal, set(asn.candidates), set(candidates)


# Basic utilities.
def check_in_list(element, alist):
    assert element in alist


def check_not_in_list(element, alist):
    assert element not in alist


def check_equal(left, right):
    assert left == right


def not_none(value):
    assert value is not None


def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.dirname(__file__)
    return os.path.join(test_dir, partial_path)


def combine_pools(pools, **pool_kwargs):
    """Combine pools into a single pool

    Parameters
    ----------
    pools: str, astropy.table.Table, [str|Table, ...]
        The pools to combine. Either a singleton is
        passed or and iterable can be passed.
        The entries themselves can be either a file path
        or an astropy.table.Table-like object.

    pool_kwargs: dict
        Other keywoard arguments to pass to AssociationPool.read

    Returns
    -------
    AssociationPool|astropy.table.Table
        The combined pool
    """
    if not is_iterable(pools):
        pools = [pools]
    just_pools = []
    for pool in pools:
        if not isinstance(pool, Table):
            pool = AssociationPool.read(pool, **pool_kwargs)
        just_pools.append(pool)
    mega_pool = vstack(just_pools)
    return mega_pool
