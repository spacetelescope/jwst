import pytest

from .. import (
    AssociationPool,
    AssociationRegistry,
    )
from .helpers import t_path

@pytest.fixture(scope='session')
def full_pool_rules(request):
    pool_fname = t_path('data/mega_pool.csv')
    pool = AssociationPool.read(pool_fname)
    rules = AssociationRegistry()

    return (pool, rules, pool_fname)
