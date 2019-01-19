"""Association test configuration"""
import pytest

from .. import (AssociationRegistry, AssociationPool)
from .helpers import t_path


@pytest.fixture(scope='session')
def full_pool_rules(request):
    """Setup to use the full example pool and registry"""
    pool_fname = t_path('data/mega_pool.csv')
    pool = AssociationPool.read(pool_fname)
    rules = AssociationRegistry()

    return (pool, rules, pool_fname)
