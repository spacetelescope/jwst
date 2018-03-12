"""Test for duplication/missing associations"""

from .helpers import (
    level3_rule_path,
    registry_level3_only,
    t_path,
)

from .. import (AssociationPool, generate)
from ..main import (Main, constrain_on_candidates)


def test_duplicate_generate():
    """Test for duplicate/overwrite association

    The pool has two exposures, one without a valid `asn_candidate`,
    and one with a valid observation `asn_candidate`.
    When set with the "all candidates" constraint, only one association
    should be made.

    The prompt for this test was that three associations were being created,
    two of which were the observation candidate, with the second
    being a duplicate of the first. The third was an extraneous
    discovered candidate.
    """
    pool = AssociationPool.read(t_path('data/pool_duplicate.csv'))
    constrain_all_candidates = constrain_on_candidates(None)
    rules = registry_level3_only(global_constraints=constrain_all_candidates)
    asns = generate(pool, rules)
    assert len(asns) == 1
    asn = asns[0]
    assert asn['asn_type'] == 'image3'
    assert asn['asn_id'] == 'o029'


def test_duplicate_main():
    """Test the same but with Main
    """
    cmd_args = [
        t_path('data/pool_duplicate.csv'),
        '--dry-run',
        '-r', level3_rule_path(),
        '--ignore-default'
    ]

    generated = Main(cmd_args)
    asns = generated.associations
    assert len(asns) == 2
    asn_types = set([
        asn['asn_type']
        for asn in asns
    ])
    assert len(asn_types) == 1
    assert 'image3' in asn_types
    asn_ids = set([
        asn['asn_id']
        for asn in asns
    ])
    assert asn_ids == set(('a3001', 'o029'))
