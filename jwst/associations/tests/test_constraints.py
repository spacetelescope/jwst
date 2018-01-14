"""Constraint Tests"""
import pytest

from ..lib.constraint import (
    Constraint,
    SimpleConstraint,
    SimpleConstraintABC,
)


def test_abc():
    """Test ABC istelf"""
    with pytest.raises(TypeError):
        SimpleConstraintABC()


def test_simpleconstraint():
    """Test initialization"""

    # Basic initialization
    c = SimpleConstraint()
    assert c.value is None
    assert c.force_unique
    assert c.test == c.eq
    assert c._reprocess == []

    # Parameter initialization
    c = SimpleConstraint(value='my_value')
    assert c.value == 'my_value'

    # Dict initialization
    c = SimpleConstraint({'value': 'my_value'})
    assert c.value == 'my_value'


def test_simpleconstraint_checkset():
    """Test check_and_set"""

    # Check and set.
    c = SimpleConstraint()
    new_c, reprocess = c.check_and_set('my_value')
    assert c != new_c
    assert new_c.value == 'my_value'
    assert len(reprocess) == 0

    # Non-match
    c = SimpleConstraint(value='my_value')
    new_c, reprocess = c.check_and_set('bad_value')
    assert not new_c
    assert len(reprocess) == 0

    # Don't force unique
    c = SimpleConstraint(force_unique=False)
    new_c, reprocess = c.check_and_set('my_value')
    assert c == new_c
    assert new_c.value is None
    assert len(reprocess) == 0


def test_constraint_default():
    """Test constraint operations"""

    sc1 = SimpleConstraint()
    sc2 = SimpleConstraint()
    c = Constraint([sc1, sc2])
    new_c, reprocess = c.check_and_set('my_value')
    assert new_c
    assert new_c != c
    for constraint in new_c.constraints:
        assert constraint.value == 'my_value'


def test_invalid_init():
    with pytest.raises(TypeError):
        Constraint('bad init')


def test_constraint_all():
    """Test the all operation"""

    sc1 = SimpleConstraint(value='value_1')
    sc2 = SimpleConstraint(value='value_2')
    c = Constraint([sc1, sc2])
    new_c, reprocess = c.check_and_set('value_1')
    assert not new_c


def test_constraint_any_basic():
    """Test the all operation"""

    sc1 = SimpleConstraint(value='value_1')
    sc2 = SimpleConstraint(value='value_2')
    c = Constraint([sc1, sc2], reduce=Constraint.any)
    new_c, reprocess = c.check_and_set('value_1')
    assert new_c
    new_c, reprocess = c.check_and_set('value_2')
    assert new_c
    new_c, reprocess = c.check_and_set('value_3')
    assert not new_c


def test_constraint_any_remember():
    """Ensure that any doesn't forget other or propositions"""

    sc1 = SimpleConstraint(value='value_1')
    sc2 = SimpleConstraint(value='value_2')
    c = Constraint([sc1, sc2], reduce=Constraint.any)
    new_c, reprocess = c.check_and_set('value_1')
    assert new_c
    new_c, reprocess = new_c.check_and_set('value_2')
    assert new_c
    new_c, reprocess = new_c.check_and_set('value_1')
    assert new_c
    new_c, reprocess = new_c.check_and_set('value_3')
    assert not new_c


def test_iteration():
    """Test various iterations"""
    sc = SimpleConstraint()
    for idx in sc:
        assert isinstance(idx, SimpleConstraint)

    c = Constraint([sc, sc])
    count = 0
    for idx in c:
        assert isinstance(idx, SimpleConstraint)
        count += 1
    assert count == 2

    c = Constraint([
        Constraint([sc, sc]),
        Constraint([sc, sc])
    ])
    count = 0
    for idx in c:
        assert isinstance(idx, SimpleConstraint)
        count += 1
    assert count == 4  # Not 6
