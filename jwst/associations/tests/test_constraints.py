"""Constraint Tests"""
import pytest

from ..lib.constraint import (
    Constraint,
    SimpleConstraint,
    SimpleConstraintABC,
)


def test_sc_dup_names():
    """Test that SimpleConstraint returns an empty dict"""
    sc = SimpleConstraint(name='sc_name')
    dups = sc.dup_names
    assert not len(dups)

class TestDupNames:
    """Test duplicate names in Constraint"""

    # Sub constraints
    sc1 = SimpleConstraint(name='sc1')
    sc2 = SimpleConstraint(name='sc2')
    c1 = Constraint([sc1, sc2], name='c1')
    c2 = Constraint([sc1, sc1], name='c2')
    c3 = Constraint([sc1, sc2], name='sc1')

    @pytest.mark.parametrize(
        'constraints, expected', [
            ([sc1], {}),
            ([sc1, sc2], {}),
            ([sc1, sc1], {'sc1': [sc1, sc1]}),
            ([c1], {}),
            ([c2], {'sc1': [sc1, sc1]}),
            ([c3], {'sc1': [sc1, c3]}),
        ]
    )
    def test_dups(self, constraints, expected):
        c = Constraint(constraints)
        dups = c.dup_names
        assert set(dups.keys()) == set(expected.keys())
        for name, constraints in dups.items():
            assert set(dups[name]) == set(expected[name])


def test_sc_get_all_attr():
    """Get attribute value of a simple constraint"""
    name = 'my_sc'
    sc = SimpleConstraint(name=name, value='my_value')
    assert sc.get_all_attr('name') == [(sc, name)]


def test_constraint_get_all_attr():
    """Get attribute value of all constraints in a constraint"""
    names = ['sc1', 'sc2']
    constraints = [
        SimpleConstraint(name=name)
        for name in names
    ]
    c = Constraint(constraints, name='c1')

    expected = [
        (constraint, constraint.name)
        for constraint in constraints
    ]
    expected.append((c, 'c1'))
    result = c.get_all_attr('name')
    assert set(result) == set(expected)


def test_simpleconstraint_reprocess_match():
    """Test options for reprocessing"""
    sc = SimpleConstraint(
        value='my_value',
        reprocess_on_match=True
    )
    match, reprocess = sc.check_and_set('my_value')
    assert match
    assert len(reprocess)


def test_simpleconstraint_reprocess_nomatch():
    """Test options for reprocessing"""
    sc = SimpleConstraint(
        value='my_value',
        reprocess_on_fail=True
    )
    match, reprocess = sc.check_and_set('bad_value')
    assert not match
    assert len(reprocess)


def test_constraint_reprocess_match():
    """Test options for reprocessing"""
    sc = SimpleConstraint(value='my_value')
    c = Constraint([sc], reprocess_on_match=True)
    match, reprocess = c.check_and_set('my_value')
    assert match
    assert len(reprocess)


def test_constraint_reprocess_nomatch():
    """Test options for reprocessing"""
    sc = SimpleConstraint(value='my_value')
    c = Constraint([sc], reprocess_on_fail=True)
    match, reprocess = c.check_and_set('bad_value')
    assert not match
    assert len(reprocess)


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
    match, reprocess = c.check_and_set('my_value')
    assert match
    assert c.value == 'my_value'
    assert len(reprocess) == 0

    # Non-match
    c = SimpleConstraint(value='my_value')
    match, reprocess = c.check_and_set('bad_value')
    assert not match
    assert c.value == 'my_value'
    assert len(reprocess) == 0

    # Don't force unique
    c = SimpleConstraint(force_unique=False)
    match, reprocess = c.check_and_set('my_value')
    assert match
    assert c.value is None
    assert len(reprocess) == 0


def test_constraint_default():
    """Test constraint operations"""

    sc1 = SimpleConstraint()
    sc2 = SimpleConstraint()
    c = Constraint([sc1, sc2])
    match, reprocess = c.check_and_set('my_value')
    assert match
    for constraint in c.constraints:
        assert constraint.value == 'my_value'


def test_invalid_init():
    with pytest.raises(TypeError):
        Constraint('bad init')


def test_constraint_all():
    """Test the all operation"""

    sc1 = SimpleConstraint(value='value_1')
    sc2 = SimpleConstraint(value='value_2')
    c = Constraint([sc1, sc2])
    match, reprocess = c.check_and_set('value_1')
    assert not match


def test_constraint_any_basic():
    """Test the all operation"""

    sc1 = SimpleConstraint(value='value_1')
    sc2 = SimpleConstraint(value='value_2')
    c = Constraint([sc1, sc2], reduce=Constraint.any)
    match, reprocess = c.check_and_set('value_1')
    assert match
    match, reprocess = c.check_and_set('value_2')
    assert match
    match, reprocess = c.check_and_set('value_3')
    assert not match


def test_constraint_any_remember():
    """Ensure that any doesn't forget other or propositions"""

    sc1 = SimpleConstraint(value='value_1')
    sc2 = SimpleConstraint(value='value_2')
    c = Constraint([sc1, sc2], reduce=Constraint.any)
    match, reprocess = c.check_and_set('value_1')
    assert match
    match, reprocess = c.check_and_set('value_2')
    assert match
    match, reprocess = c.check_and_set('value_1')
    assert match
    match, reprocess = c.check_and_set('value_3')
    assert not match


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


def test_name_index():
    """Test for name indexing"""
    sc1 = SimpleConstraint(name='sc1', value='value1')
    sc2 = SimpleConstraint(name='sc2', value='value2')
    c1 = Constraint([sc1, sc2])
    assert c1['sc1'].value
    assert c1['sc2'].value

    sc3 = SimpleConstraint(name='sc3', value='value3')
    sc4 = SimpleConstraint(name='sc4', value='value4')
    c2 = Constraint([sc3, sc4, c1])
    assert c2['sc1'].value
    assert c2['sc2'].value
    assert c2['sc3'].value
    assert c2['sc4'].value

    with pytest.raises(KeyError):
        c2['nonexistant'].value

    with pytest.raises(AttributeError):
        c2['sc1'].nonexistant


def test_copy():
    sc1 = SimpleConstraint(name='sc1')
    sc1_copy = sc1.copy()
    assert id(sc1) != id(sc1_copy)
    sc1.check_and_set('value1')
    assert sc1.value == 'value1'
    assert sc1_copy.value is None
    sc1_copy.check_and_set('value2')
    assert sc1_copy.value == 'value2'
    assert sc1.value == 'value1'
