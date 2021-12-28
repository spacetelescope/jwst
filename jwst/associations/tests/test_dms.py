"""Test DMSBaseMixin features"""
import inspect
from os import path
import pytest
import re

from jwst.associations.tests.helpers import (
    t_path
)

from jwst.associations import AssociationRegistry


@pytest.fixture(scope='module')
def dms_registry():
    """Create the registry"""
    dms_test_rules_path = t_path(path.join('data', 'dms_rules.py'))
    dms_registry = AssociationRegistry(
        [dms_test_rules_path], include_default=False
    )
    return dms_registry


@pytest.fixture(scope='module')
def dms_asns(dms_registry):
    """Create basic associations"""
    result = dms_registry.match('item')
    return result


def test_asn_name(dms_asns):
    """Test for generating an association name"""
    asns, _ = dms_asns
    asn = asns[0]

    regex = re.compile(r'jwnoprogram-a3001_none_\d\d\d_asn')
    assert regex.match(asn.asn_name)


def test_asn_name_override(dms_asns):
    """Test for generating an association name"""
    asns, _ = dms_asns
    asn = asns[0]
    asn.asn_name = 'new_name'
    assert asn.asn_name == 'new_name'


def test_registry(dms_registry):
    """Test basic registry creation and usage"""
    assert len(dms_registry) == 1
    assert 'Asn_DMS_Base' in dms_registry


def test_asn(dms_asns):
    """Test basic association creation"""
    asns, orphaned = dms_asns
    assert len(asns) == 1
    assert len(orphaned) == 0


def test_finalize(dms_registry, dms_asns):
    """Test finalization"""
    asns, orphaned = dms_asns

    finalized = dms_registry.callback.reduce('finalize', asns)

    assert len(finalized) == 1


def test_include_bases():
    """Test for included bases"""
    dms_test_rules_path = t_path(path.join('data', 'dms_rules.py'))
    dms_registry = AssociationRegistry(
        [dms_test_rules_path], include_default=False, include_bases=True
    )

    assert len(dms_registry) > 1
    assert 'DMSBaseMixin' in dms_registry


def test_utility(dms_registry):
    """Test the utility inclusion marker"""

    names, objs = zip(*inspect.getmembers(dms_registry.Utility))
    assert 'valid_function' in names
