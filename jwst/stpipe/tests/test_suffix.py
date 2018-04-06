"""Test suffix replacement"""

import pytest

from ..step import remove_suffix
from ..suffix import (KNOW_SUFFIXES, find_suffixes)


def test_suffix_existence():
    """Generate current suffix list and compare"""
    new_suffixes = find_suffixes()
    assert new_suffixes == KNOW_SUFFIXES


@pytest.mark.parametrize(
    'suffix',
    KNOW_SUFFIXES
)
def test_suffix_removal(suffix):
    """Test suffix removal"""
    basename = 'file'
    full_fpath = basename + '_' + suffix
    removed_path, separator = remove_suffix(full_fpath)
    assert removed_path == basename
    assert separator == '_'
