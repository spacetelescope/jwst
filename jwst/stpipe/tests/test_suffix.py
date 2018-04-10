"""Test suffix replacement"""

import pytest

from ..step import remove_suffix
from .. import suffix


def test_suffix_existence():
    """Generate current suffix list and compare"""
    calculated_suffixes = suffix.find_suffixes()
    found_suffixes = suffix.combine_suffixes(
        to_add=(calculated_suffixes, suffix.SUFFIXES_TO_ADD),
        to_remove=(suffix.SUFFIXES_TO_DISCARD, )
    )
    assert set(found_suffixes) == set(suffix.KNOW_SUFFIXES)


@pytest.mark.parametrize(
    'suffix',
    suffix.KNOW_SUFFIXES
)
def test_suffix_removal(suffix):
    """Test suffix removal"""
    basename = 'file'
    full_fpath = basename + '_' + suffix
    removed_path, separator = remove_suffix(full_fpath)
    assert removed_path == basename
    assert separator == '_'
