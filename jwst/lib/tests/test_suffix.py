"""Test suffix replacement"""

import logging
import pytest

from jwst.lib.basic_utils import LoggingContext

from .. import suffix as s


@pytest.fixture(scope='module')
def enable_logging():
    with LoggingContext(
            logging.getLogger(s.__name__),
            handler=logging.StreamHandler()
    ):
        yield


def test_suffix_existence(enable_logging):
    """Generate current suffix list and compare"""

    calculated_suffixes = s.find_suffixes()
    found_suffixes = s.combine_suffixes(
        to_add=(calculated_suffixes, s.SUFFIXES_TO_ADD),
        to_remove=(s.SUFFIXES_TO_DISCARD, )
    )
    diff = set(found_suffixes).symmetric_difference(s.KNOW_SUFFIXES)
    assert not diff, f'Suffixes unaccounted for:\n{diff}'


@pytest.mark.parametrize(
    'suffix',
    s.KNOW_SUFFIXES
)
def test_suffix_removal(suffix, enable_logging):
    """Test suffix removal"""
    basename = 'file'
    full_fpath = basename + '_' + suffix
    removed_path, separator = s.remove_suffix(full_fpath)
    assert removed_path == basename
    assert separator == '_'


@pytest.mark.parametrize(
    'suffix',
    s.KNOW_SUFFIXES
)
def test_suffix_replacement(suffix, enable_logging, base='file', new='junk', sep='_'):
    """Test suffix replacement"""
    full_path = base + sep + suffix
    replaced = s.replace_suffix(full_path, new)
    assert replaced == base + sep + new
