"""Handy helpful pytest helpers helping pytest test"""

import os
from os import path
import pytest
import re

import crds


def abspath(filepath):
    """Get the absolute file path"""
    return path.abspath(path.expanduser(path.expandvars(filepath)))


# Decorator to indicate slow test
runslow = pytest.mark.skipif(
    not pytest.config.getoption("--slow"),
    reason="need --slow option to run"
)


# Decorator to indicate TEST_BIGDATA required
require_bigdata = pytest.mark.skipif(
    not pytest.config.getoption('bigdata'),
    reason='requires --bigdata'
)


# Decorator to skip test if running under a TravisCI
not_under_travis = pytest.mark.skipif(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason='Temporarily disable due to performance issues'
)


# Decorator to skip if CRDS_CONTEXT is not at lest a certain level.
def require_crds_context(required_context):
    """Ensure CRDS context is a certain level

    Parameters
    ----------
    level: int
        The minimal level required

    Returns
    -------
    pytest.mark.skipif decorator
    """
    current_context_string = crds.get_context_name('jwst')
    match = re.match('jwst_(\d\d\d\d)\.pmap', current_context_string)
    current_context = int(match.group(1))
    return pytest.mark.skipif(
        current_context < required_context,
        reason='CRDS context {} less than required context {}'.format(
            current_context_string, required_context
        )
    )


# Check strings based on words using length precision
def word_precision_check(str1, str2, length=5):
    """Check to strings word-by-word based for word length

    The strings are checked word for word, but only for the first
    `length` characters

    Parameters
    ----------
    str1, str2: str
        The strings to compare

    length: int
        The number of characters in each word to check.

    Returns
    -------
    match: boolean
        True if the strings match
    """
    words1 = str1.split()
    words2 = str2.split()
    if len(words1) != len(words2):
        return False
    for w1, w2 in zip(words1, words2):
        if w1[:length] != w2[:length]:
            break
    else:
        return True
    return False


def test_word_precision_check():
    """Test word_precision_check"""
    s1 = "a b c"
    s2 = "aa bb cc"
    s3 = "aa bb cc dd"
    s4 = "aazz bbzz cczz"

    assert word_precision_check(s1, s1)
    assert not word_precision_check(s1, s2)
    assert word_precision_check(s1, s2, length=1)
    assert not word_precision_check(s2, s3)
    assert word_precision_check(s2, s4, length=2)
