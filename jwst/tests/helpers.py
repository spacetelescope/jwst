"""Handy helpful pytest helpers helping pytest test"""
from os import path as op


def abspath(filepath):
    """Get the absolute file path"""
    return op.abspath(op.expanduser(op.expandvars(filepath)))


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
