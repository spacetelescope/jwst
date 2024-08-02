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


class LogWatcher:
    """
    The pytest caplog fixture only works for the root
    logger. We use all sorts of loggers which can lead
    to random test failures with caplog.

    This class can be monkeypatched onto a logger method to
    check for a specific message. When the logger method is
    called, a `seen` attribute is set in the watcher,
    if the calling message matches the expected value.

    The `seen` flag can be checked via `assert_seen`. When
    called, this function will reset the `seen` attribute
    to False. This allows the same watcher to be reused for
    multiple function calls.
    """
    def __init__(self, message):
        self.seen = False
        self.message = message

    def __call__(self, *args, **kwargs):
        if not args or not isinstance(args[0], str):
            return
        if self.message in args[0]:
            self.seen = True

    def assert_seen(self):
        """Check if message has been seen.

        After calling, the `seen` attribute is reset to False.
        """
        assert self.seen, f"{self.message} not in logs"

        # reset flag after check
        self.seen = False
