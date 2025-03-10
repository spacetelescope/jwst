"""Handy helpful pytest helpers helping pytest test."""

from os import path as op
from pathlib import Path


def abspath(filepath):
    """
    Get the absolute file path.

    Parameters
    ----------
    filepath: str
        The file path to get the absolute path for.

    Returns
    -------
    str
        The absolute file path.
    """
    return Path(op.expandvars(filepath)).expanduser().resolve()


# Check strings based on words using length precision
def word_precision_check(str1, str2, length=5):
    """
    Check two strings word-by-word based for word length.

    The strings are checked word for word, but only for the first `length` characters

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
    for w1, w2 in zip(words1, words2, strict=False):
        if w1[:length] != w2[:length]:
            break
    else:
        return True
    return False


class LogWatcher:
    """
    Watch a logger for a specific message.

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
        """
        Initialize the watcher with a message to watch for.

        Parameters
        ----------
        message : str
            The message of interest
        """
        self.seen = False
        self._message = message

    def __call__(self, *args):
        """Watch the logs for a specific message."""
        if not args or not isinstance(args[0], str):
            return
        if self.message in args[0]:
            self.seen = True

    @property
    def message(self):
        """
        str: The message to watch for.

        When the message is set, the `seen` flag is set to False.
        """
        return self._message

    @message.setter
    def message(self, new_message):
        self.seen = False
        self._message = new_message

    def assert_seen(self):
        """
        Check if message has been seen.

        After calling, the `seen` attribute is reset to False.
        """
        assert self.seen, f"{self.message} not in logs"  # noqa: S101

        # reset flag after check
        self.seen = False

    def assert_not_seen(self):
        """
        Check if message has not been seen.

        After calling, the `seen` attribute is reset to False.
        """
        assert not self.seen, f"{self.message} is in logs"  # noqa: S101

        # reset flag after check
        self.seen = False
