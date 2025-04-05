"""Testing utilities for steps and pipelines."""

import contextlib
import logging
from pathlib import Path
import re


class ListHandler(logging.Handler):
    """Put log messages into a list for testing."""

    def __init__(self):
        self.records = []
        logging.Handler.__init__(self)

    def emit(self, record):
        """
        Emit a record.

        Parameters
        ----------
        record : LogRecord
            The record to append to the records list.
        """
        self.records.append(record)


@contextlib.contextmanager
def capture_log():
    """
    Capture all LogRecords to a list, and return the list.

       with capture_log() as log:
           # ...do something...

    log is a list of all the captured LogRecords

    Yields
    ------
    list
        List of LogRecord objects
    """
    handler = ListHandler()
    log = logging.getLogger("stpipe")
    log.addHandler(handler)
    yield handler.records
    log.removeHandler(handler)


def pattern_to_re(pattern):
    """
    Convert a pattern containing embedded regex inside {{ }} to a Python regular expression.

    Parameters
    ----------
    pattern : str
        A pattern that may contain embedded regular expressions inside of {{ }}.

    Returns
    -------
    str
        A regular expression that matches the pattern.
    """
    regex = []
    while pattern:
        verbatim, sep, pattern = pattern.partition("{{")
        regex.append(re.escape(verbatim))
        exp, sep, pattern = pattern.partition("}}")
        regex.append(exp)
    return "^{}$".format("".join(regex))


def match_log(log, expected):
    """
    Match a log to an expected log.

    Parameters
    ----------
    log : list of LogRecord objects
        The log to match.

    expected : list of strings
        Each string may contain embedded regular expressions inside of
        {{ }}.  For example, to match on any filename::

            "Opened file: {{.*}}."

    Raises
    ------
    ValueError : when one of the entries doesn't match
    """
    for a, b in zip(log, expected, strict=False):
        msg = a
        regex = pattern_to_re(b)
        match = re.match(regex, msg)
        if not match:
            with Path("correct.txt").open("w") as fd:
                fd.write("[\n")
                for a in log:
                    fd.write(f"    {a!r},\n")
                fd.write("]\n")

            raise ValueError(f"Logs do not match.\nExpected:\n   '{b}'\nGot:\n   '{msg}'\n")
