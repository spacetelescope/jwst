"""
Contains a number of testing utilities.
"""

import contextlib
import logging
import os
import re


class ListHandler(logging.Handler):
    def __init__(self):
        self.records = []
        logging.Handler.__init__(self)

    def emit(self, record):
        self.records.append(record)


@contextlib.contextmanager
def capture_log():
    """
    Captures all LogRecords to a list, and returns the list::

       with capture_log() as log:
           # ...do something...

       # log is a list of all the captured LogRecords
    """
    handler = ListHandler()
    log = logging.getLogger('stpipe')
    log.addHandler(handler)
    yield handler.records
    log.removeHandler(handler)


def pattern_to_re(pattern):
    """
    Converts a pattern containing embedded regular expressions inside
    {{ }} to a Python regular expression.
    """
    regex = []
    while pattern:
        verbatim, sep, pattern = pattern.partition('{{')
        regex.append(re.escape(verbatim))
        exp, sep, pattern = pattern.partition('}}')
        regex.append(exp)
    return '^{0}$'.format(''.join(regex))


def match_log(log, expected):
    """
    Matches a log to an expected log.

    Parameters
    ----------
    log : list of LogRecord objects

    expected : list of strings
        Each string may contain embedded regular expressions inside of
        {{ }}.  For example, to match on any filename::

            "Opened file: {{.*}}."

    Raises
    ------
    ValueError : when one of the entries doesn't match
    """
    for a, b in zip(log, expected):
        msg = a
        regex = pattern_to_re(b)
        match = re.match(regex, msg)
        if not match:
            with open("correct.txt", 'w') as fd:
                fd.write('[\n')
                for a in log:
                    fd.write('    {0!r},\n'.format(a))
                fd.write(']\n')

            raise ValueError((
                "Logs do not match."
                "\nExpected:"
                "\n   '{0}'"
                "\nGot:"
                "\n   '{1}'\n".format(
                    b, msg
                )))


def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.dirname(__file__)
    return os.path.join(test_dir, partial_path)
