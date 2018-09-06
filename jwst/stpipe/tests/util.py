# Copyright (C) 2010 Association of Universities for Research in Astronomy(AURA)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

"""
Contains a number of testing utilities.
"""

import contextlib
import logging
import os
import pytest
import re
import tempfile

from ...tests.helpers import mk_tmp_dirs


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


@contextlib.contextmanager
def get_tempfile():
    with tempfile.NamedTemporaryFile(delete=False) as fd:
        filename = fd.name

    yield filename

    os.remove(filename)
