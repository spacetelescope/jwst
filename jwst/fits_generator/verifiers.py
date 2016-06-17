# Copyright (C) 2010-2011 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

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
Contains a number of helper functions useful for verifying header
card values.
"""

# STDLIB
import re

def is_empty(val):
    """
    Verifies that the value is empty.
    """
    return val == ''

def is_integer(val):
    """
    Verifies that the value is an integer.
    """
    return isinstance(val, int)
is_int = is_integer

def is_float(val):
    """
    Verifies that the value is a float.
    """
    return isinstance(val, float)

def is_numeric(val):
    """
    Verifies that the value is numeric (either an integer or a
    float).
    """
    return isinstance(val, (int, float))

def is_string(val):
    """
    Verifies that the value is a string.
    """
    return isinstance(val, str)

def is_bool(val):
    """
    Verifies that the value is a boolean ('T' or 'F').
    """
    return isinstance(val, bool)

def is_date(val):
    """
    Verifies that the value is a date or datetime of the form::

        YYYY-MM-DD

    or::

        YYYY-MM-DDThh:mm:ss.ssssss
    """
    return re.match(
        r"^[0-9]{4}-[0-9]{2}-[0-9]{2}(T[0-9]{2}:[0-9]{2}:[0-9]{2}(\.[0-9]+)?)?$",
        val) is not None

def is_cds(val):
    """
    Verifies that the value is a CDS number.
    """
    if not isinstance(val, int):
        raise ValueError(
            "Expected integer value for CDS")

    day = (value & 0xffff000000000000) >> (32 + 16)
    millis = (value & 0xffffffff0000) >> 16
    submillis = (value & 0xffff)

    success = True
    if millis > 86400999:
        error_collector(
            "CDS milliseconds field out of range (%d)" % millis,
            state)
        success = False

    if submillis > 999:
        error_collector(
            "CDS submilliseconds field out of range (%d)" % millis,
            state)
        success = False

    return success

def in_range(val, min=0, max=0, inclusive=True):
    """
    Verifies that the value is in the given range [min, max].  If
    *inclusive* is False, min and max will not be considered part
    of the valid range.
    """
    assert isinstance(val, (int, float))
    assert isinstance(min, (int, float))
    assert isinstance(max, (int, float))

    if ((inclusive and val < min or val > max) or
        (val <= min or val >= max)):
        raise ValueError(
            "Value outside of range %s - %s" % (min, max))
    return True

def regex(val, pattern):
    """
    Verifies that the value matches the given regular expression.
    This uses Python's regular expression syntax, except a complete
    match is always implied by surrounding the pattern in '^...$'
    """
    return re.match(val, '^%s$' % pattern)

T = True
F = False
len = len

__all__ = [
    'is_empty',
    'is_integer',
    'is_int',
    'is_float',
    'is_numeric',
    'is_string',
    'is_bool',
    'is_date',
    'is_cds',
    'in_range',
    'regex',
    'T',
    'F',
    'len'
    ]
