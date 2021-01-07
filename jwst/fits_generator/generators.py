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
Contains a number of helper functions for generating header card
values.
"""

# STDLIB
import datetime

from astropy.time import Time

# LOCAL
from . import util
from jwst import __version__

bool = __builtins__['bool']
str = __builtins__['str']
int = __builtins__['int']
float = __builtins__['float']

def date_and_time_to_cds(date_str, time_str):
    """
    Given a date and a time from different source keywords,
    *date_key* and *time_key*, generates a CDS number.
    """
    date = datetime.datetime.strptime(
        date_str, "%Y-%m-%d")
    time = datetime.datetime.strptime(
        time_str[:time_str.find('.')], "%H:%M:%S")

    dt = datetime.datetime(
        date.year, date.month, date.day,
        time.hour, time.minute, time.second, time.microsecond)

    return util.datetime2cds(dt)

def map(value, mapping={}):
    """
    Given a source value with the keyword *key*, returns the
    result of mapping it through the given dictionary-like
    *mapping*.  If *key* is None, the source key is the same as
    the destination key.
    """
    if value not in mapping:
        raise ValueError(
            "Unexpected source value '%s'" % value)

    return mapping[value]

def now(date_only=False):
    """
    Generates a date or datetime value from the current UTC time.
    """
    dt = datetime.datetime.now()
    if date_only:
        return dt.strftime("%Y-%m-%d")
    else:
        return dt.strftime("%Y-%m-%dT%H:%M:%S")

def version():
    """
    Generates the version of fits_generator software
    """
    return __version__

def substr(value, start=0, end=None):
    """
    Copies a substring of the value with the given *key* from the
    source FITS file.  If *key* is None, the source key is the
    same as the destination key.
    """
    if start >= len(value) or end >= len(value):
        raise ValueError(
            "Value shorter than expected.  Wanted [%d:%d], got %d" %
            (start, end, len(value)))

    return value[start:end]

def fromfitstime(fitstimestring):
    return datetime.datetime.strptime(fitstimestring, "%Y-%m-%dT%H:%M:%S.%f")

def tofitstime(datetimeobject):
    return datetimeobject.strftime("%Y-%m-%dT%H:%M:%S.%f")

def toMJD(timestring):
    #
    # Convert date-time to MJD for EXPSTART, EXPMID and EXPEND keywords
    oldtime = Time(timestring, format='isot', scale='utc')
    return oldtime.mjd

T = True
F = False

__all__ = [
    'substr',
    'bool',
    'date_and_time_to_cds',
    'map',
    'now',
    'version',
    'fromfitstime',
    'tofitstime',
    'toMJD',
    'T',
    'F',
    ]
