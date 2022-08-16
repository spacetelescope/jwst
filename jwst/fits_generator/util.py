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

# STDLIB
import collections
import datetime
import os
import sys
import textwrap

# THIRD-PARTY
from astropy.io import fits as pyfits
from astropy.time import Time

def get_templates_dir():
    """
    Returns a path to the directory containing definitions of the file
    types.
    """
    return os.path.join(os.path.dirname(__file__), 'templates')

def get_fitsfile(fitsfile):
    """
    Given a filename or pyfits.HDUList object, returns a (filename,
    pyfits.HDUList, opened) triple.
    """
    from . import input_file_types
    if isinstance(fitsfile, str):
        return fitsfile, pyfits.open(fitsfile), True
    elif isinstance(fitsfile, pyfits.HDUList):
        return fitsfile.filename, fitsfile, False
    elif isinstance(fitsfile, input_file_types.InputFITSFile):
        return fitsfile._hdulist.filename, fitsfile._hdulist, False

def get_error_collector(error_collector=None):
    """
    Get a usable error collecting function.  If `error_collector` is
    `None`, returns a default function that writes to `stderr`.
    """
    def error_collector_stderr(s, state):
        content = "%s: %s" % (str(state), s)
        wrapped = textwrap.fill(content, subsequent_indent='    ')
        sys.stderr.write(wrapped)
        sys.stderr.write('\n')

    if error_collector is None:
        return error_collector_stderr
    # TODO: Assert it's callable
    return error_collector

def datetime2cds(dt):
    # SPEC: What is the epoch for a CDS value?
    delta = dt - datetime.datetime(1970, 1, 1)
    cds = ((int(delta.days) & 0xffff) << (16 + 32) |
           (int(delta.seconds * 1000) + int(delta.microseconds / 1000.) & 0xffffffff) << 16 |
           (int(delta.microseconds % 1000) & 0xffff))
    return cds

def toMJD(timestring):
    #
    # Convert date-time to MJD for EXPSTART, EXPMID and EXPEND keywords
    oldtime = Time(timestring, format='isot', scale='utc')
    return oldtime.mjd

def issequence(obj):
    """
    Returns True if object quacks like a sequence.
    """
    return getattr(obj, '__iter__', False)

if sys.hexversion >= 0x02060000:
    def iscallable(obj):
        """
        Returns True if object quacks like a callable.
        """
        return isinstance(obj, collections.abc.Callable)
else:
    def iscallable(obj):
        """
        Returns True if object quacks like a callable.
        """
        return hasattr(obj, '__call__')
