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
This file contains classes for loading the different input file types
we handle, currently:

  - FITS
  - Data (a simple version of FITS headers, see template.py)
  - APT proposal files

Each type is defined in a class that allows the template parser to
grab values in appropriate ways from the input file.
"""

# STDLIB
import os
from xml.etree import ElementTree
import zipfile

# THIRD-PARTY
from astropy.io import fits as pyfits

# LOCAL
from ..fits_generator import objects
from ..fits_generator import template

class FITSGetter:
    """
    A getter used by both FITS and data files.
    """
    def __init__(self, mapping, key=None, hdu=None):
        self._mapping = mapping
        self._key = key
        self._hdu = hdu

    def __getitem__(self, item):
        return self._mapping[item]

    def __call__(self, key=None, hdu=None, default=''):
        """
        Get the keyword 'key' in the given hdu number 'hdu'.

        If *key* is None, the name of the output key will be used.

        If *hdu* is None, the same-numbered hdu as the output file
        will be used.
        """
        if key is None:
            use_key = self._key
        else:
            use_key = key
        if hdu is None:
            use_hdu = self._hdu
        else:
            use_hdu = hdu
        return self._mapping[use_hdu].header.get(use_key, default)

class InputFileType():
    """
    Base class for an input file type.
    """
    type_name = "unknown"

    def __init__(self, obj):
        """
        Given a file path or object representing that file, create a
        new file type object.
        """
        raise NotImplementedError()

    @staticmethod
    def is_file(file_id):
        """
        *file_id* is the first 6 characters of the file.  Returns
        `True` if the magic token indicates that the file is probably
        of this type.
        """
        raise NotImplementedError()

    @staticmethod
    def is_object(obj):
        """
        Returns `True` if *obj* is a Python object representing the
        contents of a file of this type.
        """
        raise NotImplementedError()

    def get_getter(self, name, hdu):
        """
        Returns a callable object used to get values from the file.
        This callable object is present in the generator template
        namespace and used to get values from the input file.
        """
        raise NotImplementedError()

    def get_name(self):
        """
        Get the logical name of the input file.  This name is used as
        the function name present in the generator template namespace.
        """
        raise NotImplementedError()

    def close(self):
        """
        Close the input file if necessary.
        """
        pass

class InputFITSFile(InputFileType):
    type_name = "FITS"

    def __init__(self, obj):
        if isinstance(obj, pyfits.HDUList):
            self._hdulist = obj
            self._close = False
        else:
            self._hdulist = pyfits.open(obj)
            self._close = True

    @staticmethod
    def is_file(file_id):
        return file_id.startswith(b'SIMPLE')

    @staticmethod
    def is_object(obj):
        return isinstance(obj, pyfits.HDUList)

    def get_getter(self, _name, _hdu):
        return FITSGetter(self._hdulist, _name, _hdu)

    def get_name(self):
        return 'input'

    def close(self):
        if self._close:
            self._hdulist.close()

class InputDataFile(InputFileType):
    type_name = "data"

    def __init__(self, obj):
        if isinstance(obj, objects.File):
            self._mapping = obj
        else:
            parser = template.DataTemplateParser()
            self._mapping = parser.parse(obj)

    @staticmethod
    def is_file(file_id):
        return file_id.startswith(b'<<file')

    @staticmethod
    def is_object(obj):
        return isinstance(obj, objects.File)

    def get_getter(self, _name, _hdu):
        return FITSGetter(self._mapping, _name, _hdu)

    def get_name(self):
        return self._mapping.name

class InputAPTFile(InputFileType):
    type_name = "APT"

    def __init__(self, obj):
        basename, ext = os.path.splitext(os.path.basename(obj))
        with zipfile.ZipFile(obj, 'r') as zip:
            with zip.open(basename + ".xml", 'r') as xml:
                self._tree = ElementTree.parse(xml)
        # Remove all namespaces for convenience
        for element in self._tree.getroot().iter():
            element.tag = element.tag[element.tag.find('}') + 1:]

    @staticmethod
    def is_file(file_id):
        # The 'magic' bytes indicating a zip file
        return file_id.startswith(b'\x50\x4b\x03\x04')

    @staticmethod
    def is_object(obj):
        return False

    def get_getter(self, _name, _hdu):
        def getter(xpath, default=None):
            result = self._tree.findtext(xpath)
            result = result.strip().replace('\n', '')
            if result:
                return result

            # Try element name
            result = self._tree.find(xpath)
            if result is not None and len(result):
                return result[0].tag

            if default is not None:
                return default
            return ''
        return getter

    def get_name(self):
        return 'apt'

def get_inputfiles(input_files):
    """
    Given a list of file paths, creates a dictionary mapping the
    file's logical name (the function name used in the generator
    template namespace) to a InputFileType object that can be used to
    get values from that file.
    """
    data_files = {}

    for f in input_files:
        obj = None
        if isinstance(f, str):
            with open(f, 'rb') as fd:
                file_id = fd.read(6)
            for file_type in input_file_types:
                if file_type.is_file(file_id):
                    obj = file_type(f)
                    break
        else:
            for file_type in input_file_types:
                if file_type.is_object(f):
                    obj = file_type(f)
                    break

        if obj is not None:
            if obj.get_name() in data_files:
                raise ValueError(
                    "Specified more than one '%s' file" % obj.get_name())
            data_files[obj.get_name()] = obj
        else:
            raise ValueError(
                "'%s' was not recognized as a %s file" %
                (f, file_type_names))

    if 'input' not in data_files:
        raise ValueError(
            "Must specify at least one input FITS file.")

    return data_files

def is_fitswriter(hdulist):
    """Returns True if the file was generated by the FITSWriter
    FITSWriter files have a HISTORY record that includes the string
    'FITSWriter'.  Note that MIRI DHAS LVL2 files also have the same
    HISTORY record...
    """
    target_string = 'FITSWriter'
    try:
        history = hdulist[0].header['HISTORY']
        for card in history:
            if card.find(target_string) != -1:
                return True
        return False
    except KeyError:
        return False

def is_ncont(hdulist):
    """Returns True if the file was generated by the NIRCAM NCont Post
    processor.  These files have a HISTORY record that includes the string
    'NCont Post processor'.
    """
    target_string = 'NCont Post processor'
    try:
        history = hdulist[0].header['HISTORY']
        for card in history:
            if card.find(target_string) != -1:
                return True
        return False
    except KeyError:
        return False

def is_swts(hdulist):
    """Returns True if the file was generated by the NIRCAM task swts2fits.
    These files have a HISTORY record that includes the string
    'swts2fits'.
    """
    target_string = 'swts2fits'
    try:
        history = hdulist[0].header['HISTORY']
        for card in history:
            if card.find(target_string) != -1:
                return True
    finally:
        return False

def is_miri(hdulist):
    """Returns True if the file is MIRI data.  These files have a
    INSTRUME keyword set to 'MIRI'."""

    keywordname = 'INSTRUME'
    target_value = 'MIRI'
    if hdulist[0].header[keywordname].strip() == target_value:
        return True
    else:
        return False

def is_miri_ifu(hdulist):
    """Returns True if the EXP_TYPE keyword is set to MIR_MRS"""
    keywordname = 'EXP_TYPE'
    target_value = 'MIR_MRS'
    try:
        if hdulist[0].header[keywordname].strip() == target_value:
            return True
        else:
            return False
    except:
        return False

def is_miri_lrs(hdulist):
    """Returns True if the filter (FWA_POS) is P750L"""
    keywordname = 'FWA_POS'
    target_value = 'P750L'
    try:
        if hdulist[0].header[keywordname].strip() == target_value:
            return True
        else:
            return False
    except:
        return False

def is_nirspec_fm1(hdulist):
    """Returns True if the file is NIRSpec FM1 data.  These files
    have a DATE-OBS in the first half of 2011"""
    instrument = hdulist[0].header['INSTRUME']
    if instrument != 'NIRSPEC': return False
    date = hdulist[0].header['DATE-OBS']
    year = date[:4]
    target_year = '2011'
    month = int(date[5:7])
    target_month_range = (0, 7)
    if year == target_year and month > target_month_range[0] and month < target_month_range[1]:
        return True
    return False

def is_nircam_fm1(hdulist):
    """Returns True is the file is NIRCAM FM1 data.  These files
    have a DATE-OBS after the first half of 2012"""
    instrument = hdulist[0].header['INSTRUME']
    if instrument != 'NIRCAM': return False
    date = hdulist[0].header['DATE-OBS']
    year = int(date[:4])
    target_year = 2012.5
    month = int(date[5:7])
    year = year + month / 12.0
    if year >= target_year:
        return True
    return False

def is_nirspec_ips(hdulist):
    """Returns True if the file is NIRSpec IPS simulation data.  These
    file have keywords that begin with 'IPS' - I choose IPSDET as a sentinel
    for these type of data"""
    keywordname = 'IPSDET'
    try:
        if hdulist[0].header[keywordname].strip():
            return True
    except:
        return False

def is_nirspec_irs2(hdulist):
    """Returns true if the file uses NIRSPEC IRS2 readout pattern.  The
    READOUT keyword will contain the string "IRS2"
    """
    keywordname = 'READOUT'
    try:
        keywordvalue = hdulist[0].header[keywordname]
    except KeyError:
        print("No keyword 'READOUT'")
        print("Returning False for is_nirspec_irs2")
        return False
    if keywordvalue.upper().find('IRS2') != -1:
        return True
    else:
        return False

def is_tfi(hdulist):
    """Returns True if the INSTRUME string is 'TFI'."""
    try:
        if hdulist[0].header['INSTRUME'].strip() == 'TFI':
            return True
        else:
            return False
    except:
        return False

def is_niriss(hdulist):
    """Returns True if the INSTRUME string is 'NIRISS'."""
    try:
        if hdulist[0].header['INSTRUME'].strip() == 'NIRISS':
            return True
        else:
            return False
    except:
        return False

def is_niriss_spec_vert(hdulist):
    """Returns True if the NIRISS spectrum is horizontal i.e. for
    PUPIL = GR700XD or FILTER = GR150R. """
    try:
        if hdulist[0].header['PWCCRPUP'].strip() == 'GR700XD' or \
                hdulist[0].header['FWCCRFIL'].strip() == 'GR150R':
            return True
        else:
            return False
    except:
        return False

def is_niriss_spec_horiz(hdulist):
    """Returns True if the NIRISS spectrum is vertical, i.e. for
    FILTER = GR150C. """
    try:
        if hdulist[0].header['FWCCRFIL'].strip() == 'GR150C':
            return True
        else:
            return False
    except:
        return False

input_file_types = [
    InputFITSFile,
    InputDataFile,
    InputAPTFile
    ]

file_type_names = [
    x.type_name for x in input_file_types]

file_type_names[-1] = 'or ' + file_type_names[-1]
file_type_names = ', '.join(file_type_names)
