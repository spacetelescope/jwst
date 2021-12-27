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
import datetime
import os
import re
import warnings

# THIRD-PARTY
from astropy.io import fits as pyfits

# LOCAL
from ..fits_generator import util

# Set DEBUG_TYPES to True to get warning information about problems in
# the type definitions themselves.
DEBUG_TYPES = True

class ParseState():
    """
    Keeps track of the current state of the parser -- that is which
    HDU and card is currently being examined.  Used to give detailed
    information in error log messages.
    """
    def __init__(self):
        self.file = ''
        self.hdu = 0
        self.card = None
        self.depth = 0

    def __repr__(self):
        val = "%s:%02d" % (self.file, self.hdu)
        if self.card is not None:
            val += ":" + self.card
        return val


class Object():
    """
    The base class of all objects that form the FITS file hierarchy.
    """
    def verify(self, value, error_collector, state):
        """
        Verifies *value* against the type.

        Parameters
        ----------
        value : object
            The value to verify.

        error_collector : callable
            A callable that accepts (*message*, *state*) that will be
            called when a verification fails.  *message* is a string
            containing a reason for the failure.  *state* is a
            `ParseState` object containing information about the
            location in the FITS file where the error occurred.

        state : `ParseState` instance
            Keeps track of the location in the source file where the
            value was found.

        Returns
        -------
        success : boolean
            Returns `True` if the value is valid, otherwise `False`.
            Extra information about the error will be logged to the
            *error_collector* function.
        """
        raise NotImplementedError(
            "This method must be overridden in the subclass.")

    def generate(self, input_files, error_collector, state):
        """
        Generates a new value of the type.

        Parameters
        ----------
        input_files : dict
            A dictionary mapping input file types to `pyfits.HDUList`
            objects.

        error_collector : callable
            A callable that accepts (*message*, state*) that will be
            called when generation fails.  *message* is a string
            containing a reason for the failure.  *state* is a
            `ParseState` object containing information about the
            location in the output FITS file where the error occurred.

        state : `ParseState` instance
            Keeps track of the location in the output file where the
            value was being output.

        Returns
        -------
        result : object
            The generated object.
        """
        raise NotImplementedError(
            "This method must be overridden in the subclass.")

    def describe(self, stream, state):
        """
        Generates a document describing the given type.

        Parameters
        ----------
        stream : a writable file-like object
            The content of the document is written to this stream.

        state : `ParseState` instance
            Keeps track of the location in the output file where the
            value was being output.

        Returns
        -------
        None
        """
        raise NotImplementedError(
            "This method must be overridden in the subclass.")

    def _write_section(self, stream, content, state):
        sections = "=-~'"

        stream.write("%s\n" % content)
        stream.write(sections[state.depth] * len(content))
        stream.write("\n\n")


class Card(Object):
    """
    Defines how to verify and generate a key/value card.

    """
    def __init__(self, name, value='#UNKNOWN', comment='',
                 verify=None, generate=None, optional=True):
        """
        Parameters
        ----------

        name : str
            The key

        value : object
            A constant value, if one exists

        comment : str
            The comment to add to the card when writing it out.

        verify : callable or None
            Each entry in the sequence should be a function taking the
            arguments (*value*, *error_collector*, *state*).

            - *value* is the input value (as a string) for the card

            - *error_collector* is a callable that should be called
               when the value is invalid in some way.  See `error
               collector functions`_.

            - *state* is a `ParserState` object containing the current
               state of the parser.

            The `Verify` class contains a number of methods for
            conveniently generating such functions.

        generate : callable or None
            A function to generate the value for the card.  It
            should be a function taking the arguments (*self*, *hdulist*,
            *error_collector*, *state*) and returning a string value.

            - *self* is the Card object

            - *hdulist* is a `pyfits.HDUList` object containing the
               content of the input FITS file.

            - *error_collector* is a callable that should be called
               when generation failed in some way.  See `error
               collector functions`_.

            - *state* is a `ParserState` object containing the current
               state of the parser.

            The `Gen` class contains a number of methods for
            conveniently generating such functions.

        optional : bool
            When `True`, the card is considered optional, and the
            verification step will not complain if the card is
            missing.  When `False`, if the card is missing, an error
            will be emitted.
        """
        Object.__init__(self)

        if not isinstance(name, str):
            raise TypeError("key must be a string")
        if not name == name.upper():
            raise TypeError(
                "key must be a string with all uppercase letters, got '%s'" %
                name)
        self.name = name

        if not isinstance(value, (str, int, float)):
            raise TypeError("value must be a string, int or float")
        self.value = value

        if not isinstance(comment, str):
            raise TypeError("comment must be a string")
        self.comment = comment

        if verify is not None and not util.iscallable(verify):
            raise TypeError(
                "verify must be a sequence of callables or a callable")
        self.verify_function = verify

        if generate is not None and not util.iscallable(generate):
            raise TypeError("generate must be a callable")
        self.generate_function = generate

        self.optional = bool(optional)

    def __repr__(self):
        return "<Card '%s'>" % self.name

    def match_key(self, key):
        return re.match(self.name, key)

    def verify(self, value, error_collector, state):
        success = True
        if self.verify_function:
            try:
                if not self.verify_function(value, error_collector, state):
                    success = False
            except Exception as e:
                error_collector(str(e), state)
                success = False
        return success
    verify.__doc__ = Object.verify.__doc__

    def generate(self, header, fitsfiles, error_collector, state):
        key = self.name
        if self.generate_function is None:
            value = '#NONE'
        else:
            try:
                value = self.generate_function(
                    self, fitsfiles, error_collector, state)
            except Exception as e:
                error_collector(str(e), state)
                value = '#ERROR'

        try:
            comment = self.get_comment(fitsfiles, error_collector, state)
        except Exception as e:
            error_collector(str(e), state)
            comment = '#ERROR'

        header.insert(len(header), (key, value, comment), after=True)
    generate.__doc__ = Object.generate.__doc__

    def get_comment(self, fitsfiles, error_collector, state):
        """
        Returns a comment for the card, given the input FITS file.
        """
        return self.comment

    def describe(self, stream, state):
        stream.write("``%s`` %s\n\n" % (self.name, self.comment))
    describe.__doc__ = Object.describe.__doc__


class Header(Object):
    """
    Defines a header, specifically the order and substance of the
    key/value cards it should contain.
    """
    def __init__(self, cards=None):
        """
        Parameters
        ----------
        cards : iterable of `Card` instances

            Defines how each card is verified and generated.
        """
        Object.__init__(self)

        if cards is None:
            self._cards = []
            self.map = {}
        else:
            self.map = {}
            for card in cards:
                if not isinstance(card, (Card, str)):
                    raise TypeError(
                        "Each element of cards sequence must be either a Card or a str")
                if isinstance(card, Card):
                    self.map[card.name] = card
            self._cards = cards[:]

    def __getitem__(self, key):
        if key in self.map:
            return self.map[key].value
        else:
            return self._cards[key]

    def get(self, key, default=None):
        if key in self.map:
            return self.map[key].value
        else:
            return default

    def __contains__(self, key):
        return key in self.map

    def append(self, item):
        if isinstance(item, str):
            self._cards.append(item)
        elif isinstance(item, Card):
            self._cards.append(item)
            self.map[item.name] = item
        else:
            raise TypeError(
                "Only strings and Card objects may be appended to a header")

    def __iter__(self):
        return iter(self._cards)

    def __len__(self):
        return len(self._cards)

    def verify(self, value, error_collector, state):
        success = True
        found_cards = {}
        state.header = value
        for key, val in value.items():
            if key in self.map:
                state.card = key
                if not self.map[key].verify(val, error_collector, state):
                    success = False
                state.card = None
                found_cards[key] = None
        for key, val in self.map.items():
            if not val.optional and key not in found_cards:
                error_collector(
                    "Missing required card '%s'" % key,
                    state)
        state.header = None
        return success
    verify.__doc__ = Object.verify.__doc__

    def create_header(self, fitsfiles, error_collector, state):
        return pyfits.Header()

    def generate(self, fitsfiles, error_collector, state):
        header = self.create_header(fitsfiles, error_collector, state)

        for card in self:
            if isinstance(card, Card):
                state.card = card.name
                card.generate(header, fitsfiles, error_collector, state)
                state.card = None
            elif isinstance(card, str):
                card_obj = pyfits.Card(' ', card)
                header.append(card_obj, useblanks=False)

        return header
    generate.__doc__ = Object.generate.__doc__

    def describe(self, stream, state):
        for card in self:
            if isinstance(card, Card):
                state.card = card.name
                card.describe(stream, state)
                state.card = None
            elif isinstance(card, str):
                stream.write('**%s**\n\n' % card)
    describe.__doc__ = Object.describe.__doc__


class Data(Object):
    """
    Describes how a data section is verified and generated.
    """

    def __init__(self):
        Object.__init__(self)
        self.generate_function = None

    def verify(self, value, error_collector, state):
        return True
    verify.__doc__ = Object.verify.__doc__

    def generate(self, fitsfiles, error_collector, state):
        if self.generate_function is not None:
            state.card = "<data>"
            try:
                return self.generate_function(fitsfiles, error_collector, state)
            except Exception as e:
                error_collector(str(e), state)
            state.card = None
    generate.__doc__ = Object.generate.__doc__

    def describe(self, stream, state):
        self._write_section(stream, 'Data', state)
    describe.__doc__ = Object.describe.__doc__


class HDU(Object):
    """
    Describes a single HDU (how a header and its data are bundled together.
    """
    def __init__(self, header=None, data=None, name=None):
        """
        Parameters
        ----------
        header : `Header` instance

        data : `Data` instance
        """
        Object.__init__(self)

        if not isinstance(header, (Header, type(None))):
            raise TypeError("arg 1 must be a Header instance")

        if not isinstance(data, (Data, type(None))):
            raise TypeError("arg 2 must be a Data instance")

        self.header = header
        self.data = data
        self.name = name

    def __repr__(self):
        return "<HDU %s>" % self.name

    def verify(self, value, error_collector, state):
        header = True
        if self.header is not None:
            assert isinstance(self.header, Header)
            header = self.header.verify(value.header, error_collector, state)

        data = True
        if self.data is not None:
            assert isinstance(self.data, Data)
            data = self.data.verify(value.data, error_collector, state)

        return header and data
    verify.__doc__ = Object.verify.__doc__

    def create_hdu(self, fitsfiles, error_collector, state):
        return pyfits.PrimaryHDU()

    def generate(self, fitsfiles, error_collector, state):
        hdu = self.create_hdu(fitsfiles, error_collector, state)
        fitsfiles['output'].append(hdu)
        # Append may have changed the class of the hdu so get it again
        # from pyfits
        hdu = fitsfiles['output'][-1]

        if self.header is not None:
            assert isinstance(self.header, Header)
            hdu.header = self.header.generate(fitsfiles, error_collector, state)

        if self.data is not None:
            assert isinstance(self.data, Data)
            hdu.data = self.data.generate(fitsfiles, error_collector, state)

        return hdu
    generate.__doc__ = Object.generate.__doc__

    def describe(self, stream, state):
        if state.hdu == 0:
            self._write_section(stream, 'Primary HDU', state)
        else:
            self._write_section(stream, 'Extension HDU %d' % state.hdu, state)

        state.depth += 1

        if self.header is not None:
            assert isinstance(self.header, Header)
            self.header.describe(stream, state)

        if self.data is not None:
            assert isinstance(self.data, Data)
            self.data.describe(stream, state)

        state.depth -= 1
    describe.__doc__ = Object.describe.__doc__


class File(list, Object):
    """
    Describes the order of the HDUs in a complete FITS file.
    """

    def __init__(self, hdus=[], name='FITS File'):
        """
        Parameters
        ----------
        hdus : list of `HDU` objects

            A list of HDUs as they should appear in the FITS file.

        name : string
            The name of the FITS file type
        """
        Object.__init__(self)

        for hdu in hdus:
            if not isinstance(hdu, HDU):
                raise TypeError("arg 1 must be list of HDU instances")

        list.__init__(self, hdus)
        self.name = name

    def verify(self, value, error_collector, state):
        success = True
        for i, (hdudef, hdu) in enumerate(zip(self, value)):
            assert isinstance(hdudef, HDU)
            state.hdu = i
            if not hdudef.verify(hdu, error_collector, state):
                success = False
        return success
    verify.__doc__ = Object.verify.__doc__

    def generate(self, fitsfiles, error_collector, state):
        output_hdulist = pyfits.HDUList()
        fitsfiles['output'] = output_hdulist

        for i, hdudef in enumerate(self):
            assert isinstance(hdudef, HDU)
            state.hdu = i
            hdudef.generate(fitsfiles, error_collector, state)

        state.hdu = None

        return output_hdulist
    generate.__doc__ = Object.generate.__doc__

    def describe(self, stream, state):
        self._write_section(stream, self.name, state)

        for i, hdudef in enumerate(self):
            assert isinstance(hdudef, HDU)
            state.hdu = i
            state.depth += 1
            hdu = hdudef.describe(stream, state)
            state.depth -= 1
    describe.__doc__ = Object.describe.__doc__
