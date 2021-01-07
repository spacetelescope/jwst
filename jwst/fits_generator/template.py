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
from operator import itemgetter
import os
import re

# THIRD-PARTY
import numpy as np

# LOCAL
from . import objects
from . import generators
from . import util
from . import verifiers
from pyparsing import *

class TemplateParserBase:
    # This parser is a hybrid between hand-written parsing and
    # pyparsing LBNF.  Since pyparsing does not handle #include files
    # or line-oriented languages very well, the line-level parsing is
    # hand-written, and each type of line uses a separate
    # pyparsing-based grammar.

    def __init__(self):
        self._create_parser()

    def _raise(self, message):
        return SyntaxError(
            '%s:%d: %s' % (self._filename, self._lineno, message))

    def _create_parser(self):
        # Newlines are meaningful in this grammar, so don't treat them
        # as ordinary whitespace
        ParserElement.setDefaultWhitespaceChars(' \t')

        identifier = Regex(r"[A-Za-z_][A-Za-z0-9_]*")
        keyword = Regex(r"[A-Z0-9_\-]{1,8}")
        comment = (
            Suppress('/') +
            Regex(".*", re.DOTALL)
            )
        section_start = Suppress('<<')
        section_end = Suppress('>>')

        self._file_section_line = (
            section_start +
            Keyword('file') +
            Optional(identifier, None) +
            section_end
            )
        self._inherit_section_line = (
            section_start +
            Keyword('inherit') +
            quotedString.setParseAction(removeQuotes) +
            section_end
            )
        self._header_section_line = (
            section_start +
            Keyword('header') +
            Optional(identifier, None) +
            section_end
            )
        self._data_section_line = (
            section_start +
            Keyword('data') +
            section_end
            )

        self._definition_line = (
            Group(
                Optional(
                    keyword +
                    Optional('?', '+') +
                    Suppress('=') +
                    # Read until the comment "/", except ignore "/"
                    # inside of quotes.
                    Regex("([^/\"']*(\"[^\"]*\")?('[^']*')?)*", re.DOTALL)
                    )
                ) +
            Optional(
                comment, ''
                )
            )
        self._function_line = (
            Group(
                Optional(
                    Regex(".*", re.DOTALL)
                    )
                ) +
            Optional(
                comment, ''
                )
            )

        self._include_line = (
            Suppress(
                '#' +
                Keyword('include')
                ) +
            quotedString.setParseAction(removeQuotes)
            )

    def _line_iter(self, filename):
        """
        Iterates through the lines of a file, transparently handling
        continuation lines and #include directives.
        """
        old_filename = self._filename
        self._filename = filename

        returned_line = []
        with open(filename, 'rt') as f:
            i = 0
            for line in f:
                line = line.rstrip()
                i += 1

                # See if the line looks like an include directive
                try:
                    include_filename, = self._include_line.parseString(
                        line, True)
                except ParseException:
                    # Handle ignored comments
                    if line.startswith('#'):
                        continue

                    # It wasn't an include line, so yield it
                    self._lineno = i
                    if len(line) and line[-1] == '\\':
                        # Continuation line
                        returned_line.append(line[:-1])
                    else:
                        returned_line.append(line)
                        yield '\n'.join(returned_line)
                        returned_line = []
                else:
                    # Recurse into another included file and return
                    # its lines
                    for x in self._line_iter(
                        os.path.join(os.path.dirname(filename),
                                     include_filename)):
                        yield x

        self._filename = old_filename

    def parse(self, filename):
        self._filename = None

        lines = self._line_iter(filename)

        # <<file>> line
        line = next(lines)
        try:
            _, name = self._file_section_line.parseString(line, True)
        except ParseException:
            raise self._raise("file does not begin with <<file>> line")
        self._parse_file(name)

        line = next(lines)
        if isinstance(self, VerificationTemplateParser):
            while True:
                try:
                    _, inherit = self._inherit_section_line.parseString(line, True)
                except ParseException:
                    break
                else:
                    inherit = os.path.join(os.path.dirname(filename), inherit)
                    self._parse_inherit(inherit)
                    line = next(lines)

        stop = False
        while True:
            # <<header>> line
            try:
                _, name = self._header_section_line.parseString(line, True)
            except ParseException:
                raise self._raise("Expected <<header>> line")
            self._parse_header(name)

            # keyword definition lines
            while True:
                try:
                    line = next(lines)
                except StopIteration:
                    stop = True
                    break
                if line.startswith('<<'):
                    break
                try:
                    definition, comment = self._definition_line.parseString(
                        line, True)
                except ParseException:
                    raise self._raise("Expected keyword definition")
                else:
                    self._parse_card(definition, comment)

            if stop:
                break

            if isinstance(self, DataTemplateParser):
                continue

            # <<data>> line
            try:
                self._data_section_line.parseString(line, True)
            except ParseException:
                # Data sections are optional, so loop around to
                # see if it's a header
                continue

            self._parse_data()
            # data function lines
            while True:
                try:
                    line = next(lines)
                except StopIteration:
                    stop = True
                    break
                if line.startswith('<<'):
                    break
                try:
                    function, comment = self._function_line.parseString(
                        line, True)
                except ParseException:
                    raise self._raise("Expected function call")
                function = function[0].strip()
                if function != '':
                    self._parse_data_function(function)
            if stop:
                break

        return self._return_parse_value()

    def _parse_file(self, name):
        self._file_obj = objects.File(name=name)

    def _parse_header(self, name):
        self._header = objects.Header()
        self._hdu = objects.HDU(name=name, header=self._header)
        self._file_obj.append(self._hdu)

    def _parse_card(self, definition, comment):
        self._add_card(self._header, definition, comment)

    def _parse_data(self):
        self._data = objects.Data()
        self._data.generate_function = None
        self._hdu.data = self._data

    def _parse_data_function(self, function):
        self._add_data_function(self._data, function)

    def _return_parse_value(self):
        return self._file_obj


class GeneratorTemplateParser(TemplateParserBase):
    class _KeywordGeneratorNamespace:
        def __init__(self, card, fitsfiles, state):
            self._card = card
            self._fitsfiles = fitsfiles
            self._state = state

        def __getitem__(self, item):
            if item in self._fitsfiles:
                return self._fitsfiles[item].get_getter(
                    self._card.name, self._state.hdu)
            elif hasattr(generators, item):
                return getattr(generators, item)
            else:
                raise NameError("Unknown identifier '%s'" % item)

    def _make_generator(self, s):
        code = compile('(%s)' % s, '%s:%d' % (self._filename, self._lineno), 'eval')

        def run(self, fitsfiles, error_collector, state):
            return eval(
                code, {},
                GeneratorTemplateParser._KeywordGeneratorNamespace(
                    self, fitsfiles, state))
        return run

    def _add_card(self, header, definition, comment):
        if len(definition):
            optional = definition[1] == '?'
            expr = definition[2]
            try:
                generator = self._make_generator(expr)
            except SyntaxError as e:
                raise self._raise(
                    "Invalid syntax: \n%s\n%s^\n" %
                    (expr, ' ' * (e.offset - 1)))
            key = definition[0]
            card = objects.Card(key, comment=comment, generate=generator,
                                optional=optional)
        else:
            card = comment
        header.append(card)

    class _DataGeneratorNamespace:
        def __init__(self, fitsfiles, state):
            self._fitsfiles = fitsfiles
            self._state = state

        def __getitem__(self, item):
            if item in self._fitsfiles:
                return self._fitsfiles[item].get_getter(None, None)
            elif item == 'np':
                return np
            else:
                raise NameError("Unknown identifier '%s'" % item)

    def _add_data_function(self, data, function):
        code = compile('(%s)' % function, '%s:%d' % (self._filename, self._lineno), 'eval')

        def run(fitsfiles, error_collector, state):
            return eval(
                code, {},
                GeneratorTemplateParser._DataGeneratorNamespace(
                    fitsfiles, state))
        data.generate_function = run


class VerificationTemplateParser(TemplateParserBase):
    class _VerifierNamespace:
        def __init__(self, val, header):
            self._val = val
            self._header = header

        def __getitem__(self, item):
            if item == 'x':
                return self._val
            elif item == 'output':
                def getter(kwd):
                    return self._header[kwd]
                return getter
            elif hasattr(verifiers, item):
                return getattr(verifiers, item)
            else:
                raise NameError("Unknown identifier '%s'" % item)

    def _make_verifier(self, s):
        code = compile('(%s)' % s, '%s:%d' % (self._filename, self._lineno), 'eval')

        def run(val, error_collector, state):
            try:
                result = eval(
                    code, {},
                    VerificationTemplateParser._VerifierNamespace(val, state.header))
            except ValueError as e:
                error_collector("'%s': %s in '%s'" % (val, str(e), s), state)
                return False
            if result is False:
                error_collector("'%s' did not match '%s'" % (val, s), state)
            return result
        return run

    def _add_card(self, header, definition, comment):
        if len(definition):
            optional = (definition[1] == '?')
            expr = definition[2]
            try:
                verifier = self._make_verifier(expr)
            except SyntaxError as e:
                raise self._raise(
                    "Invalid syntax: \n%s\n%s^\n" %
                    (expr, ' ' * (e.offset - 1)))
            key = definition[0]
            card = objects.Card(key, comment=comment, verify=verifier,
                                optional=optional)
        else:
            card = comment
        header.append(card)

    def _add_data_function(self, data, function):
        pass

    def _parse_file(self, name):
        TemplateParserBase._parse_file(self, name)
        self._file_obj.inherit = []

    def _parse_inherit(self, inherit):
        print(inherit)
        self._file_obj.inherit.append(self.__class__().parse(inherit))


class DataTemplateParser(TemplateParserBase):
    def _add_card(self, header, definition, comment):
        if len(definition):
            try:
                value = eval(definition[2])
            except Exception as e:
                raise self._raise(
                    "Can't parse as literal '%s'" % definition[2])
            key = definition[0]
            card = objects.Card(key, value=value)
        else:
            card = comment
        header.append(card)


def merge_cards(gen, val, debug=False):
    # Deal with inherited cards first, that way all further- nested
    # cards will override them
    for inherit in val.inherit:
        merge_cards(gen, inherit)

    for i, gen_hdu in enumerate(gen):
        # If there are more generator HDUs than validator HDUs, just
        # reuse the last validator HDU indefinitely.
        if i >= len(val):
            val_hdu = val[-1]
        else:
            val_hdu = val[i]
        for card in gen_hdu.header:
            if isinstance(card, objects.Card):
                if card.name in val_hdu.header.map:
                    card.verify_function = val_hdu.header.map[card.name].verify_function
    return gen

def get_generator(filetype, debug=False):
    """
    Loads a file template and returns an instance of its top-level
    `File` class.
    """
    if os.path.exists(filetype):
        gen_filename = filetype
    else:
        gen_filename = os.path.join(util.get_templates_dir(), filetype + ".gen.txt")
        if not os.path.exists(gen_filename):
            raise ValueError("No generator template available for '%s'" % filetype)
    gen = GeneratorTemplateParser().parse(gen_filename)

    val_filename = gen_filename[:-8] + ".val.txt"
    if os.path.exists(val_filename):
        val = VerificationTemplateParser().parse(val_filename)
        return merge_cards(gen, val, debug)
    else:
        return gen
