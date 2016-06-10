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
from cStringIO import StringIO
import glob
import os
import re
import textwrap

# LOCAL
from jwst_tools.fits_generator import template
from jwst_tools.fits_generator import util

def write_if_different(filename, data):
    data = data.encode('ascii')

    if os.path.exists(filename):
        fd = open(filename, 'rb')
        original_data = fd.read()
        fd.close()
    else:
        original_data = None

    if original_data != data:
        fd = open(filename, 'wb')
        fd.write(data)
        fd.close()

class TemplateToRst(template.TemplateParserBase):
    def _print(self, string=''):
        self._output.write(string)
        self._output.write('\n')

    def _parse_file(self, name):
        self._output = StringIO()
        self._print(name)
        self._print('=' * len(name))
        self._print()

    def _parse_header(self, name):
        self._print('%s' % name)
        self._print('-' * len(name))
        self._print()
        self._print('Header')
        self._print('``````')
        self._print()
        self._print('.. csv-table::')
        self._print('   :widths: 7, 2, 2, 42, 42')
        self._print('   :quote: |')
        self._print()
        self._print()

    def _parse_card(self, definition, comment):
        if len(definition):
            key, optional, expr = definition
            if optional == '?':
                optional = '*?*'
            else:
                optional = ''
            if comment:
                comment = '/ ' + comment
            expr = expr.strip()
            expr = re.sub('\n\\s*', ' ', expr)
            self._print('   |``%s``|, |%s|, |=|, |``%s``|, |%s|' %
                        (key, optional, expr, comment))
        elif comment:
            self._print('   ||, ||, ||, |**%s**|, ||' % comment)
        else:
            self._print()

    def _parse_data(self):
        self._print('Data')
        self._print('````')
        self._print()

    def _parse_data_function(self, function):
        self._print('.. code-block:: python')
        self._print()
        self._print('  (')
        self._print('  ' + function.replace('\n', '\n  '))
        self._print('  )')
        self._print()

    def _return_parse_value(self):
        return self._output.getvalue()

def convert_templates(outdir):
    templates_dir = util.get_templates_dir()

    files = set()
    parser = TemplateToRst()
    for generator in glob.glob(os.path.join(templates_dir, "*.gen.txt")):
        try:
            data = parser.parse(generator)
        except IOError:
            pass
        else:
            outfile = os.path.basename(generator) + '.rst'
            write_if_different(os.path.join(outdir, outfile), data)
            files.add(outfile)

    files = list(files)
    files.sort()
    with open(os.path.join(outdir, 'templates.rst'), 'wb') as fd:
        fd.write('Templates\n')
        fd.write('=========\n')
        fd.write('\n')
        fd.write('.. toctree::\n\n')
        for outfile in files:
            fd.write('   %s\n' % outfile)

if __name__ == '__main__':
    convert_templates('source')
