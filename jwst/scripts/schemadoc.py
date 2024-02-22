#!/usr/bin/env python

# Copyright (C) 2018 Association of Universities for Research in Astronomy (AURA)

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

import argparse
import sys

from stdatamodels.schema import build_docstring


def get_docstrings(template, model_names, all=False):
    # Get the docstring for every model class
    from stdatamodels.jwst.datamodels import _defined_models as defined_models

    if all:
        klasses = defined_models
    else:
        klasses = {}
        for model_name in model_names:
            klasses[model_name] = defined_models[model_name]

    for klassname, klass in klasses.items():
        try:
            docstring = build_docstring(klass, template)
        except Exception as err:
            print(klassname, ':', str(err), file=sys.stderr)
        else:
            print('.. ' + klassname + ' ..')
            print(docstring, end='')


def main():
    long_description = """
    Create documentation from the schema file of a datamodel class
    """

    parser = argparse.ArgumentParser(description=long_description)
    parser.add_argument('-a', '--all', action='store_true',
                        help='generate docstring for all models')
    parser.add_argument('-t', '--template', type=str,
                        help='input file containing templates')
    parser.add_argument('models', nargs='*',
                        help='Models to generate docstrings for')
    args = parser.parse_args()

    if args.template is None:
        template = '{path} : {title} ({datatype})\n'
    else:
        with open(args.template, 'r') as fd:
            template = fd.readlines()
            template = ''.join(template)

    get_docstrings(template, args.models, all=args.all)


if __name__ == '__main__':
    main()
