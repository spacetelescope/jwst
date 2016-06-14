# -*- coding: utf-8 -*-
# Copyright (c) 2010, 2011, Sebastian Wiesner <lunaryorn@googlemail.com>
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


"""
    sphinxcontrib.programoutput
    ===========================

    This extension provides a directive to include the output of commands as
    literal block while building the docs.

    .. moduleauthor::  Sebastian Wiesner  <lunaryorn@googlemail.com>
"""

__version__ = '0.4.1'

import sys
import shlex
from subprocess import Popen, CalledProcessError, PIPE, STDOUT
from collections import defaultdict

from docutils import nodes
from docutils.parsers import rst
from docutils.parsers.rst.directives import flag, unchanged


class program_output(nodes.Element):
    pass


def _slice(value):
    parts = [int(v.strip()) for v in value.split(',')]
    if len(parts) > 2:
        raise ValueError('too many slice parts')
    return tuple((parts + [None]*2)[:2])


class ProgramOutputDirective(rst.Directive):
    has_content = False
    final_argument_whitespace = True
    required_arguments = 1

    option_spec = dict(shell=flag, prompt=flag, nostderr=flag,
                       ellipsis=_slice, extraargs=unchanged)

    def run(self):
        node = program_output()
        node['command'] = self.arguments[0]

        if self.name == 'command-output':
            node['show_prompt'] = True
        else:
            node['show_prompt'] = 'prompt' in self.options

        node['hide_standard_error'] = 'nostderr' in self.options
        node['extraargs'] = self.options.get('extraargs', '')
        node['use_shell'] = 'shell' in self.options
        if 'ellipsis' in self.options:
            node['strip_lines'] = self.options['ellipsis']
        return [node]


class ProgramOutputCache(defaultdict):
    """
    :class:`collections.defaultdict` sub-class, which caches program output.

    If a program's output is not contained in this cache, the program is
    executed, and its output is placed in the cache.
    """

    def __missing__(self, key):
        """
        Called, if a command was not found in the cache.

        ``key`` is a triple of ``(cmd, shell, hide_stderr)``.  ``cmd`` is
        the command tuple.  If ``shell`` is ``True``, the command is
        executed in the shell, otherwise it is executed directly.  If
        ``hide_stderr`` is ``True``, the standard error of the program is
        discarded, otherwise it is included in the output.
        """
        cmd, shell, hide_stderr = key
        proc = Popen(cmd, shell=shell, stdout=PIPE,
                     stderr=PIPE if hide_stderr else STDOUT)
        stdout = proc.communicate()[0].decode(
            sys.getfilesystemencoding()).rstrip()
        if proc.returncode != 0:
            raise CalledProcessError(proc.returncode, cmd)
        self[key] = stdout
        return stdout


def run_programs(app, doctree):
    """
    Execute all programs represented by ``program_output`` nodes in
    ``doctree``.  Each ``program_output`` node in ``doctree`` is then
    replaced with a node, that represents the output of this program.

    The program output is retrieved from the cache in
    ``app.env.programoutput_cache``.
    """
    if app.config.programoutput_use_ansi:
        # enable ANSI support, if requested by config
        from sphinxcontrib.ansi import ansi_literal_block
        node_class = ansi_literal_block
    else:
        node_class = nodes.literal_block

    cache = app.env.programoutput_cache

    for node in doctree.traverse(program_output):
        command = node['command']
        cmd_bytes = command.encode(sys.getfilesystemencoding())

        extra_args = node.get('extraargs', '').encode(
            sys.getfilesystemencoding())
        if node['use_shell']:
            cmd = cmd_bytes
            if extra_args:
                cmd += ' ' + extra_args
        else:
            cmd = shlex.split(cmd_bytes)
            if extra_args:
                cmd.extend(shlex.split(extra_args))
            cmd = tuple(cmd)

        cache_key = (cmd, node['use_shell'], node['hide_standard_error'])
        output = cache[cache_key]

        # replace lines with ..., if ellipsis is specified
        if 'strip_lines' in node:
            lines = output.splitlines()
            start, stop = node['strip_lines']
            lines[start:stop] = ['...']
            output = '\n'.join(lines)

        if node['show_prompt']:
            tmpl = app.config.programoutput_prompt_template
            output = tmpl % dict(command=command, output=output)

        new_node = node_class(output, output)
        new_node['language'] = 'text'
        node.replace_self(new_node)


def init_cache(app):
    """
    Initialize the cache for program output at
    ``app.env.programoutput_cache``, if not already present (e.g. being
    loaded from a pickled environment).

    The cache is of type :class:`ProgramOutputCache`.
    """
    if not hasattr(app.env, 'programoutput_cache'):
        app.env.programoutput_cache = ProgramOutputCache()


def setup(app):
    app.add_config_value('programoutput_use_ansi', False, 'env')
    app.add_config_value('programoutput_prompt_template',
                         '$ %(command)s\n%(output)s', 'env')
    app.add_directive('program-output', ProgramOutputDirective)
    app.add_directive('command-output', ProgramOutputDirective)
    app.connect('builder-inited', init_cache)
    app.connect('doctree-read', run_programs)
