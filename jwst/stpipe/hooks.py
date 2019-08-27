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
Pre- and post-hooks
"""
import types

from . import Step
from . import utilities

def hook_from_string(step, type, num, command):
    name = '{0}_hook{1:d}'.format(type, num)

    step_class = None
    try:
        step_class = utilities.import_class(
            command, Step, step.config_file)
    except Exception:
        pass

    if step_class is not None:
        return step_class(
            name, parent=step, config_file=step.config_file)

    step_func = None
    try:
        step_func = utilities.import_class(
            command, types.FunctionType, step.config_file)
    except Exception:
        pass

    if step_func is not None:
        from . import function_wrapper
        return function_wrapper.FunctionWrapper(
            step_func, parent=step, config_file=step.config_file)

    from .subproc import SystemCall

    return SystemCall(name, parent=step, command=command)


def get_hook_objects(step, type, hooks):
    return [hook_from_string(step, type, i, hook)
            for i, hook in enumerate(hooks)]
