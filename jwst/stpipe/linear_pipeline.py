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
LinearPipeline

"""
import gc

from .pipeline import Pipeline


class _LinearPipelineMetaclass(type):
    def __init__(cls, name, bases, dct):
        super(_LinearPipelineMetaclass, cls).__init__(name, bases, dct)
        pipeline_steps = cls.pipeline_steps
        if pipeline_steps is not None and len(pipeline_steps) == 0:
            raise ValueError(
                "{0!r} LinearPipeline subclass defines no pipeline_steps"
                .format(name))
        if pipeline_steps is None:
            pipeline_steps = []
        cls.step_defs = dict(pipeline_steps)

# Since the pipeline_steps member needs to be converted to a step_defs
# at the class level, we need to use a metaclass.

class LinearPipeline(Pipeline, metaclass=_LinearPipelineMetaclass):
    """
    A LinearPipeline is a way of combining a number of steps together
    in a simple linear order.
    """

    spec = """
    # start_step and end_step allow only a part of the pipeline to run
    start_step = string(default=None)  # Start the pipeline at this step
    end_step = string(default=None)    # End the pipeline right before this step

    # [steps] section is implicitly added by the Pipeline class.
    """

    # To be overridden by subclasses
    pipeline_steps = None

    def _check_start_and_end_steps(self):
        """
        Given the start_step and end_step members (which are strings
        or None), find the actual step objects they correspond to.
        """
        start_step = end_step = None

        if self.start_step is not None:
            if hasattr(self, self.start_step):
                start_step = getattr(self, self.start_step)
            else:
                raise ValueError(
                    "start_step {0!r} not found".format(
                        self.start_step))

        if self.end_step is not None:
            if hasattr(self, self.end_step):
                end_step = getattr(self, self.end_step)
            else:
                raise ValueError(
                    "end_step {0!r} not found".format(
                        self.end_step))

        return start_step, end_step

    def process(self, input_file):
        """
        Run the pipeline.
        """
        self._check_start_and_end_steps()

        do_caching = (
            self.end_step is not None and
            self.end_step != self.pipeline_steps[-1][0])

        if self.start_step is None:
            mode = 'RUN'
        else:
            mode = 'BEFORE'

        # It would be easiest to do this in a loop,
        # but we use recursion instead to make the "with" statements
        # work correctly

        def recurse(mode, input_file, pipeline_steps):
            gc.collect()
            if pipeline_steps == []:
                if (hasattr(self, 'output_file') and
                    self.output_file is not None):
                    input_file.save(self.output_file)
                return input_file

            name, cls = pipeline_steps[0]
            step = getattr(self, name)
            filename = '{0}.fits'.format(self.qualified_name)
            if name == self.start_step:
                mode = 'RUN'

            if mode == 'BEFORE':
                from .. import datamodels

                try:
                    with datamodels.open(filename) as dm:
                        pass
                except (ValueError, TypeError, IOError):
                    return recurse(mode, filename, pipeline_steps[1:])
                else:
                    dm = datamodels.open(filename)
                    return recurse(mode, dm, pipeline_steps[1:])

            elif mode == 'RUN':
                dm = step(input_file)
                if do_caching:
                    dm.save(filename)
                if name == self.end_step:
                    return None
                return recurse(mode, dm, pipeline_steps[1:])

        result = recurse(mode, input_file, self.pipeline_steps)
        gc.collect()
        return result

    def set_input_filename(self, path):
        for name, cls in self.pipeline_steps:
            getattr(self, name).set_input_filename(path)
