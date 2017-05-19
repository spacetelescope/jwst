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
Pipeline

"""

from __future__ import absolute_import, division, print_function

import gc
from os.path import dirname, join, split, splitext
import re

from .configobj.configobj import Section

from . import config_parser
from . import Step

REMOVE_SUFFIX = '(.+?)(_(rate|cal)(ints)?)?$'


class Pipeline(Step):
    """
    A Pipeline is a way of combining a number of steps together.
    """

    # Configuration
    spec = """
    output_basename = string(default=None)    # Output base name
    output_ext = string(default=".fits")      # Output extension
    suffix = string(default=None)             # Suffix for output file name
    output_use_model = boolean(default=False) # force use `meta.filename` as the output name
    """
    # A set of steps used in the Pipeline.  Should be overridden by
    # the subclass.
    step_defs = {}

    def __init__(self, *args, **kwargs):
        """
        See `Step.__init__` for the parameters.
        """
        Step.__init__(self, *args, **kwargs)

        # Configure all of the steps
        for key, val in self.step_defs.items():
            cfg = self.steps.get(key)
            if cfg is not None:
                new_step = val.from_config_section(
                    cfg, parent=self, name=key,
                    config_file=self.config_file)
            else:
                new_step = val(
                    key, parent=self, config_file=self.config_file,
                    **kwargs.get(key, {}))

            setattr(self, key, new_step)

    @classmethod
    def merge_config(cls, config, config_file):
        steps = config.get('steps', {})

        # Configure all of the steps
        for key, val in cls.step_defs.items():
            cfg = steps.get(key)
            if cfg is not None:
                # If a config_file is specified, load those values and
                # then override them with our values.
                if cfg.get('config_file'):
                    cfg2 = config_parser.load_config_file(
                        join(dirname(config_file or ''), cfg.get('config_file')))
                    del cfg['config_file']
                    config_parser.merge_config(cfg2, cfg)
                    steps[key] = cfg2

        return config

    @classmethod
    def load_spec_file(cls, preserve_comments=False):
        spec = config_parser.get_merged_spec_file(
            cls, preserve_comments=preserve_comments)

        spec['steps'] = Section(spec, spec.depth + 1, spec.main, name="steps")
        steps = spec['steps']
        for key, val in cls.step_defs.items():
            if not issubclass(val, Step):
                raise TypeError(
                    "Entry {0!r} in step_defs is not a Step subclass"
                    .format(key))
            stepspec = val.load_spec_file(preserve_comments=preserve_comments)
            steps[key] = Section(steps, steps.depth + 1, steps.main, name=key)

            config_parser.merge_config(steps[key], stepspec)

            # Also add a key that can be used to specify an external
            # config_file
            step = spec['steps'][key]
            step['config_file'] = 'string(default=None)'
            step['name'] = "string(default='')"
            step['class'] = "string(default='')"

        return spec

    def _precache_reference_files(self, input_file):
        """
        Precache all of the expected reference files in this Pipeline
        and all of its constituent Steps process method is called.
        """
        from .. import datamodels
        gc.collect()
        if self._is_association_file(input_file):
            return
        try:
            with datamodels.open(input_file) as model:
                pass
        except (ValueError, TypeError, IOError):
            self.log.info(
                'First argument {0} does not appear to be a '
                'model'.format(input_file))
        else:
            super(Pipeline, self)._precache_reference_files(input_file)
            for name in self.step_defs.keys():
                step = getattr(self, name)
                step._precache_reference_files(input_file)
        gc.collect()

    def set_input_filename(self, path):
        self._input_filename = path
        for key, val in self.step_defs.items():
            getattr(self, key).set_input_filename(path)

    @staticmethod
    def make_output_path(
            step, data,
            basepath=None, suffix=None, ext=None,
            result_id=None, ignore_use_model=False
    ):
        """Make up a path based on data and user specification

        Parameters
        ----------
        step: Step
            The step which produced the data

        data: obj
            Unused by this routine

        basepath: str or None
            The output file name. If `None` or empty string, create
            a filename based on the data.

        suffix: str or None
            The suffix to append to the basename.

        ext: str or None
            The file format extension

        result_id: str or None
            If a suffix cannot be determined, use this as the suffix.
            If the result is still None, raise ValueError

        ignore_use_model: bool
            Ignore configuration parameter `output_use_model`

        Returns
        -------
        output_path: str
            The fully qualified output path
        """
        from ..datamodels import DataModel

        has_basepath = basepath is not None and len(basepath) > 0
        use_model_name = not ignore_use_model and getattr(step, 'output_use_model', False)
        output_path = basepath

        # If a basepath was specified, use that, but adding the
        # result_id if necessary
        if has_basepath:
            path, filename = split(output_path)
            name, filename_ext = splitext(filename)
            output_name = [name]
            suffix = _get_suffix(suffix, step, default_suffix=result_id)
            if suffix is not None:
                output_name.append('_' + suffix)
            output_name.append(filename_ext)
            output_name = ''.join(output_name)

        # Otherwise, construct a name
        else:

            # Make names based on DataModels
            if isinstance(data, DataModel):

                # If using what is in the model, just retrieve that.
                if use_model_name:
                    basepath = data.meta.filename
                    path, output_name = split(basepath)

                # Otherwise, create a fully qualified pipeline outputname
                else:
                    basepath = step.search_attr('output_basename')
                    if basepath is None:
                        basepath = data.meta.filename

                    # Breakdown the components
                    path, filename = split(basepath)
                    name, filename_ext = splitext(filename)

                    # Remove any known, previous suffixes.
                    match = re.match(REMOVE_SUFFIX, name)
                    name = match.group(1)

                    # Rebuild the path.
                    output_name = [name]
                    suffix = _get_suffix(suffix, step, default_suffix=result_id)
                    if suffix is not None:
                        output_name.append('_' + suffix)
                    if ext is None:
                        ext = step.search_attr('output_ext')
                        if ext is None:
                            ext = filename_ext
                    if ext is not None:
                        output_name.append(ext)
                    output_name = ''.join(output_name)

        output_path = output_name
        output_dir = step.search_attr('output_dir')
        if output_dir is not None:
            output_path = join(output_dir, output_path)
        return output_path


# #########
# Utilities
# #########
def _get_suffix(suffix, step, default_suffix=None):
    """Retrieve either specified or pipeline-supplied suffix

    Parameters
    ----------
    suffix: str or None
        Suffix to use if specified.

    step: Step or Pipeline
        The step to retrieve the suffux.

    default_suffix: str
        If the pipeline does not supply a suffix,
        use this.

    Returns
    -------
    suffix: str or None
        Suffix to use
    """
    if suffix is None:
        suffix = step.search_attr('suffix')
    if suffix is None:
        suffix = default_suffix
    return suffix
