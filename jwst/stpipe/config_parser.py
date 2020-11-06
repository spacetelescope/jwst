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
Our configuration files are ConfigObj/INI files.
"""
from inspect import isclass
import logging
import os
import os.path
import textwrap

from asdf import open as asdf_open
from asdf import ValidationError as AsdfValidationError
from stdatamodels import DataModel
from stdatamodels import s3_utils

from ..extern.configobj.configobj import (
    ConfigObj, Section, flatten_errors, get_extra_values)
from ..extern.configobj.validate import Validator, ValidateError, VdtTypeError

from ..datamodels import StepParsModel

from . import utilities


# Configure logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class ValidationError(Exception):
    pass


def _get_input_file_check(root_dir):
    from . import cmdline

    root_dir = root_dir or ''

    def _input_file_check(path):
        if not isinstance(path, cmdline.FromCommandLine):
            try:
                path = str(path)
            except ValueError:
                pass

            path = os.path.join(root_dir, path)

        path = os.path.abspath(path)
        if not os.path.exists(path):
            raise ValidateError(
                "Path {0!r} does not exist".format(path))

        return path

    return _input_file_check


def _get_output_file_check(root_dir):
    from . import cmdline

    root_dir = root_dir or ''

    def _output_file_check(path):
        if not isinstance(path, cmdline.FromCommandLine):
            try:
                path = str(path)
            except ValueError:
                pass

            path = os.path.join(root_dir, path)

        path = os.path.abspath(path)
        dir = os.path.dirname(path)
        if dir and not os.path.exists(dir):
            os.makedirs(dir)

        return path

    return _output_file_check


def _is_datamodel(value, default=None):
    """Verify that value is either is a DataModel."""
    if isinstance(value, DataModel):
        return value
    else:
        raise VdtTypeError(value)


def _is_string_or_datamodel(value, default=None):
    """Verify that value is either a string (nominally a reference file path)
    or a DataModel (possibly one with no corresponding file.)
    """
    if isinstance(value, DataModel):
        return value
    elif isinstance(value, str):
        return value
    else:
        raise VdtTypeError(value)


def load_config_file(config_file):
    """
    Read the file `config_file` and return the parsed configuration.
    """
    if s3_utils.is_s3_uri(config_file):
        return _load_config_file_s3(config_file)
    else:
        return _load_config_file_filesystem(config_file)


def _load_config_file_filesystem(config_file):
    if not os.path.isfile(config_file):
        raise ValueError("Config file {0} not found.".format(config_file))
    try:
        cfg = asdf_open(config_file)
    except (AsdfValidationError, ValueError):
        logger.debug('Config file did not parse as ASDF. Trying as ConfigObj: %s', config_file)
        return ConfigObj(config_file, raise_errors=True)

    # Seems to be ASDF. Create the configobj from that.
    return _config_obj_from_asdf(cfg)


def _load_config_file_s3(config_file):
    if not s3_utils.object_exists(config_file):
        raise ValueError("Config file {0} not found.".format(config_file))

    content = s3_utils.get_object(config_file)
    try:
        cfg = asdf_open(content)
    except (AsdfValidationError, ValueError):
        logger.debug('Config file did not parse as ASDF. Trying as ConfigObj: %s', config_file)
        content.seek(0)
        return ConfigObj(content, raise_errors=True)

    # Seems to be ASDF. Create the configobj from that.
    return _config_obj_from_asdf(cfg)


def _config_obj_from_asdf(cfg):
    configobj = ConfigObj()
    configobj.merge(cfg['parameters'])
    configobj.pars_model = StepParsModel(cfg)
    cfg.close()
    return configobj


def get_merged_spec_file(cls, preserve_comments=False):
    """
    Creates a merged spec file for a Step class and all of its
    subclasses.

    Parameters
    ----------
    cls: `Step`-derived or `Step` instance
        A class or instance of a `Step`-based class.

    preserve_comments : bool, optional
        When True, preserve the comments in the spec file
    """
    if not isclass(cls):
        cls = cls.__class__
    subclasses = cls.mro()
    subclasses.reverse()

    config = ConfigObj()
    cfg = None
    for subclass in subclasses:
        cfg = load_spec_file(subclass, preserve_comments=preserve_comments)
        if cfg:
            merge_config(config, cfg)

    if cfg is not None:
        config.initial_comment = cfg.initial_comment
        config.final_comment = cfg.final_comment

    return config


def load_spec_file(cls, preserve_comments=False):
    """
    Load the spec file corresponding to the given class.

    Parameters
    ----------
    cls: `Step`-derived or `Step` instance
        A class or instance of a `Step`-based class.

    preserve_comments: bool
        True to keep comments in the resulting `ConfigObj`

    Returns
    -------
    spec_file: ConfigObj
        The resulting configuration object
    """
    # Don't use 'hasattr' here, because we don't want to inherit spec
    # from the base class.
    if not isclass(cls):
        cls = cls.__class__
    if 'spec' in cls.__dict__:
        spec = cls.spec.strip()
        spec_file = textwrap.dedent(spec)
        spec_file = spec_file.split('\n')
        encoded = []
        for line in spec_file:
            if isinstance(line, str):
                encoded.append(line.encode('utf8'))
            else:
                encoded.append(line)
        spec_file = encoded
    else:
        spec_file = utilities.find_spec_file(cls)
    if spec_file:
        return ConfigObj(spec_file, _inspec=not preserve_comments,
                         raise_errors=True)
    return


def merge_config(into, new):
    """
    Merges a configuration tree into another one.

    Unlike merge in configobj itself, this one updates inline
    comments.

    Parameters
    ----------
    into : `configobj.Section`
        The configuration tree to merge into

    new : `configobj.Section`
        The source of new configuration values
    """
    for key, val in new.items():
        if isinstance(val, Section):
            if key not in into:
                section = Section(
                    into, into.depth + 1, into.main, name=key)
                into[key] = section
            merge_config(into[key], val)
        elif key not in new.defaults:
            into[key] = val
            into.inline_comments[key] = new.inline_comments[key]
            into.comments[key] = new.comments[key]


def config_from_dict(d, spec=None, root_dir=None, allow_missing=False):
    """
    Create a ConfigObj from a dict.

    Parameters
    ----------
    d: dict
        The dictionary to merge into the resulting ConfigObj.

    spec: ConfigObj
        The specification to validate against.
        If None, just convert dictionary into a ConfigObj.

    root_dir: str
        The base directory to use for file-based parameters.

    allow_missing: bool
        If a parameter is not defined and has no default in the spec,
        set that parameter to its specification.
    """
    config = ConfigObj()

    config.update(d)

    if spec:
        validate(config, spec, root_dir=root_dir, allow_missing=allow_missing)
    else:
        config.walk(string_to_python_type)

    return config


def validate(config, spec, section=None, validator=None, root_dir=None, allow_missing=False):
    """
    Parse config_file, in INI format, and do validation with the
    provided specfile.

    Parameters
    ----------
    config: ConfigObj
        The configuration to validate.

    spec: ConfigObj
        The specification to validate against.

    section: ConfigObj or None
        The specific section of config to validate.
        If None, then all sections are validated.

    validator: extern.configobj.validator.Validator or None
        The validator to use. If None, the default will be used.

    root_dir: str
        The directory to use as the basis for any file-based parameters.

    allow_missing: bool
        If a parameter is not defined and has no default in the spec,
        set that parameter to its specification.
    """
    if spec is None:
        config.walk(string_to_python_type)
        return config

    if validator is None:
        validator = Validator()
        validator.functions['input_file'] = _get_input_file_check(root_dir)
        validator.functions['output_file'] = _get_output_file_check(root_dir)
        validator.functions['is_datamodel'] = _is_datamodel
        validator.functions['is_string_or_datamodel'] = _is_string_or_datamodel

    orig_configspec = config.main.configspec
    config.main.configspec = spec

    try:
        if config.main != config:
            section = config
        else:
            section = None

        errors = config.main.validate(
            validator, preserve_errors=True,
            section=section)

        messages = []
        if errors is not True:
            for section_list, key, err in flatten_errors(config, errors):
                if key is not None:
                    section_list.append(key)
                else:
                    section_list.append('[missing section]')
                section_string = '/'.join(section_list)
                if err == False:
                    if allow_missing:
                        config[key] = spec[key]
                        continue
                    else:
                        err = 'missing'

                messages.append(
                    "Config parameter {0!r}: {1}".format(
                        section_string, err))

        extra_values = get_extra_values(config)
        if extra_values:
            sections, name = extra_values[0]
            if len(sections) == 0:
                sections = 'root'
            else:
                sections = '/'.join(sections)
            messages.append(
                "Extra value {0!r} in {1}".format(
                    name, sections))

        if len(messages):
            raise ValidationError('\n'.join(messages))
    finally:
        config.main.configspec = orig_configspec

    return config


def string_to_python_type(section, key):
    """
    Do blind type inferring.
    """
    # We parse scalars and lists.
    val = section[key]
    if isinstance(val, list):
        typed_val = [_parse(x) for x in val]
    else:
        typed_val = _parse(val)
    section[key] = typed_val
    return


def _parse(val):
    """
    Parse scalar strings into scalar python types.
    """
    if val.lower() == 'true':
        return True
    elif val.lower() == 'false':
        return False
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    return str(val)
