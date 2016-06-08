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
Step
"""
from __future__ import absolute_import, division, print_function

import contextlib
from astropy.extern import six
from os.path import dirname, join, basename, splitext, abspath, split
import sys
import gc

try:
    from astropy.io import fits
    DISCOURAGED_TYPES = (fits.HDUList,)
except ImportError:
    DISCOURAGED_TYPES = None

from . import config_parser
from . import crds_client
from . import log
from . import utilities


class Step(object):
    """
    Step
    """
    spec = """
    pre_hooks = string_list(default=list())
    post_hooks = string_list(default=list())

    output_dir = string(default=None) # Directory path for output files
    output_file = output_file(default=None) # File to save output to.
    skip = boolean(default=False) # Skip this step
    """

    reference_file_types = []

    @classmethod
    def merge_config(cls, config, config_file):
        return config

    @classmethod
    def load_spec_file(cls, preserve_comments=False):
        spec = config_parser.get_merged_spec_file(
            cls, preserve_comments=preserve_comments)
        # Add arguments for all of the expected reference files
        for reference_file_type in cls.reference_file_types:
            override_name = crds_client.get_override_name(reference_file_type)
            spec[override_name] = 'string(default=None)'
            spec.inline_comments[override_name] = (
                '# Override the {0} reference file'.format(
                    reference_file_type))
        return spec

    @classmethod
    def print_configspec(cls, stream=sys.stdout):
        specfile = cls.load_spec_file(preserve_comments=True)
        specfile.write(stream)

    @classmethod
    def from_config_file(cls, config_file, parent=None, name=None):
        """
        Create a step from a configuration file.

        Parameters
        ----------
        config_file : path or readable file-like object
            The config file to load parameters from

        parent : Step instance, optional
            The parent step of this step.  Used to determine a
            fully-qualified name for this step, and to determine
            the mode in which to run this step.

        name : str, optional
            If provided, use that name for the returned instance.
            If not provided, the following are tried (in order):
                - The `name` parameter in the config file
                - The filename of the config file
                - The name of returned class

        Returns
        -------
        step : Step instance
            If the config file has a `class` parameter, the return
            value will be as instance of that class.  The `class`
            parameter in the config file must specify a subclass of
            `cls`.  If the configuration file has no `class`
            parameter, then an instance of `cls` is returned.

            Any parameters found in the config file will be set
            as member variables on the returned `Step` instance.
        """
        config = config_parser.load_config_file(config_file)

        # If a file object was passed in, pass the file name along
        if hasattr(config_file, 'name'):
            config_file = config_file.name

        step_class, name = cls._parse_class_and_name(
            config, parent, name, config_file)

        return step_class.from_config_section(
            config, parent=parent, name=name, config_file=config_file)

    @staticmethod
    def from_cmdline(args):
        """
        Create a step from a configuration file.

        Parameters
        ----------
        args : list of str
            Commandline arguments

        Returns
        -------
        step : Step instance
            If the config file has a `class` parameter, the return
            value will be as instance of that class.

            Any parameters found in the config file will be set
            as member variables on the returned `Step` instance.
        """
        from . import cmdline
        return cmdline.step_from_cmdline(args)

    @classmethod
    def _parse_class_and_name(
        cls, config, parent=None, name=None, config_file=None):
        if 'class' in config:
            step_class = utilities.import_class(config['class'],
                                                config_file=config_file)
            if not issubclass(step_class, cls):
                raise TypeError(
                    "Configuration file does not match the "
                    "expected step class.  Expected {0}, "
                    "got {1}".format(cls, step_class))
        else:
            step_class = cls

        if not name:
            name = config.get('name')
            if not name:
                if isinstance(config_file, six.string_types):
                    name = splitext(basename(config_file))[0]
                else:
                    name = step_class.__name__

        if 'name' in config:
            del config['name']
        if 'class' in config:
            del config['class']

        return step_class, name

    @classmethod
    def from_config_section(cls, config, parent=None, name=None,
                            config_file=None):
        """
        Create a step from a configuration file fragment.

        Parameters
        ----------
        config : configobj.Section instance
            The config file fragment containing parameters for this
            step only.

        parent : Step instance, optional
            The parent step of this step.  Used to determine a
            fully-qualified name for this step, and to determine
            the mode in which to run this step.

        name : str, optional
            If provided, use that name for the returned instance.
            If not provided, try the following (in order):
                - The `name` parameter in the config file fragment
                - The name of returned class

        config_file : str, optional
            The path to the config file that created this step, if
            any.  This is used to resolve relative file name
            parameters in the config file.

        Returns
        -------
        step : instance of cls
            Any parameters found in the config file fragment will be
            set as member variables on the returned `Step` instance.
        """
        if not name:
            if config.get('name'):
                name = config['name']
            else:
                name = cls.__name__

        if 'name' in config:
            del config['name']
        if 'class' in config:
            del config['class']
        if 'config_file' in config:
            del config['config_file']

        spec = cls.load_spec_file()
        config = cls.merge_config(config, config_file)
        config_parser.validate(
            config, spec, root_dir=dirname(config_file or ''))

        if 'config_file' in config:
            del config['config_file']
        if 'name' in config:
            del config['name']

        return cls(
            name=name,
            parent=parent,
            config_file=config_file,
            _validate_kwds=False,
            **config)

    def __init__(self, name=None, parent=None, config_file=None,
                 _validate_kwds=True, **kws):
        """
        Create a `Step` instance.

        Parameters
        ----------
        name : str, optional
            The name of the Step instance.  Used in logging messages
            and in cache filenames.  If not provided, one will be
            generated based on the class name.

        parent : Step instance, optional
            The parent step of this step.  Used to determine a
            fully-qualified name for this step, and to determine
            the mode in which to run this step.

        config_file : str path, optional
            The path to the config file that this step was initialized
            with.  Use to determine relative path names.

        **kws : dict
            Additional parameters to set.  These will be set as member
            variables on the new Step instance.
        """
        if _validate_kwds:
            spec = self.load_spec_file()
            kws = config_parser.config_from_dict(
                kws, spec, root_dir=dirname(config_file or ''))

        if name is None:
            name = self.__class__.__name__
        self.name = name
        if parent is None:
            self.qualified_name = '.'.join([
                log.STPIPE_ROOT_LOGGER, self.name])
        else:
            self.qualified_name = '.'.join([
                parent.qualified_name, self.name])
        self.parent = parent

        # Set the parameters as member variables
        for (key, val) in kws.items():
            setattr(self, key, val)

        # Create a new logger for this step
        self.log = log.getLogger(self.qualified_name)

        self.log.setLevel(log.logging.DEBUG)

        # Log the fact that we have been init-ed.
        self.log.info('{0} instance created.'.format(self.__class__.__name__))

        # Store the config file path so filenames can be resolved
        # against it.
        self.config_file = config_file

        if len(self.pre_hooks) or len(self.post_hooks):
            from . import hooks
            self._pre_hooks = hooks.get_hook_objects(self, 'pre', self.pre_hooks)
            self._post_hooks = hooks.get_hook_objects(self, 'post', self.post_hooks)
        else:
            self._pre_hooks = []
            self._post_hooks = []

        self._reference_files_used = []

        self._input_filename = None

    def _check_args(self, args, discouraged_types, msg):
        if discouraged_types is None:
            return

        if type(args) not in (list, tuple):
            args = [args]

        for i, arg in enumerate(args):
            if isinstance(arg, discouraged_types):
                self.log.error(
                    "{0} {1} object.  Use jwst_lib.models instead.".format(
                        msg, i))

    def run(self, *args):
        """
        Run handles the generic setup and teardown that happens with
        the running of each step.  The real work that is unique to
        each step type is done in the `process` method.
        """
        from jwst import datamodels
        gc.collect()

        # Make generic log messages go to this step's logger
        orig_log = log.delegator.log
        log.delegator.log = self.log

        result = None

        try:
            if len(args):
                self._precache_reference_files(args[0])

            self.log.info(
                'Step {0} running with args {1}.'.format(
                    self.name, args))

            for pre_hook in self._pre_hooks:
                pre_hook.run(*args)

            self._reference_files_used = []

            # Warn if passing in objects that should be
            # discouraged.
            self._check_args(args, DISCOURAGED_TYPES, "Passed")

            # Run the Step-specific code.
            if self.skip:
                self.log.info('Step skipped.')
                result = args[0]
            else:
                try:
                    result = self.process(*args)
                except TypeError as e:
                    if "process() takes exactly" in e.message:
                        raise TypeError("Incorrect number of arguments to step")
                    raise

            # Warn if returning a discouraged object
            self._check_args(result, DISCOURAGED_TYPES, "Returned")

            if not isinstance(result, (list, tuple)):
                results = [result]
            else:
                results = result

            if len(self._reference_files_used) and not self._is_association_file(args[0]):
                for result in results:
                    if isinstance(result, models.DataModel):
                        for ref_name, filename in self._reference_files_used:
                            if hasattr(result.meta.ref_file, ref_name):
                                getattr(result.meta.ref_file, ref_name).name = filename
                        result.meta.ref_file.crds.sw_version = crds_client.get_svn_version()
                        result.meta.ref_file.crds.context_used = crds_client.get_context_used()
                self._reference_files_used = []

            for post_hook in self._post_hooks:
                post_hook.run(result)

            # Save the output file if one was specified
            for i, result in enumerate(results):
                if isinstance(result, models.DataModel):
                    result.meta.calibration_software_revision = __svn_revision__
                    result.meta.calibration_software_version = __version__
                    if self.output_file is not None:
                        if len(results) > 1:
                            base, ext = splitext(self.output_file)
                            output_file_name = "{0}.{1}{2}".format(base, i, ext)
                        else:
                            output_file_name = self.output_file

                        # If the user specified an output_dir, replace the
                        # default path with output_dir
                        if self.output_dir is not None:
                            dirname, filename = split(output_file_name)
                            output_file_name = join(self.output_dir, filename)

                        self.log.info('Saving file {0}'.format(output_file_name))
                        result.save(output_file_name, clobber=True)

            self.log.info(
                'Step {0} done'.format(self.name))
        finally:
            log.delegator.log = orig_log

        return result

    __call__ = run

    def process(self, *args):
        """
        This is where real work happens. Every Step subclass has to
        override this method. The default behaviour is to raise a
        NotImplementedError exception.
        """
        raise NotImplementedError('Steps have to override process().')

    def resolve_file_name(self, file_name):
        """
        Resolve a file name expressed relative to this Step's
        configuration file.
        """
        return join(dirname(self.config_file or ''), file_name)

    @classmethod
    def call(cls, *args, **kwargs):
        """
        Make the step more conveniently callable from Python.

        To set configuration parameters, pass a `config_file` path or
        keyword arguments.

        Any positional `*args` will be passed along to the step's
        `process` method.
        """
        if 'config_file' in kwargs:
            config_file = kwargs['config_file']
            del kwargs['config_file']
            config = config_parser.load_config_file(config_file)
            auto_cls, name = cls._parse_class_and_name(config)
            config.update(kwargs)
            instance = cls.from_config_section(
                config, name=name, config_file=config_file)
        else:
            instance = cls(**kwargs)
        return instance.run(*args)

    @classmethod
    def _is_association_file(cls, input_file):
        """Return True IFF `input_file` is an association file."""
        from jwst import datamodels
        return (isinstance(input_file, str) and input_file.endswith((".asn",".json"))) or \
               isinstance(input_file, models.ModelContainer)

    def _precache_reference_files(self, input_file):
        """
        Precache all of the expected reference files before the Step's
        process method is called.
        """
        gc.collect()
        if self.skip:
            return
        if self._is_association_file(input_file):
            return
        if len(self.reference_file_types):
            from jwst import datamodels
            try:
                model = models.open(input_file)
            except (ValueError, TypeError, IOError):
                self.log.info(
                    'First argument {0} does not appear to be a '
                    'model'.format(input_file))
            else:
                self._precache_reference_files_impl(model)
                model.close()
        gc.collect()

    @classmethod
    def list_reference_files(cls, input_file):
        """
        List reference types and files name for a particular input file.
        """
        if cls._is_association_file(input_file):
            return
        return crds_client.get_multiple_reference_paths(
            input_file, cls.reference_file_types)

    def _precache_reference_files_impl(self, model):
        """Given open data `model`,  determine and cache reference files for
        any reference types which are not overridden on the command line.

        Verify that all CRDS and overridden reference files are readable.
        """
        if self.skip:
            return
        if self._is_association_file(model):
            return
        ovr_refs = {
            reftype : self._get_ref_override(reftype)
            for reftype in self.reference_file_types
            if self._get_ref_override(reftype) is not None
            }
        fetch_types = sorted(set(self.reference_file_types) - set(ovr_refs.keys()))
        crds_refs = crds_client.get_multiple_reference_paths(model, fetch_types)
        ref_path_map = dict(list(crds_refs.items()) + list(ovr_refs.items()))
        for (reftype, refpath) in sorted(ref_path_map.items()):
            how = "Override" if reftype in ovr_refs else "Prefetch"
            self.log.info("{0} for {1} reference file is '{2}'.".format(how, reftype.upper(), refpath))
            crds_client.check_reference_open(refpath)

    def _get_ref_override(self, reference_file_type):
        """Determine and return any override for `reference_file_type`.

        Returns
        -------
        override_filepath or None.
        """
        override_name = crds_client.get_override_name(reference_file_type)
        path = getattr(self, override_name, None)
        return abspath(path) if path else path

    def get_reference_file(self, input_file, reference_file_type):
        """
        Get a reference file from CRDS.  If the configuration file or
        commandline parameters override the reference file, it will be
        automatically used when calling this function.

        Parameters
        ----------
        input_file : jwst_lib.models.ModelBase instance
            A model of the input file.  Metadata on this input file
            will be used by the CRDS "bestref" algorithm to obtain a
            reference file.

        reference_file_type : string
            The type of reference file to retrieve.  For example, to
            retrieve a flat field reference file, this would be
            'flat_field'.

        Returns
        -------
        reference_file : readable file-like object
            A readable file-like object with the contents of the
            reference file.
        """
        override = self._get_ref_override(reference_file_type)
        if override is not None:
            if override.strip() != "":
                self._reference_files_used.append(
                    (reference_file_type, abspath(override)))
                reference_name = override
            else:
                return ""
        else:
            reference_name = crds_client.get_reference_file(
                input_file, reference_file_type)
            gc.collect()
            if reference_name != "N/A":
                hdr_name = "crds://" + basename(reference_name)
            else:
                hdr_name = "N/A"
            self._reference_files_used.append(
                (reference_file_type, hdr_name))
        return crds_client.check_reference_open(reference_name)

    @contextlib.contextmanager
    def get_reference_file_model(self, input_file, reference_file_type):
        """
        Get a reference file from CRDS as a jwst_lib.models.ModelBase
        object.  If the configuration file or commandline parameters
        override the reference file, it will be automatically used
        when calling this function.

        Parameters
        ----------
        input_file : jwst_lib.models.ModelBase instance
            A model of the input file.  Metadata on this input file
            will be used by the CRDS "bestref" algorithm to obtain a
            reference file.

        reference_file_type : string
            The type of reference file to retrieve.  For example, to
            retrieve a flat field reference file, this would be
            'flat_field'.

        Returns
        -------
        reference_file_model : jwst_lib.models.ModelBase instance
            A model to access the contents of the reference file.
        """
        from jwst import datamodels

        filename = self.get_reference_file(input, reference_file_type)
        with models.open(filename) as model:
            yield model
        gc.collect()

    def set_input_filename(self, path):
        """
        Sets the name of the master input file.  Used to generate output
        file names.
        """
        self._input_filename = path

    def save_model(self, model, name, *args, **kwargs):
        """
        Saves the given model using a filename based on the model's name
        or the input filename.

        The "root" or "source" filename is determined first by looking
        at `meta.filename` on the passed-in model.  If that entry
        doesn't exist, the input filename passed to stpipe on the
        command line (or in the config file) is used.

        To generate an output filename from the root filename, the
        part of the input filename following the last underscore is
        removed and replaced with `name`, and the same filename
        extension is used.

        Parameters
        ----------
        model : jwst_lib.models.Model instance
            The model to save.

        name : str
            The suffix to add to the filename.

        """
        if model.meta.filename is not None:
            root = model.meta.filename
        elif self._input_filename is not None:
            root = self._input_filename
        else:
            raise ValueError(
                "Model has no filename, and step has no input filename")

        dirname, filename = split(root)
        base, ext = splitext(filename)
        new_filename = base[:base.rfind('_')] + '_' + name + ext

        # If the user specified an output_dir, replace the
        # original dirname with output_dir
        if self.output_dir is not None:
            dirname = self.output_dir

        new_path = join(dirname, new_filename)
        model.save(new_path, *args, **kwargs)
