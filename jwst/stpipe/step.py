"""
Step
"""
from contextlib import contextmanager
from functools import partial
import gc
import inspect
import os
from os.path import (
    abspath,
    basename,
    dirname,
    expanduser,
    expandvars,
    isfile,
    join,
    split,
    splitext,
)
import sys

try:
    from astropy.io import fits
    DISCOURAGED_TYPES = (fits.HDUList,)
except ImportError:
    DISCOURAGED_TYPES = None
from stdatamodels import DataModel

from . import config_parser
from . import crds_client
from . import log
from . import utilities
from .. import __version_commit__, __version__
from ..associations.load_as_asn import (LoadAsAssociation, LoadAsLevel2Asn)
from ..associations.lib.format_template import FormatTemplate
from ..associations.lib.update_path import update_key_value
from ..datamodels import (ModelContainer, StepParsModel)
from ..datamodels import open as dm_open
from ..lib.class_property import ClassInstanceMethod
from ..lib.suffix import remove_suffix

class Step():
    """
    Step
    """
    spec = """
    pre_hooks          = string_list(default=list())
    post_hooks         = string_list(default=list())
    output_file        = output_file(default=None)   # File to save output to.
    output_dir         = string(default=None)        # Directory path for output files
    output_ext         = string(default='.fits')     # Default type of output
    output_use_model   = boolean(default=False)      # When saving use `DataModel.meta.filename`
    output_use_index   = boolean(default=True)       # Append index.
    save_results       = boolean(default=False)      # Force save results
    skip               = boolean(default=False)      # Skip this step
    suffix             = string(default=None)        # Default suffix for output files
    search_output_file = boolean(default=True)       # Use outputfile define in parent step
    input_dir          = string(default=None)        # Input directory
    """

    # Correction parameters. These store and use whatever information a Step
    # may need to perform its operations without re-calculating, or to use
    # from a previous run of the Step.  The structure is up to each Step.
    correction_pars = None
    use_correction_pars = False

    # Reference types for both command line override
    # definition and reference prefetch
    reference_file_types = []

    # Set to False in subclasses to skip prefetch,
    # but by default attempt to prefetch
    prefetch_references = True

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
            spec[override_name] = 'is_string_or_datamodel(default=None)'
            spec.inline_comments[override_name] = (
                '# Override the {0} reference file'.format(
                    reference_file_type))
        return spec

    @classmethod
    def print_configspec(cls):
        specfile = cls.load_spec_file(preserve_comments=True)
        specfile.write(sys.stdout.buffer)

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
    def _parse_class_and_name(cls, config, parent=None, name=None, config_file=None):
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
                if isinstance(config_file, str):
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
            - The ``name`` parameter in the config file fragment
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

        step = cls(
            name=name,
            parent=parent,
            config_file=config_file,
            _validate_kwds=False,
            **config)

        return step

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
            with.  Use to determine relative path names of other config files.

        **kws : dict
            Additional parameters to set.  These will be set as member
            variables on the new Step instance.
        """
        self._pars_model = None
        self._reference_files_used = []
        self._input_filename = None
        self._input_dir = None
        self._keywords = kws
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

        # Store the config file path so config filenames can be resolved
        # against it.
        self.config_file = config_file

        # Setup the hooks
        if len(self.pre_hooks) or len(self.post_hooks):
            from . import hooks
            self._pre_hooks = hooks.get_hook_objects(
                self, 'pre', self.pre_hooks
            )
            self._post_hooks = hooks.get_hook_objects(
                self, 'post', self.post_hooks
            )
        else:
            self._pre_hooks = []
            self._post_hooks = []

    def _check_args(self, args, discouraged_types, msg):
        if discouraged_types is None:
            return

        if type(args) not in (list, tuple):
            args = [args]

        for i, arg in enumerate(args):
            if isinstance(arg, discouraged_types):
                self.log.error(
                    "{0} {1} object.  Use jwst.datamodels instead.".format(
                        msg, i))

    def run(self, *args):
        """
        Run handles the generic setup and teardown that happens with
        the running of each step.  The real work that is unique to
        each step type is done in the `process` method.
        """
        from .. import datamodels
        gc.collect()

        # Make generic log messages go to this step's logger
        orig_log = log.delegator.log
        log.delegator.log = self.log

        step_result = None

        self.log.info(
            'Step {0} running with args {1}.'.format(
                self.name, args))

        self.log.info(
            f'Step {self.name} parameters are: {self.get_pars()}'
        )

        if len(args):
            self.set_primary_input(args[0])

        try:
            # Default output file configuration
            if self.output_file is not None:
                self.save_results = True

            if self.suffix is None:
                self.suffix = self.default_suffix()

            hook_args = args
            for pre_hook in self._pre_hooks:
                hook_results = pre_hook.run(*hook_args)
                if hook_results is not None:
                    hook_args = hook_results
            args = hook_args

            self._reference_files_used = []

            # Warn if passing in objects that should be
            # discouraged.
            self._check_args(args, DISCOURAGED_TYPES, "Passed")

            # Run the Step-specific code.
            if self.skip:
                self.log.info('Step skipped.')
                step_result = args[0]
            else:
                if self.prefetch_references:
                    self.prefetch(*args)
                try:
                    step_result = self.process(*args)
                except TypeError as e:
                    if "process() takes exactly" in str(e):
                        raise TypeError(
                            "Incorrect number of arguments to step"
                        )
                    raise

            # Warn if returning a discouraged object
            self._check_args(step_result, DISCOURAGED_TYPES, "Returned")

            # Run the post hooks
            for post_hook in self._post_hooks:
                hook_results = post_hook.run(step_result)
                if hook_results is not None:
                    step_result = hook_results

            # Update meta information
            if not isinstance(
                    step_result, (list, tuple, datamodels.ModelContainer)
            ):
                results = [step_result]
            else:
                results = step_result

            if len(self._reference_files_used):
                for result in results:
                    if isinstance(result, DataModel):
                        for ref_name, filename in self._reference_files_used:
                            if hasattr(result.meta.ref_file, ref_name):
                                getattr(result.meta.ref_file, ref_name).name = filename
                        result.meta.ref_file.crds.sw_version = crds_client.get_svn_version()
                        result.meta.ref_file.crds.context_used = \
                            crds_client.get_context_used(result.meta.telescope)
                self._reference_files_used = []

            # Mark versions
            for result in results:
                if isinstance(result, DataModel):
                    result.meta.calibration_software_revision = __version_commit__ or 'RELEASE'
                    result.meta.calibration_software_version = __version__

            # Save the output file if one was specified
            if not self.skip and self.save_results:

                # Setup the save list.
                if not isinstance(step_result, (list, tuple)):
                    results_to_save = [step_result]
                else:
                    results_to_save = step_result

                for idx, result in enumerate(results_to_save):
                    if len(results_to_save) <= 1:
                        idx = None
                    if isinstance(result, DataModel):
                        self.save_model(result, idx=idx)
                    elif hasattr(result, 'save'):
                        try:
                            output_path = self.make_output_path(idx=idx)
                        except AttributeError:
                            self.log.warning(
                                '`save_results` has been requested,'
                                ' but cannot determine filename.'
                            )
                            self.log.warning(
                                'Specify an output file with `--output_file`'
                                ' or set `--save_results=false`'
                            )
                        else:
                            self.log.info(
                                'Saving file {0}'.format(output_path)
                            )
                            result.save(output_path, overwrite=True)

            self.log.info(
                'Step {0} done'.format(self.name))
        finally:
            log.delegator.log = orig_log

        return step_result

    __call__ = run

    def prefetch(self, *args):
        """Prefetch reference files,  nominally called when
        self.prefetch_references is True.  Can be called explictly
        when self.prefetch_refences is False.
        """
        # prefetch truly occurs at the Pipeline (or subclass) level.
        if len(args) and len(self.reference_file_types) and not self.skip:
            self._precache_references(args[0])

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
        Creates and runs a new instance of the class.

        Gets a config file from CRDS if one is available

        To set configuration parameters, pass a `config_file` path or
        keyword arguments.  Keyword arguments override those in the
        specified `config_file`.

        Any positional `*args` will be passed along to the step's
        `process` method.

        Note: this method creates a new instance of `Step` with the given
        `config_file` if supplied, plus any extra `*args` and `**kwargs`.
        If you create an instance of a Step, set parameters, and then use
        this `call()` method, it will ignore previously-set parameters, as
        it creates a new instance of the class with only the `config_file`,
        `*args` and `**kwargs` passed to the `call()` method.

        If not used with a `config_file` or specific `*args` and `**kwargs`,
        it would be better to use the `run` method, which does not create
        a new instance but simply runs the existing instance of the `Step`
        class.
        """
        logger_name = cls.get_pars_model().instance['parameters']['name']
        log_cls = log.getLogger(logger_name)
        if len(args) > 0:
            filename = args[0]
            crds_config = cls.get_config_from_reference(filename)
        else:
            log_cls.info("No filename given, cannot retrieve config from CRDS")
            crds_config = config_parser.ConfigObj()
        if 'config_file' in kwargs:
            config_file = kwargs['config_file']
            del kwargs['config_file']
            config_from_file = config_parser.load_config_file(config_file)
            config_parser.merge_config(crds_config, config_from_file)
        else:
            config_file = None

        crds_config.update(kwargs)

        if 'class' in crds_config:
            del crds_config['class']

        name = crds_config.get('name', None)
        instance = cls.from_config_section(crds_config,
            name=name, config_file=config_file)

        return instance.run(*args)

    @property
    def input_dir(self):
        return self.search_attr('_input_dir', '')

    @input_dir.setter
    def input_dir(self, input_dir):
        self._input_dir = input_dir

    def default_output_file(self, input_file=None):
        """Create a default filename based on the input name"""
        output_file = input_file
        if output_file is None or not isinstance(output_file, str):
                output_file = self.search_attr('_input_filename')
        if output_file is None:
            output_file = 'step_{}{}'.format(
                self.name,
                self.output_ext
            )
        return output_file

    def default_suffix(self):
        """Return a default suffix based on the step"""
        return self.name.lower()

    def search_attr(self, attribute, default=None, parent_first=False):
        """Return first non-None attribute in step heirarchy

        Parameters
        ----------
        attribute : str
            The attribute to retrieve

        default : obj
            If attribute is not found, the value to use

        parent_first : bool
            If `True`, allow parent definition to override step version

        Returns
        -------
        value : obj
            Attribute value or `default` if not found
        """
        if parent_first:
            try:
                value = self.parent.search_attr(
                    attribute, parent_first=parent_first
                )
            except AttributeError:
                value = None
            if value is None:
                value = getattr(self, attribute, default)
            return value
        else:
            value = getattr(self, attribute, None)
            if value is None:
                try:
                    value = self.parent.search_attr(attribute)
                except AttributeError:
                    pass
            if value is None:
                value = default
            return value

    def _precache_references(self, input_file):
        """Because Step precaching precedes calls to get_reference_file() almost
        immediately, true precaching has been moved to Pipeline where the
        interleaving of precaching and Step processing is more of an
        issue. This null method is intended to be overridden in Pipeline by
        true precache operations and avoids having to override the more complex
        Step.run() instead.
        """
        pass

    def get_ref_override(self, reference_file_type):
        """Determine and return any override for `reference_file_type`.

        Returns
        -------
        override_filepath or None.
        """
        override_name = crds_client.get_override_name(reference_file_type)
        path = getattr(self, override_name, None)
        if isinstance(path, DataModel):
            return path
        else:
            return abspath(path) if path and path != 'N/A' else path

    def get_reference_file(self, input_file, reference_file_type):
        """
        Get a reference file from CRDS.

        If the configuration file or commandline parameters override the
        reference file, it will be automatically used when calling this
        function.

        Parameters
        ----------
        input_file : jwst.datamodels.ModelBase instance
            A model of the input file.  Metadata on this input file
            will be used by the CRDS "bestref" algorithm to obtain a
            reference file.

        reference_file_type : string
            The type of reference file to retrieve.  For example, to
            retrieve a flat field reference file, this would be 'flat'.

        Returns
        -------
        reference_file : path of reference file,  a string
        """
        override = self.get_ref_override(reference_file_type)
        if override is not None:
            if isinstance(override, DataModel):
                self._reference_files_used.append(
                    (reference_file_type, override.override_handle))
                return override
            elif override.strip() != "":
                self._reference_files_used.append(
                    (reference_file_type, basename(override)))
                reference_name = override
            else:
                return ""
        else:
            reference_name = crds_client.get_reference_file(
                input_file, reference_file_type)
            if reference_name != "N/A":
                hdr_name = "crds://" + basename(reference_name)
            else:
                hdr_name = "N/A"
            self._reference_files_used.append(
                (reference_file_type, hdr_name))
        return crds_client.check_reference_open(reference_name)

    @classmethod
    def get_config_from_reference(cls, dataset, observatory=None, disable=None):
        """Retrieve step parameters from reference database

        Parameters
        ----------
        cls : `jwst.stpipe.step.Step`
            Either a class or instance of a class derived
            from `Step`.
        dataset : `jwst.datamodels.ModelBase`
            A model of the input file.  Metadata on this input file will
            be used by the CRDS "bestref" algorithm to obtain a reference
            file.
        observatory : str
            telescope name used with CRDS,  e.g. 'jwst'.

        disable: bool or None
            Do not retrieve parameters from CRDS. If None, check global settings.

        Returns
        -------
        step_parameters : configobj
            The parameters as retrieved from CRDS. If there is an issue, log as such
            and return an empty config obj.
        """

        # Get the root logger, since the following operations
        # are happening in the surrounding architecture.
        logger = log.delegator.log
        pars_model = cls.get_pars_model()

        # If the dataset is not an operable DataModel, log as such and return
        # an empty config object
        try:
            model = dm_open(dataset)
        except (IOError, TypeError, ValueError):
            logger.warning('Input dataset is not a DataModel.')
            disable = True

        # Check if retrieval should be attempted.
        if disable is None:
            disable = get_disable_crds_steppars()
        if disable:
            logger.info(f'{pars_model.meta.reftype.upper()}: CRDS parameter reference retrieval disabled.')
            return config_parser.ConfigObj()

        # Retrieve step parameters from CRDS
        logger.debug(f'Retrieving step {pars_model.meta.reftype.upper()} parameters from CRDS')
        exceptions = crds_client.get_exceptions_module()
        try:
            ref_file = crds_client.get_reference_file(model,
                                                      pars_model.meta.reftype,
                                                      observatory=observatory,
                                                      asn_exptypes=['science'])
        except (AttributeError, exceptions.CrdsError, exceptions.CrdsLookupError):
            logger.debug(f'{pars_model.meta.reftype.upper()}: No parameters found')
            return config_parser.ConfigObj()
        if ref_file != 'N/A':
            logger.info(f'{pars_model.meta.reftype.upper()} parameters found: {ref_file}')
            ref = config_parser.load_config_file(ref_file)

            ref_pars = {
                par: value
                for par, value in ref.items()
                if par not in ['class', 'name']
            }
            logger.info(f'{pars_model.meta.reftype.upper()} parameters are {ref_pars}')

            return ref
        else:
            logger.debug(f'No {pars_model.meta.reftype.upper()} reference files found.')
            return config_parser.ConfigObj()

    @classmethod
    def reference_uri_to_cache_path(cls, reference_uri):
        """Convert an abstract CRDS reference URI to an absolute file path in the CRDS
        cache.  Reference URI's are typically output to dataset headers to record the
        reference files used.

        e.g. 'crds://jwst_miri_flat_0177.fits'  -->
            '/grp/crds/cache/references/jwst/jwst_miri_flat_0177.fits'

        The CRDS cache is typically located relative to env var CRDS_PATH
        with default value /grp/crds/cache.   See also https://jwst-crds.stsci.edu
        """
        return crds_client.reference_uri_to_cache_path(reference_uri)

    def set_primary_input(self, obj, exclusive=True):
        """
        Sets the name of the master input file and input directory.
        Used to generate output file names.

        Parameters
        ----------
        obj : str or DataModel
            The object to base the name on. If a datamodel,
            use Datamodel.meta.filename.

        exclusive : bool
            If True, only set if an input name is not already used
            by a parent Step. Otherwise, always set.
        """
        self._set_input_dir(obj, exclusive=exclusive)

        err_message = (
            'Cannot set master input file name from object'
            ' {}'.format(obj)
        )
        parent_input_filename = self.search_attr('_input_filename')
        if not exclusive or parent_input_filename is None:
            if isinstance(obj, str):
                self._input_filename = obj
            elif isinstance(obj, DataModel):
                try:
                    self._input_filename = obj.meta.filename
                except AttributeError:
                    self.log.debug(err_message)
            else:
                self.log.debug(err_message)

    def save_model(self,
                   model,
                   suffix=None,
                   idx=None,
                   output_file=None,
                   force=False,
                   format=None,
                   **components):
        """
        Saves the given model using the step/pipeline's naming scheme

        Parameters
        ----------
        model : jwst.datamodels.Model instance
            The model to save.

        suffix : str
            The suffix to add to the filename.

        idx : object
            Index identifier.

        output_file : str
            Use this file name instead of what the Step
            default would be.

        force : bool
            Regardless of whether `save_results` is `False`
            and no `output_file` is specified, try saving.

        format : str
            The format of the file name.  This is a format
            string that defines where `suffix` and the other
            components go in the file name. If False,
            it will be presumed `output_file` will have
            all the necessary formatting.

        components : dict
            Other components to add to the file name.

        Returns
        -------
        output_paths : [str[, ...]]
            List of output file paths the model(s) were saved in.
        """
        if output_file is None or output_file == '':
            output_file = self.output_file

        # Check if saving is even specified.
        if not force and \
           not self.save_results and \
           not output_file:
            return

        if isinstance(model, ModelContainer):
            save_model_func = partial(
                self.save_model,
                suffix=suffix,
                force=force,
                format=format,
                **components
            )
            output_path = model.save(
                path=output_file,
                save_model_func=save_model_func)
        else:

            # Search for an output file name.
            if (
                    self.output_use_model or
                    (output_file is None and not self.search_output_file)
            ):
                output_file = model.meta.filename
                idx = None
            output_path = model.save(
                self.make_output_path(
                    basepath=output_file,
                    suffix=suffix,
                    idx=idx,
                    name_format=format,
                    **components
                )
            )
            self.log.info('Saved model in {}'.format(output_path))

        return output_path

    @property
    def make_output_path(self):
        """Return function that creates the output path"""
        make_output_path = self.search_attr(
            '_make_output_path'
        )
        return partial(make_output_path, self)

    @staticmethod
    def _make_output_path(
            step,
            basepath=None,
            ext=None,
            suffix=None,
            name_format=None,
            component_format='',
            separator='_',
            **components
    ):
        """Create the output path

        Parameters
        ----------
        step : Step
            The `Step` in question.

        basepath : str or None
            The basepath to use. If None, `output_file`
            is used. Only the basename component of the path
            is used.

        ext : str or None
            The extension to use. If none, `output_ext` is used.
            Can include the leading period or not.

        suffix : str or None or False
            Suffix to append to the filename.
            If None, the `Step` default will be used.
            If False, no suffix replacement will be done.

        name_format : str or None
            The format string to use to form the base name.
            If False, it will be presumed that `basepath`
            has all the necessary formatting.

        component_format : str
            Format to use for the components

        separator : str
            Separator to use between replacement components

        components : dict
            dict of string replacements.

        Returns
        -------
        The fully qualified path name.

        Notes
        -----
        The values found in the `components` dict are placed in the string
        where the "{components}" replacement field is specified. If there are
        more than one component, the components are separated by the `separator`
        string.
        """
        if basepath is None and step.search_output_file:
            basepath = step.search_attr('output_file')
        if basepath is None:
            basepath = step.default_output_file()

        basename, basepath_ext = splitext(split(basepath)[1])
        if ext is None:
            ext = step.output_ext
        if ext is None and len(basepath_ext):
            ext = basepath_ext
        if ext.startswith('.'):
            ext = ext[1:]

        # Suffix check. An explicit check on `False` is necessary
        # because `None` is also allowed.
        suffix = _get_suffix(suffix, step=step)
        if suffix is not False:
            default_name_format = '{basename}{components}{suffix_sep}{suffix}.{ext}'
            suffix_sep = None
            if suffix is not None:
                basename, suffix_sep = remove_suffix(basename)
            if suffix_sep is None:
                suffix_sep = separator
        else:
            default_name_format = '{basename}{components}.{ext}'
            suffix = None
            suffix_sep = None

        # Setup formatting
        if name_format is None:
            name_format = default_name_format
        elif not name_format:
            name_format = basename + '.{ext}'
            basename = ''
            suffix_sep = ''
            separator = ''
        formatter = FormatTemplate(
            separator=separator,
            remove_unused=True
        )

        if len(components):
            component_str = formatter(component_format, **components)
        else:
            component_str = ''

        basename = formatter(
            name_format,
            basename=basename,
            suffix=suffix,
            suffix_sep=suffix_sep,
            ext=ext,
            components=component_str
        )

        output_dir = step.search_attr('output_dir', default='')
        output_dir = expandvars(expanduser(output_dir))
        full_output_path = join(output_dir, basename)

        return full_output_path

    def closeout(self, to_close=None, to_del=None):
        """Close out step processing

        Parameters
        ----------
        to_close : [object(, ...)]
            List of objects with a `close` method to execute
            The objects will also be deleted

        to_del : [object(, ...)]
            List of objects to simply delete

        Notes
        -----
        Other operations, such as forced garbage collection
        will also be done.
        """
        if to_close is None:
            to_close = []
        if to_del is None:
            to_del = []
        to_del += to_close
        for item in to_close:
            try:
                if hasattr(item, 'close'):
                    item.close()
            except Exception as exception:
                self.log.debug(
                    'Could not close "{}"'
                    'Reason:\n{}'.format(item, exception)
                )
        for item in to_del:
            try:
                del item
            except NameError as error:
                self.log.debug("An error has occurred: %s", error)
        gc.collect()

    def open_model(self, obj):
        """Open a datamodel

        Primarily a wrapper around `DataModel.open` to
        handle `Step` peculiarities

        Parameters
        ----------
        obj : object
            The object to open

        Returns
        -------
        datamodel : DataModel
            Object opened as a datamodel
        """
        return dm_open(self.make_input_path(obj))

    def make_input_path(self, file_path):
        """Create an input path for a given file path

        If `file_path` has no directory path, use `self.input_dir`
        as the directory path.

        Parameters
        ----------
        file_path : str or obj
            The supplied file path to check and modify.
            If anything other than `str`, the object
            is simply passed back.

        Returns
        -------
        full_path : str or obj
            File path using `input_dir` if the input
            had no directory path.
        """
        full_path = file_path
        if isinstance(file_path, str):
            original_path, file_name = split(file_path)
            if not len(original_path):
                full_path = join(self.input_dir, file_name)

        return full_path

    def load_as_level2_asn(self, obj):
        """Load object as an association

        Loads the specified object into a Level2 association.
        If necessary, prepend `Step.input_dir` to all members.

        Parameters
        ----------
        obj : object
            Object to load as a Level2 association

        Returns
        -------
        association : jwst.associations.lib.rules_level2_base.DMSLevel2bBase
            Association
        """
        asn = LoadAsLevel2Asn.load(obj, basename=self.output_file)
        update_key_value(asn, 'expname', (), mod_func=self.make_input_path)
        return asn

    def load_as_level3_asn(self, obj):
        """Load object as an association

        Loads the specified object into a Level3 association.
        If necessary, prepend `Step.input_dir` to all members.

        Parameters
        ----------
        obj : object
            Object to load as a Level3 association

        Returns
        -------
        association : jwst.associations.lib.rules_level3_base.DMS_Level3_Base
            Association
        """
        asn = LoadAsAssociation.load(obj)
        update_key_value(asn, 'expname', (), mod_func=self.make_input_path)
        return asn

    def _set_input_dir(self, input, exclusive=True):
        """Set the input directory

        If sufficient information is at hand, set a value
        for the attribute `input_dir`.

        Parameters
        ----------
        input : str
            Input to determine path from.

        exclusive : bool
            If True, only set if an input directory is not already
            defined by a parent Step. Otherwise, always set.

        """
        if not exclusive or self.search_attr('_input_dir') is None:
            try:
                if isfile(input):
                    self.input_dir = split(input)[0]
            except Exception:
                # Not a file-checkable object. Ignore.
                pass

    def record_step_status(self, datamodel, cal_step, success=True):
        """Record whether or not a step completed in meta.cal_step

        Parameters
        ----------
        datamodel : `~jwst.datamodels.Datamodel` instance
            This is the datamodel or container of datamodels to modify in place

        cal_step : str
            The attribute in meta.cal_step for recording the status of the step

        success : bool
            If True, then 'COMPLETE' is recorded.  If False, then 'SKIPPED'
        """
        if success:
            status = 'COMPLETE'
        else:
            status = 'SKIPPED'
            self.skip = True

        if isinstance(datamodel, ModelContainer):
            for model in datamodel:
                model.meta.cal_step._instance[cal_step] = status
        else:
            datamodel.meta.cal_step._instance[cal_step] = status

        # TODO: standardize cal_step naming to point to the offical step name

    @ClassInstanceMethod
    def get_pars(step, full_spec=True):
        """Retrieve the configuration parameters of a step

        Parameters
        ----------
        step : `Step`-derived class or instance
            The class or instance to retrieve the parameters for.

        full_spec : bool
            Return all parameters, including parent-specified parameters.
            If `False`, return only parameters specific to the class/instance.

        Returns
        -------
        pars : dict
            Keys are the parameters and values are the values.
        """
        from . import cmdline

        if full_spec:
            spec_file_func = config_parser.get_merged_spec_file
        else:
            spec_file_func = config_parser.load_spec_file
        spec = spec_file_func(step)
        if spec is None:
            return {}
        instance_pars = {}
        for key in spec:
            if hasattr(step, key):
                value = getattr(step, key)
                if not isinstance(value, property):
                    instance_pars[key] = value
        pars = config_parser.config_from_dict(instance_pars, spec, allow_missing=True)

        # Convert the config to a pure dict.
        pars_dict = {}
        for key, value in pars.items():
            if isinstance(value, cmdline.FromCommandLine):
                pars_dict[key] = str(value)
            else:
                pars_dict[key] = value
        return pars_dict

    @ClassInstanceMethod
    def get_pars_model(step, full_spec=True):
        """Return Step parameters as StepParsModel

        Parameters
        ----------
        step : `Step`-derived class or instance
            The `Step` or `Step` instance to retrieve the parameters model for.

        full_spec : bool
            Return all parameters, including parent-specified parameters.
            If `False`, return only parameters specific to the class/instance.

        Returns
        -------
        model : `StepParsModel`
            The `StepParsModel`.
        """
        pars_model = StepParsModel()
        pars_model.parameters.instance.update(step.get_pars(full_spec=full_spec))

        # Update class and name.
        full_class_name = _full_class_name(step)
        pars_model.parameters.instance.update({
            'class': full_class_name,
            'name': getattr(step, 'name', full_class_name.split('.')[-1])
        })
        pars_model.meta.reftype = 'pars-' + pars_model.parameters.name.lower()

        return pars_model

    def update_pars(self, parameters):
        """Update step parameters

        Only existing parameters are updated. Otherwise, new keys
        found in `parameters` are ignored.

        Parameters
        ----------
        parameters : dict
            Parameters to update.

        Notes
        -----
        `parameters` is presumed to have been produced by the
        `Step.get_pars` method. As such, the "steps" key is treated
        special in that it is a dict whose keys are the steps assigned
        directly as parameters to the current step. This is standard
        practice for `Pipeline`-based steps.
        """
        existing = self.get_pars().keys()
        for parameter, value in parameters.items():
            if parameter in existing:
                if parameter != 'steps':
                    setattr(self, parameter, value)
                else:
                    for step_name, step_parameters in value.items():
                        getattr(self, step_name).update_pars(step_parameters)
            else:
                self.log.debug(f'Parameter {parameter} is not valid for step {self}. Ignoring.')


# #########
# Utilities
# #########
def _full_class_name(obj):
    """Return the fully qualified class name

    Parameters
    ----------
    obj : object
        The object in question. Can be a class

    Returns
    class_name : str
        The full name
    """
    cls = obj if inspect.isclass(obj) else obj.__class__
    module = cls.__module__
    if module is None or module == str.__class__.__module__:
        return cls.__name__  # Avoid reporting __builtin__
    else:
        return module + '.' + cls.__name__


def _get_suffix(suffix, step=None, default_suffix=None):
    """Retrieve either specified or pipeline-supplied suffix

    Parameters
    ----------
    suffix : str or None
        Suffix to use if specified.

    step : Step or None
        The step to retrieve the suffux.

    default_suffix : str
        If the pipeline does not supply a suffix,
        use this.

    Returns
    -------
    suffix : str or None
        Suffix to use
    """
    if suffix is None and step is not None:
        suffix = step.search_attr('suffix')
    if suffix is None:
        suffix = default_suffix
    if suffix is None and step is not None:
        suffix = step.name.lower()
    return suffix


def get_disable_crds_steppars(default=None):
    """Return either the explicit default flag or retrieve from the environment

    If a default is not specified, retrieve the value from the environmental variable
    `STPIPE_DISABLE_CRDS_STEPPARS`.

    Parameters
    ----------
    default: str, bool, or None
        Flag to use. If None, the environmental is used.

    Returns
    -------
    flag: bool
        True to disable CRDS STEPPARS retrieval.
    """
    truths =  ('true', 'True', 't', 'yes', 'y')
    if default:
        if isinstance(default, bool):
            return default
        elif isinstance(default, str):
            return default in truths
        raise ValueError(f'default must be string or boolean: {default}')

    flag = os.environ.get('STPIPE_DISABLE_CRDS_STEPPARS', '')
    return flag in truths


@contextmanager
def preserve_step_pars(step):
    """Context manager to preserve step parameters

    Ensure step parameters are not modified during a block
    of operations. Allows local re-use of a Step instance without
    having to worry about side-effects on that Step. If used with
    a `Pipeline`, all substep parameters are also restored.

    Yields
    ------
    saved_pars: dict
        The saved parameters.
    """
    saved_pars = step.get_pars()
    try:
        yield saved_pars
    finally:
        step.update_pars(saved_pars)
