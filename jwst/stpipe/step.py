"""
Step
"""
from functools import partial
import gc
from os.path import (
    abspath,
    basename,
    dirname,
    expanduser,
    expandvars,
    join,
    split,
    splitext,
)
import re
import sys

try:
    from astropy.io import fits
    DISCOURAGED_TYPES = (fits.HDUList,)
except ImportError:
    DISCOURAGED_TYPES = None

from . import config_parser
from . import crds_client
from . import log
from . import utilities
from .. import __version_commit__, __version__
from ..datamodels import DataModel

SUFFIX_LIST = [
    'cal', 'calints', 'crf', 'crfints',
    'dark',
    'i2d',
    'jump',
    'psfalign', 'psfstack', 'psfsub',
    'ramp', 'rate', 'rateints',
    's2d', 's3d',
    'uncal',
    'x1d', 'x1dints',
]
REMOVE_SUFFIX = '^(?P<root>.+?)((?P<separator>_|-)(' \
                + '|'.join(SUFFIX_LIST) + '))?$'


class Step():
    """
    Step
    """
    spec = """
    pre_hooks        = string_list(default=list())
    post_hooks       = string_list(default=list())
    output_file      = output_file(default=None)   # File to save output to.
    output_dir       = string(default=None)        # Directory path for output files
    output_ext       = string(default='.fits')     # Default type of output
    output_use_model = boolean(default=False)      # When saving use `DataModel.meta.filename`
    output_use_index = boolean(default=True)       # Append index.
    save_results     = boolean(default=False)      # Force save results
    skip             = boolean(default=False)      # Skip this step
    suffix           = string(default=None)        # Default suffix for output files
    """

    # Reference types for both command line override definition and reference prefetch
    reference_file_types = []

    # Set to False in subclasses to skip prefetch,  but by default attempt to prefetch
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
            spec[override_name] = 'string(default=None)'
            spec.inline_comments[override_name] = (
                '# Override the {0} reference file'.format(
                    reference_file_type))
        return spec

    @classmethod
    def print_configspec(cls, stream=sys.stdout):

        # Python2/3 issue: Python3 doesn't like bytes
        # going to stdout directly.
        if stream == sys.stdout:
            try:
                stream = sys.stdout.buffer
            except AttributeError:
                pass

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

        try:
            if (len(args) and len(self.reference_file_types) and not self.skip
                and self.prefetch_references):
                self._precache_reference_files(args[0])

            self.log.info(
                'Step {0} running with args {1}.'.format(
                    self.name, args))

            if self.output_file is not None:
                self.save_results = True

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
                try:
                    step_result = self.process(*args)
                except TypeError as e:
                    if "process() takes exactly" in str(e):
                        raise TypeError("Incorrect number of arguments to step")
                    raise

            # Warn if returning a discouraged object
            self._check_args(step_result, DISCOURAGED_TYPES, "Returned")

            # Run the post hooks
            for post_hook in self._post_hooks:
                hook_results = post_hook.run(step_result)
                if hook_results is not None:
                    step_result = hook_results

            # Update meta information
            if not isinstance(step_result, (list, tuple, datamodels.ModelContainer)):
                results = [step_result]
            else:
                results = step_result

            if len(self._reference_files_used):
                for result in results:
                    if isinstance(result, datamodels.DataModel):
                        for ref_name, filename in self._reference_files_used:
                            if hasattr(result.meta.ref_file, ref_name):
                                getattr(result.meta.ref_file, ref_name).name = filename
                        result.meta.ref_file.crds.sw_version = crds_client.get_svn_version()
                        result.meta.ref_file.crds.context_used = crds_client.get_context_used()
                self._reference_files_used = []

            # Mark versions
            for result in results:
                if isinstance(result, datamodels.DataModel):
                    result.meta.calibration_software_revision = __version_commit__
                    result.meta.calibration_software_version = __version__

            # Save the output file if one was specified
            if not self.skip and (
                    self.save_results or self.output_file is not None
            ):
                # Setup the save list.
                if not isinstance(step_result, (list, tuple)):
                    results_to_save = [step_result]
                else:
                    results_to_save = step_result

                result_id = _make_result_id(
                    self.output_file, len(results_to_save), self.name
                )
                make_output_path = self.search_attr(
                    'make_output_path', parent_first=True
                )
                for idx, result in enumerate(results_to_save):
                    if len(results_to_save) <= 1:
                        idx = None
                    if isinstance(result, DataModel):
                        self.save_model(result, idx=idx)
                    elif hasattr(result, 'save'):
                        try:
                            out_dir, out_name = make_output_path(
                                self, result,
                                basepath=self.output_file,
                                result_id=result_id(idx)
                            )
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
                            out_path = join(out_dir, out_name)
                            self.log.info(
                                'Saving file {0}'.format(out_path)
                            )
                            result.save(out_path, overwrite=True)

            self.log.info(
                'Step {0} done'.format(self.name))
        finally:
            log.delegator.log = orig_log

        return step_result

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
        Creates and runs a new instance of the class.

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

    def search_attr(self, attribute, parent_first=False, default=None):
        """Return first non-None attribute in step heirarchy

        Parameters
        ----------
        attribute: str
            The attribute to retrieve

        parent_first: bool
            If `True`, allow parent definition to override step version

        default: obj
            If attribute is not found, the value to use

        Returns
        -------
        value: obj
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

    @classmethod
    def _is_container(cls, input_file):
        """Return True IFF `input_file` is a ModelContainer or successfully
        loads as an association.
        """
        from ..associations import load_asn
        from .. import datamodels
        if isinstance(input_file, datamodels.ModelContainer):
            return True
        try:
            with open(input_file, 'r') as input_file_fh:
                asn = load_asn(input_file_fh)
        except:
            return False
        return True

    def _precache_reference_files(self, input_file):
        """
        Precache all of the expected reference files before the Step's
        process method is called.

        Perform's garbage collection before and after prefetching.

        Handles opening `input_file` as a model if it is a filename.

        input_file:  filename, model container, or model

        returns:  None
        """
        gc.collect()
        from .. import datamodels
        try:
            with datamodels.open(input_file) as model:
                self._precache_reference_files_opened(model)
        except (ValueError, TypeError, IOError):
            self.log.info(
                'First argument {0} does not appear to be a '
                'model'.format(input_file))
        gc.collect()

    def _precache_reference_files_opened(self, model_or_container):
        """Pre-fetches references for `model_or_container`.

        Handles recursive pre-fetches for any models inside a container,
        or just a single model.

        Assumes model_or_container is an open model or container object,
        not a filename.

        No garbage collection.
        """
        if self._is_container(model_or_container):
            # recurse on each contained model
            for contained_model in model_or_container:
                self._precache_reference_files_opened(contained_model)
        else:
            # precache a single model object
            self._precache_reference_files_impl(model_or_container)

    def _precache_reference_files_impl(self, model):
        """Given open data `model`,  determine and cache reference files for
        any reference types which are not overridden on the command line.

        Verify that all CRDS and overridden reference files are readable.

        model:  An open Model object;  not a filename, ModelContainer, etc.
        """
        ovr_refs = {
            reftype: self._get_ref_override(reftype)
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

    def reference_uri_to_cache_path(self, reference_uri):
        """Convert an abstract CRDS reference URI to an absolute file path in the CRDS
        cache.  Reference URI's are typically output to dataset headers to record the
        reference files used.

        e.g. 'crds://jwst_miri_flat_0177.fits'  -->
            '/grp/crds/cache/references/jwst/jwst_miri_flat_0177.fits'

        The CRDS cache is typically located relative to env var CRDS_PATH
        with default value /grp/crds/cache.   See also https://jwst-crds.stsci.edu
        """
        return crds_client.reference_uri_to_cache_path(reference_uri)

    def set_input_filename(self, path):
        """
        Sets the name of the master input file.  Used to generate output
        file names.
        """
        self._input_filename = path

    def save_model(self, model, suffix=None, idx=None, output_file=None, force=False):
        """
        Saves the given model using the step/pipeline's naming scheme

        Parameters
        ----------
        model : jwst.datamodels.Model instance
            The model to save.

        suffix : str
            The suffix to add to the filename.

        idx: object
            Index identifier.

        output_file: str
            Use this file name instead of what the Step
            default would be.

        force: bool
            Regardless of whether `save_results` is `False`
            and no `output_file` is specified, try saving.
        """
        if output_file is None:
            output_file = self.output_file

        # Check if saving is even specified.
        if not force and \
           not self.save_results and \
           not output_file:
            return

        # Get the output path as defined by the current step.
        make_output_path_partial = partial(
            self.make_output_path,
            self,
            basepath=output_file,
            suffix=suffix,
            idx=idx
        )
        model.save(make_output_path_partial)

    @property
    def make_output_path(self):
        """Return function that creates the output path"""
        make_output_path = self.search_attr(
            '_make_output_path', parent_first=True
        )
        return make_output_path

    @staticmethod
    def _make_output_path(
            step, data,
            basepath=None, suffix=None, ext=None,
            idx=None, ignore_use_model=False
    ):
        """Make up a path based on data and user specification

        Parameters
        ----------
        step: Step
            The step which produced the data

        data: obj
            Get filename info from the data to be saved.

        basepath: str or None
            The output file name. If `None` or empty string, create
            a filename based on the data.

        suffix: str or None
            The suffix to append to the basename.

        ext: str or None
            The file format extension

        idx: object or None
            An index to use. Index is not used if `output_use_model`
            is True or `output_use_index` is False

        ignore_use_model: bool
            Ignore configuration parameter `output_use_model`

        Returns
        -------
        output_path, file_name, path_name: str, str, str
            2-tuple consisting of:
            - path_name: Directory path
            - file_name: File name
        """
        from ..datamodels import DataModel

        # Determine source of the basepath
        output_path = basepath
        has_basepath = not step.output_use_model and \
                       basepath is not None and \
                       len(basepath) > 0
        use_model_name = (not ignore_use_model) and \
                         getattr(step, 'output_use_model', False)

        # Determine if the index should be used
        if idx is not None and (
                use_model_name or not getattr(step, 'output_use_index', True)
        ):
            idx = None

        # Get the suffix
        suffix = _get_suffix(suffix, step=step)

        # If not given a name, construct one.
        if not has_basepath:

            # Make names based on DataModels
            if isinstance(data, DataModel):

                # If using what is in the model, just retrieve that.
                if use_model_name:
                    output_path = data.meta.filename

                # Otherwise, create a fully qualified pipeline outputname
                else:
                    output_path = step.search_attr('output_file')
                    if output_path is None:
                        output_path = data.meta.filename

        # Construct the output name
        path, filename = split(output_path)
        name, ext = splitext(filename)

        # Remove any known, previous suffixes.
        separator = None
        if suffix is not None:
            match = re.match(REMOVE_SUFFIX, name)
            try:
                name = match.group('root')
                separator = match.group('separator')
            except AttributeError:
                pass
        if separator is None:
            separator = '_'

        output_name = [name]
        if idx is not None:
            output_name.append(separator + str(idx))
        if suffix is not None:
            output_name.append(separator + suffix)
        if ext is None or len(ext) == 0:
            ext = step.search_attr('output_ext')
        output_name.append(ext)
        output_name = ''.join(output_name)

        output_dir = step.search_attr('output_dir', default='')
        output_dir = expandvars(expanduser(output_dir))
        step.log.info('Saving file {}'.format(
            join(output_dir, output_name)
        ))
        return output_dir, output_name

    def closeout(self, to_close=None, to_del=None):
        """Close out step processing

        Parameters
        ----------
        to_close: [object(, ...)]
            List of objects with a `close` method to execute
            The objects will also be deleted

        to_del: [object(, ...)]
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
                item.close()
            except Exception as exception:
                self.logger.debug(
                    'Could not close "{}"'
                    'Reason:\n{}'.format(item, exception)
                )
        for item in to_del:
            try:
                del item
            except Exception as exception:
                self.logger.debug(
                    'Could not delete "{}"'
                    'Reason:\n{}'.format(item, exception)
                )
        gc.collect()


# #########
# Utilities
# #########
def _make_result_id(output_file, n_results, default_name):
    """Create function the constructs a result identifier

    Parameters
    ----------
    output_file: str or None
        output file name

    n_results: int
        The number of results to be saved.

    default_name: str
        The name to use if `output_file` is not defined

    Returns
    -------
    result_id: func
        A function that takes an int and produces a str
        the represents a result identifier
    """
    id_format = []
    if output_file is None:
        id_format.append(default_name.lower())
    if n_results > 1:
        id_format.append('{idx}')
    if len(id_format):
        id_format = '_'.join(id_format)
        result_id = lambda idx: id_format.format(idx=str(idx))
    else:
        result_id = lambda idx: None
    return result_id


def _get_suffix(suffix, step=None, default_suffix=None):
    """Retrieve either specified or pipeline-supplied suffix

    Parameters
    ----------
    suffix: str or None
        Suffix to use if specified.

    step: Step or None
        The step to retrieve the suffux.

    default_suffix: str
        If the pipeline does not supply a suffix,
        use this.

    Returns
    -------
    suffix: str or None
        Suffix to use
    """
    if suffix is None and step is not None:
        suffix = step.search_attr('suffix')
    if suffix is None:
        suffix = default_suffix
    if suffix is None and step is not None:
        suffix = step.name.lower()
    return suffix
