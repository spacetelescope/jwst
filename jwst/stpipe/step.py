"""
Step
"""
from __future__ import absolute_import, division, print_function

from astropy.extern import six
import contextlib
import gc
from os.path import dirname, join, basename, splitext, abspath, split
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

REMOVE_SUFFIX = '(.+?)(_(rate|cal)(ints)?)?$'


class Step(object):
    """
    Step
    """
    spec = """
    pre_hooks = string_list(default=list())
    post_hooks = string_list(default=list())

    output_dir = string(default=None)       # Directory path for output files
    output_file = output_file(default=None) # File to save output to.
    skip = boolean(default=False)           # Skip this step
    save_results = boolean(default=False)   # Force save results
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
                    if "process() takes exactly" in str(e):
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
                    if isinstance(result, datamodels.DataModel):
                        for ref_name, filename in self._reference_files_used:
                            if hasattr(result.meta.ref_file, ref_name):
                                getattr(result.meta.ref_file, ref_name).name = filename
                        result.meta.ref_file.crds.sw_version = crds_client.get_svn_version()
                        result.meta.ref_file.crds.context_used = crds_client.get_context_used()
                self._reference_files_used = []

            for post_hook in self._post_hooks:
                post_hook.run(result)

            # Mark versions
            for result in results:
                if isinstance(result, datamodels.DataModel):
                    result.meta.calibration_software_revision = __version_commit__
                    result.meta.calibration_software_version = __version__

            # Save the output file if one was specified
            if self.save_results or self.output_file is not None:
                result_id = _make_result_id(
                    self.output_file, len(results), self.name
                )
                make_output_path = self.search_attr(
                    'make_output_path', parent_first=True
                )
                for idx, result in enumerate(results):
                    if hasattr(result, 'save'):
                        output_path = make_output_path(
                            self, result,
                            basepath=self.output_file,
                            result_id=result_id(idx)
                        )
                        self.log.info('Saving file {0}'.format(output_path))
                        result.save(output_path, overwrite=True)

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

    def search_attr(self, attribute, parent_first=False):
        """Return first non-None attribute in step heirarchy

        Parameters
        ----------
        attribute: str
            The attribute to retrieve

        parent_first: bool
            If `True`, allow parent definition to override step version

        Returns
        -------
        value: obj
            Attribute value or None if not found
        """
        if parent_first:
            try:
                value = self.parent.search_attr(
                    attribute, parent_first=parent_first
                )
            except AttributeError:
                value = None
            if value is None:
                value = getattr(self, attribute, None)
            return value
        else:
            value = getattr(self, attribute, None)
            if value is None:
                try:
                    value = self.parent.search_attr(attribute)
                except AttributeError:
                    pass
            return value

    @classmethod
    def _is_association_file(cls, input_file):
        """Return True IFF `input_file` is an association file."""
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
        """
        gc.collect()
        if self.skip:
            return
        if self._is_association_file(input_file):
            return
        if len(self.reference_file_types):
            from .. import datamodels
            try:
                model = datamodels.open(input_file)
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

    @contextlib.contextmanager
    def get_reference_file_model(self, input_file, reference_file_type):
        """
        Get a reference file from CRDS as a jwst.datamodels.ModelBase
        object.  If the configuration file or commandline parameters
        override the reference file, it will be automatically used
        when calling this function.

        Parameters
        ----------
        input_file : jwst.datamodels.ModelBase instance
            A model of the input file.  Metadata on this input file
            will be used by the CRDS "bestref" algorithm to obtain a
            reference file.

        reference_file_type : string
            The type of reference file to retrieve.  For example, to
            retrieve a flat field reference file, this would be
            'flat_field'.

        Returns
        -------
        reference_file_model : jwst.datamodels.ModelBase instance
            A model to access the contents of the reference file.
        """
        from .. import datamodels

        filename = self.get_reference_file(input, reference_file_type)
        with datamodels.open(filename) as model:
            yield model
        gc.collect()

    def set_input_filename(self, path):
        """
        Sets the name of the master input file.  Used to generate output
        file names.
        """
        self._input_filename = path

    def save_model(self, model, suffix, *args, **kwargs):
        """
        Saves the given model using the step/pipeline's naming scheme

        Parameters
        ----------
        model : jwst.datamodels.Model instance
            The model to save.

        suffix : str
            The suffix to add to the filename.

        Notes
        -----
        This routine is used to save data outside of the normal step
        cycle, where results are saved at the end of a step.

        Serious consideration should be given to the step design
        if this call is needed. In particular, if such data
        is important to save, consider making the producing code
        its own step, such that it is subject to the output controls
        of the step infrastructure.
        """

        # Get the output path as defined by the current step.
        make_output_path = self.search_attr(
            'make_output_path', parent_first=True
        )
        output_path = make_output_path(
            self, model, suffix=suffix, ignore_use_model=True
        )

        model.save(output_path, *args, **kwargs)

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
            suffix = _get_suffix(suffix, default_suffix=result_id)
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
                    suffix = _get_suffix(
                        suffix, step=step, default_suffix=result_id
                    )
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
    return suffix
