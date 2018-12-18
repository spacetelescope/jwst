import copy
from collections import OrderedDict
import os.path as op

from asdf import AsdfFile

from ..associations import (
    Association, AssociationError,
    AssociationNotValidError, load_asn)

from . import model_base
from .util import open as datamodel_open
from .util import is_association

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__doctest_skip__ = ['ModelContainer']

__all__ = ['ModelContainer']


class ModelContainer(model_base.DataModel):
    """
    A container for holding DataModels.

    This functions like a list for holding DataModel objects.  It can be
    iterated through like a list, DataModels within the container can be
    addressed by index, and the datamodels can be grouped into a list of
    lists for grouped looping, useful for NIRCam where grouping together
    all detectors of a given exposure is useful for some pipeline steps.

    Parameters
    ----------
    init : file path, list of DataModels, or None

        - file path: initialize from an association table

        - list: a list of DataModels of any type

        - None: initializes an empty `ModelContainer` instance, to which
          DataModels can be added via the ``append()`` method.

    persist : bool
        If True, do not close model after opening it.


    Examples
    --------
    >>> container = ModelContainer('example_asn.json')
    >>> for dm in container:
    ...     print(dm.meta.filename)

    Say the association was a NIRCam dithered dataset. The `models_grouped`
    attribute is a list of lists, the first index giving the list of exposure
    groups, with the second giving the individual datamodels representing
    each detector in the exposure (2 or 8 in the case of NIRCam).

    >>> total_exposure_time = 0.0
    >>> for group in container.models_grouped:
    ...     total_exposure_time += group[0].meta.exposure.exposure_time

    >>> c = ModelContainer()
    >>> m = datamodels.open('myfile.fits')
    >>> c.append(m)
    """

    # This schema merely extends the 'meta' part of the datamodel, and
    # does not describe the data contents of the container.
    schema_url = "container.schema.yaml"

    def __init__(self, init=None, persist=True, **kwargs):

        super(ModelContainer, self).__init__(init=None, **kwargs)
        self._persist = persist

        if init is None:
            self._models = []
        elif isinstance(init, list):
            self._validate_model(init)
            self._models = init[:]
        elif isinstance(init, self.__class__):
            instance = copy.deepcopy(init._instance)
            self._schema = init._schema
            self._shape = init._shape
            self._asdf = AsdfFile(instance, extensions=init._extensions)
            self._instance = instance
            self._ctx = self
            self.__class__ = init.__class__
            self._models = init._models
        elif is_association(init):
            self.from_asn(init)
        elif isinstance(init, str):
            try:
                init_from_asn = self.read_asn(init)
            except (IOError):
                raise IOError('Cannot open files.')
            except AssociationError:
                raise AssociationError('{0} must be an ASN file'.format(init))
            self.from_asn(init_from_asn, asn_file_path=init)
        else:
            raise TypeError('Input {0!r} is not a list of DataModels or '
                            'an ASN file'.format(init))

    def _open_model(self, index):
        model = self._models[index]
        if isinstance(model, str):
            model = datamodel_open(
                model,
                extensions=self._extensions,
                pass_invalid_values=self._pass_invalid_values
            )
            self._models[index] = model

        return model

    def _close_model(self, filename, index):
        if not self._persist:
            self._models[index].close()
            self._models[index] = filename

    def _validate_model(self, models):
        if not isinstance(models, list):
            models = [models]
        for model in models:
            if isinstance(model, ModelContainer):
                raise ValueError(
                    "ModelContainer cannot contain ModelContainer"
                )
            if not isinstance(model, (str, model_base.DataModel)):
                raise ValueError('model must be string or DataModel')

    def __len__(self):
        return len(self._models)

    def __getitem__(self, index):
        return self._open_model(index)

    def __setitem__(self, index, model):
        self._validate_model(model)
        self._models[index] = model

    def __delitem__(self, index):
        del self._models[index]

    def __iter__(self):
        return ModelContainerIterator(self)

    def insert(self, index, model):
        self._validate_model(model)
        self._models.insert(index, model)

    def append(self, model):
        self._validate_model(model)
        self._models.append(model)

    def extend(self, models):
        self._validate_model(models)
        self._models.extend(models)

    def pop(self, index=-1):
        self._open_model(index)
        return self._models.pop(index)

    def copy(self, memo=None):
        """
        Returns a deep copy of the models in this model container.
        """
        result = self.__class__(init=None,
                                extensions=self._extensions,
                                pass_invalid_values=self._pass_invalid_values,
                                strict_validation=self._strict_validation)
        instance = copy.deepcopy(self._instance, memo=memo)
        result._asdf = AsdfFile(instance, extensions=self._extensions)
        result._instance = instance
        result._iscopy = self._iscopy
        result._schema = result._schema
        result._ctx = result
        for m in self._models:
            if isinstance(m, model_base.DataModel):
                result.append(m.copy())
            else:
                result.append(m)
        return result

    def read_asn(self, filepath):
        """
        Load fits files from a JWST association file.

        Parameters
        ----------
        filepath : str
            The path to an association file.
        """

        filepath = op.abspath(op.expanduser(op.expandvars(filepath)))
        try:
            with open(filepath) as asn_file:
                asn_data = load_asn(asn_file)
        except AssociationNotValidError:
            raise IOError("Cannot read ASN file.")
        return asn_data

    def from_asn(self, asn_data, asn_file_path=None):
        """
        Load fits files from a JWST association file.

        Parameters
        ----------
        asn_data : Association
            An association dictionary

        asn_file_path: str
            Filepath of the association, if known.
        """
        # make a list of all the input files
        infiles = [member['expname'] for member
                   in asn_data['products'][0]['members']]
        self._models = infiles

        # Pull the whole association table into meta.asn_table
        self.meta.asn_table = {}
        model_base.properties.merge_tree(
            self.meta.asn_table._instance, asn_data
        )

        self.meta.resample.output = asn_data['products'][0]['name']
        if asn_file_path is None:
            self.meta.table_name = 'not specified'
        else:
            self.meta.table_name = op.basename(asn_file_path)
        self.meta.pool_name = asn_data['asn_pool']

    def save(self,
             path=None,
             dir_path=None,
             save_model_func=None,
             *args, **kwargs):
        """
        Write out models in container to FITS or ASDF.

        Parameters
        ----------
        path : str or func or None
            - If None, the `meta.filename` is used for each model.
            - If a string, the string is used as a root and an index is
              appended.
            - If a function, the function takes the two arguments:
              the value of model.meta.filename and the
              `idx` index, returning constructed file name.

        dir_path : str
            Directory to write out files.  Defaults to current working dir.
            If directory does not exist, it creates it.  Filenames are pulled
            from `.meta.filename` of each datamodel in the container.

        save_model_func: func or None
            Alternate function to save each model instead of
            the models `save` method. Takes one argument, the model,
            and keyword argument `idx` for an index.

        Returns
        -------
        output_paths: [str[, ...]]
            List of output file paths of where the models were saved.
        """
        output_paths = []
        if path is None:
            path = lambda filename, idx: filename
        elif not callable(path):
            path = make_file_with_index

        for idx, model in enumerate(self):
            if len(self) <= 1:
                idx = None
            if save_model_func is None:
                outpath, filename = op.split(
                    path(model.meta.filename, idx=idx)
                )
                if dir_path:
                    outpath = dir_path
                save_path = op.join(outpath, filename)
                try:
                    output_paths.append(
                        model.save(save_path, *args, **kwargs)
                    )
                except IOError as err:
                    raise err

            else:
                output_paths.append(save_model_func(model, idx=idx))

        return output_paths

    @property
    def models_grouped(self):
        """
        Returns a list of a list of datamodels grouped by exposure.

        Data from different detectors of the same exposure will have the
        same group id, which allows grouping by exposure.  The following
        metadata is used for grouping:

        meta.observation.program_number
        meta.observation.observation_number
        meta.observation.visit_number
        meta.observation.visit_group
        meta.observation.sequence_id
        meta.observation.activity_id
        meta.observation.exposure_number
        """
        unique_exposure_parameters = [
            'program_number',
            'observation_number',
            'visit_number',
            'visit_group',
            'sequence_id',
            'activity_id',
            'exposure_number'
            ]
        group_dict = OrderedDict()

        for i in range(len(self)):
            model = self._open_model(i)
            params = []
            for param in unique_exposure_parameters:
                params.append(getattr(model.meta.observation, param))
            try:
                group_id = ('jw' + '_'.join([''.join(params[:3]),
                        ''.join(params[3:6]), params[6]]))
                model.meta.group_id = group_id
            except TypeError:
                params_dict = dict(zip(unique_exposure_parameters, params))
                bad_params = {'meta.observation.'+k:v for k, v in params_dict.items() if not v}
                log.warning(
                    'Cannot determine grouping of exposures: '
                    '{}'.format(bad_params)
                    )
                model.meta.group_id = 'exposure{0:04d}'.format(i + 1)

            group_id = model.meta.group_id
            if group_id in group_dict:
                group_dict[group_id].append(model)
            else:
                group_dict[group_id] = [model]

        return group_dict.values()

    @property
    def group_names(self):
        """
        Return list of names for the DataModel groups by exposure.
        """
        result = []
        for group in self.models_grouped:
            result.append(group[0].meta.group_id)
        return result

    def __get_recursively(self, field, search_dict):
        """
        Takes a dict with nested lists and dicts, and searches all dicts for
        a key of the field provided.
        """
        values_found = []
        for key, value in search_dict.items():
            if key == field:
                values_found.append(value)
            elif isinstance(value, dict):
                results = self.__get_recursively(field, value)
                for result in results:
                    values_found.append(result)
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        more_results = self.__get_recursively(field, item)
                        for another_result in more_results:
                            values_found.append(another_result)
        return values_found

    def get_recursively(self, field):
        """
        Returns a list of values of the specified field from meta.
        """
        return self.__get_recursively(field, self.meta._instance)


class ModelContainerIterator:
    """
    An iterator for model containers that opens one model at a time
    """
    def __init__(self, container):
        self.index = -1
        self.open_filename = None
        self.container = container

    def __iter__(self):
        return self

    def __next__(self):
        if self.open_filename is not None:
            self.container._close_model(self.open_filename, self.index)
            self.open_filename = None

        self.index += 1
        if self.index < len(self.container._models):
            model = self.container._models[self.index]
            if isinstance(model, str):
                name = model
                model = self.container._open_model(self.index)
                self.open_filename = name
            return model
        else:
            raise StopIteration


# #########
# Utilities
# #########
def make_file_with_index(file_path, idx):
    """Append an index to a filename

    Parameters
    ----------
    file_path: str
        The file to append the index to.
    idx: int
        An index to append


    Returns
    -------
    file_path: str
        Path with index appended
    """
    # Decompose path
    path_head, path_tail = op.split(file_path)
    base, ext = op.splitext(path_tail)
    if idx is not None:
        base = base + str(idx)
    return op.join(path_head, base + ext)
