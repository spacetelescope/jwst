import copy
from collections import OrderedDict
import os.path as op
import warnings

from asdf import AsdfFile
from astropy.io import fits

from ..associations import (
    AssociationNotValidError,
    load_asn)

from . import model_base
from .util import open as datamodel_open
from .util import is_association

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

    def __init__(self, init=None, **kwargs):

        super(ModelContainer, self).__init__(init=None, **kwargs)

        self._models = []

        if init is None:
            # Don't populate the container with models
            pass
        elif isinstance(init, fits.HDUList):
            self._models.append([datamodel_open(init)])
        elif isinstance(init, list):
            if all(isinstance(x, (str, fits.HDUList)) for x in init):
                # Try opening the list of files as datamodels
                try:
                    init = [datamodel_open(m) for m in init]
                except (FileNotFoundError, ValueError):
                    raise
            elif not all(isinstance(x, model_base.DataModel) for x in init):
                raise TypeError('list must contain DataModels')
            self._models = init
        elif isinstance(init, self.__class__):
            instance = copy.deepcopy(init._instance)
            self._schema = init._schema
            self._shape = init._shape
            self._asdf = AsdfFile(instance)
            self._instance = instance
            self._ctx = self
            self.__class__ = init.__class__
            self._models = init._models
        elif is_association(init):
            self.from_asn(init)
        elif isinstance(init, str):
            init_from_asn = self.read_asn(init)
            self.from_asn(init_from_asn, asn_file_path=init)
        else:
            raise TypeError('Input {0!r} is not a list of DataModels or '
                            'an ASN file'.format(init))

    def __len__(self):
        return len(self._models)

    def __getitem__(self, index):
        return self._models[index]

    def __setitem__(self, index, model):
        self._models[index] = model

    def __delitem__(self, index):
        del self._models[index]

    def __iter__(self):
        for model in self._models:
            yield model

    def insert(self, index, model):
        self._models.insert(index, model)

    def append(self, model):
        self._models.append(model)

    def extend(self, model):
        self._models.extend(model)

    def pop(self, index=-1):
        self._models.pop(index)

    def copy(self, memo=None):
        """
        Returns a deep copy of the models in this model container.
        """
        result = self.__class__(init=None,
                                pass_invalid_values=self._pass_invalid_values,
                                strict_validation=self._strict_validation)
        instance = copy.deepcopy(self._instance, memo=memo)
        result._asdf = AsdfFile(instance)
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
        if asn_file_path:
            asn_dir = op.dirname(asn_file_path)
            infiles = [op.join(asn_dir, f) for f in infiles]
        try:
            self._models = [datamodel_open(infile) for infile in infiles]
        except IOError:
            raise IOError('Cannot open {}'.format(infiles))

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

    def _assign_group_ids(self):
        """
        Assign an ID grouping by exposure.

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

        for i, model in enumerate(self._models):
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
                warnings.warn(
                    'Cannot determine grouping of exposures: '
                    '{}'.format(bad_params)
                    )
                model.meta.group_id = 'exposure{0:04d}'.format(i + 1)

    @property
    def models_grouped(self):
        """
        Returns a list of a list of datamodels grouped by exposure.
        """
        self._assign_group_ids()
        group_dict = OrderedDict()
        for model in self._models:
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
