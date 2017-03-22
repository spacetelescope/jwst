from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os.path as op
import os
import copy
import warnings
from collections import OrderedDict

from asdf import AsdfFile
from astropy.extern import six

from ..associations import (
    AssociationError,
    AssociationNotValidError, load_asn)
from . import model_base
from .util import open as datamodel_open


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
    >>> container = datamodels.ModelContainer('example_asn.json')
    >>> for dm in container:
    ...     print(dm.meta.filename)

    Say the association was a NIRCam dithered dataset. The `models_grouped`
    attribute is a list of lists, the first index giving the list of exposure
    groups, with the second giving the individual datamodels representing
    each detector in the exposure (2 or 8 in the case of NIRCam).

    >>> total_exposure_time = 0.0
    >>> for group in container.models_grouped:
    ...     total_exposure_time += group[0].meta.exposure.exposure_time

    >>> c = datamodels.ModelContainer()
    >>> m = datamodels.open('myfile.fits')
    >>> c.append(m)
    """

    # This schema merely extends the 'meta' part of the datamodel, and
    # does not describe the data contents of the container.
    schema_url = "container.schema.yaml"

    def __init__(self, init=None, **kwargs):

        super(ModelContainer, self).__init__(init=None, **kwargs)

        if init is None:
            self._models = []
        elif isinstance(init, list):
            for item in init:
                if not isinstance(item, model_base.DataModel):
                    raise ValueError('list must contain only DataModels')
            self._models = init
        elif isinstance(init, self.__class__):
            instance = copy.deepcopy(init._instance)
            self._schema = init._schema
            self._shape = init._shape
            self._asdf = AsdfFile(instance, extensions=self._extensions)
            self._instance = instance
            self._ctx = self
            self.__class__ = init.__class__
            self._models = init._models
        elif isinstance(init, six.string_types):
            try:
                self.from_asn(init, **kwargs)
            except (IOError):
                raise IOError('Cannot open files.')
            except AssociationError:
                raise AssociationError('{0} must be an ASN file'.format(init))
        else:
            raise TypeError('Input {0!r} is not a list of DataModels or '
                            'an ASN file'.format(init))

        self.__assign_group_ids()


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


    def pop(self, index=None):
        if not index:
            self._models.pop(-1)
        else:
            self._models.pop(index)


    def copy(self):
        """
        Returns a deep copy of the models in this model container.
        """

        models_copy = [m.copy() for m in self._models]
        return self.__class__(init=models_copy)


    def from_asn(self, filepath, **kwargs):
        """
        Load fits files from a JWST association file.

        Parameters
        ----------
        filepath : str
            The path to an association file.
        """

        filepath = op.abspath(op.expanduser(op.expandvars(filepath)))
        basedir = op.dirname(filepath)
        try:
            with open(filepath) as asn_file:
                asn_data = load_asn(asn_file)
        except AssociationNotValidError:
            raise IOError("Cannot read ASN file.")

        # make a list of all the input FITS files
        infiles = [op.join(basedir, member['expname']) for member
                   in asn_data['products'][0]['members']]
        try:
            self._models = [datamodel_open(infile, **kwargs) for infile in infiles]
        except IOError:
            raise IOError('Cannot open {}'.format(infiles))

        # Pull the whole association table into meta.asn_table
        self.meta.asn_table = {}
        model_base.properties.merge_tree(self.meta.asn_table._instance, asn_data)

        # populate the output metadata with the output file from the ASN file
        # Should remove the following lines eventually
        self.meta.resample.output = str(asn_data['products'][0]['name'])
        self.meta.table_name = str(filepath)
        self.meta.pool_name = str(asn_data['asn_pool'])
        self.meta.targname = str(asn_data['target'])
        self.meta.program = str(asn_data['program'])
        self.meta.asn_type = str(asn_data['asn_type'])
        self.meta.asn_rule = str(asn_data['asn_rule'])


    def save(self, path=None, *args, **kwargs):
        """
        Write out models in container to FITS or ASDF.

        Parameters
        ----------

        path : string
            Directory to write out files.  Defaults to current working dir.
            If directory does not exist, it creates it.  Filenames are pulled
            from `.meta.filename` of each datamodel in the container.
        """
        if path is None:
            path = os.getcwd()
        try:
            for model in self._models:
                outpath = op.join(path, model.meta.filename)
                model.save(outpath, *args, **kwargs)
        except IOError as err:
            raise err


    def __assign_group_ids(self):
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
        meta.instrument.name
        meta.instrument.channel
        """
        for model in self._models:
            try:
                model_attrs = []
                model_attrs.append(model.meta.observation.program_number)
                model_attrs.append(model.meta.observation.observation_number)
                model_attrs.append(model.meta.observation.visit_number)
                model_attrs.append(model.meta.observation.visit_group)
                model_attrs.append(model.meta.observation.sequence_id)
                model_attrs.append(model.meta.observation.activity_id)
                model_attrs.append(model.meta.observation.exposure_number)
                model_attrs.append(model.meta.instrument.name)
                model_attrs.append(model.meta.instrument.channel)
                group_id = ('jw' + '_'.join([
                                ''.join(model_attrs[:3]),
                                ''.join(model_attrs[3:6]),
                                model_attrs[6], model_attrs[7].lower(),
                                model_attrs[8].lower()]))
                model.meta.group_id = group_id
            except:
                w = '`{}` is missing'.format(model.meta.filename) + \
                    ' metadata. Grouping by exposure may not be correct.'
                warnings.warn(w, RuntimeWarning, stacklevel=2)
                model.meta.group_id = None


    @property
    def models_grouped(self):
        """
        Returns a list of a list of datamodels grouped by exposure.
        """
        self.__assign_group_ids()
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
