from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os.path as op
import os
import copy
from collections import OrderedDict

from asdf import AsdfFile
from astropy.extern import six

from ..associations import Association
from . import model_base
from .util import open as datamodel_open


__all__ = ['ModelContainer']


class ModelContainer(model_base.DataModel):
    """
    A container for holding DataModels.

    Parameters
    ----------
    init : file path, list of DataModels, or None

        - file path: initialize from an association table

        - list: a list of any DataModel models

        - None: initializes an empty `ModelContainer` instance, to which
          DataModel models can added via the ``append()`` method.

    Examples
    --------
    >>> c = ModelContainer('example_asn.json')
    >>> c[0]    # the first DataModel in the container
    >>> c.models_grouped    # a list of the DataModels grouped by exposure
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
            except ValueError:
                raise ValueError('{0} must be an ASN file'.format(init))
        else:
            raise TypeError('Input {0!r} is not a list of DataModels or '
                            'an ASN file'.format(init))

        self.assign_group_ids()


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

    def assign_group_ids(self):
        for model in self._models:
            model_attrs = [model.meta.observation.program_number,
                           model.meta.observation.observation_number,
                           model.meta.observation.visit_number,
                           model.meta.observation.visit_group,
                           model.meta.observation.sequence_id,
                           model.meta.observation.activity_id,
                           model.meta.observation.exposure_number,
                           model.meta.instrument.name,
                           model.meta.instrument.channel]
            if None not in model_attrs:
                group_id = ('jw' + "_".join([
                                ''.join(model_attrs[:3]),
                                ''.join(model_attrs[3:6]),
                                model_attrs[6], model_attrs[7].lower(),
                                model_attrs[8].lower()]))
            else:
                root, ext = os.path.splitext(model.meta.filename)
                group_id = "_".join([root, 'group{}'.format(ext)])

            model.meta.group_id = group_id

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
            The path to an ASN file.
        """

        filepath = op.abspath(op.expanduser(op.expandvars(filepath)))
        basedir = op.dirname(filepath)
        try:
            with open(filepath) as asn_file:
                asn_data = Association.load(asn_file)
        except IOError:
            raise IOError("Cannot read ASN file.")

        # make a list of all the input FITS files
        infiles = [op.join(basedir, member['expname']) for member
                   in asn_data['products'][0]['members']]
        try:
            self._models = [datamodel_open(infile, **kwargs) for infile in infiles]
        except IOError:
            raise IOError('Cannot open data models.')

        # Pull the whole association table into meta.asn_table
        self.meta.asn_table = {}
        model_base.properties.merge_tree(self.meta.asn_table._instance, asn_data)

        # populate the output metadata with the output file from the ASN file
        # Should generalize this in the future
        self.meta.resample.output = str(asn_data['products'][0]['name'])
        self.meta.table_name = str(filepath)
        self.meta.pool_name = str(asn_data['asn_pool'])
        self.meta.targname = str(asn_data['target'])
        self.meta.program = str(asn_data['program'])
        self.meta.asn_type = str(asn_data['asn_type'])
        self.meta.asn_rule = str(asn_data['asn_rule'])

    @property
    def models_grouped(self):
        """
        Return a list of lists of DataModels grouped by exposure.
        """

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
        Return a list of names for the DataModel groups by exposure.

        Note that this returns the DataModel "group_id"s appended with
        "_resamp.fits".
        """

        groups = self.models_grouped
        result = []
        for group in groups:
            filename = '{0}_resamp.fits'.format(group[0].meta.group_id)
            result.append(filename)
        return result

    def save(self, filename_not_used, path=None, *args, **kwargs):
        """
        Write out models in container to FITS or ASDF.

        Parameters
        ----------

        filename_not_used : string
            this first argument is ignored in this implementation of the
            save method.  It is used by the pipeline steps to save individual
            files, but that is not applicable here.  Instead, we use the path
            arg below and read the filename output from the meta tag in each
            file in the container.

        path : string
            directory to write out files.  Defaults to current working dir.
            If directory does not exist, it creates it.
        """
        if path is None:
            path = os.getcwd()
        try:
            for model in self._models:
                outpath = op.join(path, model.meta.filename)
                model.save(outpath, *args, **kwargs)
        except IOError as err:
            raise err

    def _get_recursively(self, field, search_dict):
        """
        Takes a dict with nested lists and dicts, and searches all dicts for
        a key of the field provided.
        """
        values_found = []
        for key, value in search_dict.items():
            if key == field:
                values_found.append(value)
            elif isinstance(value, dict):
                results = self._get_recursively(field, value)
                for result in results:
                    values_found.append(result)
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        more_results = self._get_recursively(field, item)
                        for another_result in more_results:
                            values_found.append(another_result)
        return values_found

    def get_recursively(self, field):
        """
        Returns a list of values of the specified field from meta.
        """
        return self._get_recursively(field, self.meta._instance)
