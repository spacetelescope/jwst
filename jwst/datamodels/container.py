import copy
from collections import OrderedDict
from collections.abc import Sequence
import json
import os
import re
import logging

import numpy as np

from asdf import AsdfFile
from astropy.io import fits
from stdatamodels import properties

from stdatamodels.jwst.datamodels.model_base import JwstDataModel
from stdatamodels.jwst.datamodels.util import open as datamodel_open
from stdatamodels.jwst.datamodels.util import is_association

__doctest_skip__ = ['ModelContainer']

__all__ = ['ModelContainer']

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class ModelContainer(JwstDataModel, Sequence):
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

        - None: initializes an empty `ModelContainer` instance

    iscopy : bool
        Presume this model is a copy. Members will not be closed
        when the model is closed/garbage-collected.

    Examples
    --------
    >>> container = ModelContainer('example_asn.json')
    >>> for model in container:
    ...     print(model.meta.filename)

    Say the association was a NIRCam dithered dataset. The `models_grouped`
    attribute is a list of lists, the first index giving the list of exposure
    groups, with the second giving the individual datamodels representing
    each detector in the exposure (2 or 8 in the case of NIRCam).

    >>> total_exposure_time = 0.0
    >>> for group in container.models_grouped:
    ...     total_exposure_time += group[0].meta.exposure.exposure_time

    >>> c = ModelContainer()
    >>> m = datamodels.open('myfile.fits')

    Notes
    -----
        When ASN table's members contain attributes listed in
        :py:data:`RECOGNIZED_MEMBER_FIELDS`, :py:class:`ModelContainer` will
        read those attribute values and update the corresponding attributes
        in the ``meta`` of input models.

        .. code-block::
            :caption: Example of ASN table with additional model attributes \
to supply custom catalogs.

            "products": [
                {
                    "name": "resampled_image",
                    "members": [
                        {
                            "expname": "input_image1_cal.fits",
                            "exptype": "science",
                            "tweakreg_catalog": "custom_catalog1.ecsv",
                            "group_id": "custom_group_id_number_1",
                        },
                        {
                            "expname": "input_image2_cal.fits",
                            "exptype": "science",
                            "tweakreg_catalog": "custom_catalog2.ecsv",
                            "group_id": 2
                        },
                        {
                            "expname": "input_image3_cal.fits",
                            "exptype": "science",
                            "tweakreg_catalog": "custom_catalog3.ecsv",
                            "group_id": Null
                        },
                    ]
                }
            ]

        .. warning::
            Input files will be updated in-place with new ``meta`` attribute
            values when ASN table's members contain additional attributes.

        .. warning::
            Custom ``group_id`` affects how models are grouped **both** for
            ``tweakreg`` and ``skymatch`` steps. If one wants to group models
            in one way for the ``tweakreg`` step and in a different way for the
            ``skymatch`` step, one will need to run each step separately with
            their own ASN tables.

        .. note::
            ``group_id`` can be an integer, a string, or Null. When ``group_id``
            is `Null`, it is converted to `None` in Python and it will be
            assigned a group ID based on various exposure attributes - see
            ``models_grouped`` property for more details.

    """
    schema_url = None

    def __init__(self, init=None, asn_exptypes=None, asn_n_members=None,
                 iscopy=False, **kwargs):

        super().__init__(init=None, **kwargs)

        # if True, keep models open and in memory
        self._in_memory = kwargs.get("in_memory", True)

        # _models will be the same length as _members
        # and possibly contain:
        # - None: if the model is not loaded
        # - DataModel: if the model is in-memory
        self._models = None
        self._members = None

        self._iscopy = iscopy

        self.asn_table_name = None
        self.asn_pool_name = None

        self._memmap = kwargs.get("memmap", False)

        # if None, assume an empty list of models
        if init is None:
            init = []

        if isinstance(init, list):
            # might be a list of filenames or a list of models
            self._models = []
            for item in init:
                #if isinstance(item, str):
                #    self._models.append(datamodel_open(item, memmap=self._memmap))
                if isinstance(item, JwstDataModel):
                    self._models.append(item)
                else:
                    raise TypeError("list must contain items that can be opened "
                                    "with jwst.datamodels.open()")

            # since these models are already in memory...
            self._in_memory = True
            # don't close them when the container closes
            self._iscopy = True

            # generate a fake association table
            asn_data = _models_to_association(self._models)
            self._from_asn(asn_data)
        elif isinstance(init, self.__class__):
            instance = copy.deepcopy(init._instance)
            self._schema = init._schema
            self._shape = init._shape
            self._asdf = AsdfFile(instance)
            self._instance = instance
            # TODO what if some models were not loaded?
            self._models = init._models
            self._members = init._members
            self.asn_table_name = init.asn_table_name
            self.asn_pool_name = init.asn_pool_name
            self._iscopy = True
            self._inmemory = init._in_memory
        elif is_association(init):
            self._from_asn(init, asn_exptypes=asn_exptypes, asn_n_members=asn_n_members)
        elif isinstance(init, str):
            # assume an input json association file
            from ..associations import AssociationNotValidError, load_asn

            filepath = os.path.abspath(os.path.expanduser(os.path.expandvars(init)))
            try:
                with open(filepath) as asn_file:
                    asn_data = load_asn(asn_file)
            except AssociationNotValidError as e:
                raise IOError("Cannot read ASN file.") from e

            self._from_asn(asn_data, asn_file_path=init, asn_exptypes=asn_exptypes, asn_n_members=asn_n_members)
        else:
            raise TypeError('Input {0!r} is not a list of JwstDataModels or '
                            'an ASN file'.format(init))

        # FIXME stpipe reaches into _models[0], make sure it's loaded
        if asn_n_members == 1 and self._models:
            self._models[0] = self[0]

    def __len__(self):
        # TODO encapsulate these so they are always the same length
        if len(self._models) != len(self._members):
            raise Exception()
        return len(self._models)

    def __getitem__(self, index):
        m = self._models[index]
        if m is None:
            m = self._load_member(index)
            # if _in_memory is True, save the model for later use
            if self._in_memory:
                self._models[index] = m
        return m

    #def __setitem__(self, index, model):
    #    self._models[index] = model

    #def __delitem__(self, index):
    #    del self._models[index]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    #def insert(self, index, model):
    #    self._models.insert(index, model)

    #def append(self, model):
    #    self._models.append(model)

    #def extend(self, model):
    #    self._models.extend(model)

    #def pop(self, index=-1):
    #    self._models.pop(index)

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
        # TODO shouldn't result._iscopy be True?
        result._iscopy = self._iscopy
        result._in_memory = self._in_memory
        result._schema = self._schema
        # TODO is there a cleaner way to copy the models?
        for m in self._models:
            if isinstance(m, JwstDataModel):
                # if the model is copied here, it will get a new FileReference
                # and result._iscopy should be False...
                result._models.append(m.copy())
            else:
                result._models.append(m)
        result._members = copy.deepcopy(self._members)
        result.asn_table_name = self.asn_table_name
        result.asn_pool_name = self.asn_pool_name
        return result

    def _from_asn(self, asn_data, asn_file_path=None, asn_exptypes=None, asn_n_members=None):
        """
        Load fits files from a JWST association file.

        Parameters
        ----------
        asn_data : ~jwst.associations.Association
            An association dictionary

        asn_file_path: str
            Filepath of the association, if known.

        asn_exptypes: str
            list of exposure types from the asn file to read
            into the ModelContainer, if None read all the given files.

        asn_n_members : int
            Open only the first N qualifying members.
        """
        # copy asn_data so it can be modified below
        asn_data = copy.deepcopy(asn_data)

        # match the asn_exptypes to the exptype in the association and retain
        # only those file that match, as a list, if asn_exptypes is set to none
        # grab all the files
        if asn_exptypes:
            members = []
            logger.debug('Filtering datasets based on allowed exptypes {}:'
                         .format(asn_exptypes))
            for member in asn_data['products'][0]['members']:
                if any([x for x in asn_exptypes if re.match(member['exptype'],
                                                                 x, re.IGNORECASE)]):
                    members.append(member)
                    logger.debug('Files accepted for processing {}:'.format(member['expname']))
        else:
            members = [member for member
                       in asn_data['products'][0]['members']]

        if asn_file_path:
            asn_dir = os.path.dirname(asn_file_path)
        else:
            asn_dir = ''

        # Only handle the specified number of members.
        if asn_n_members:
            members = members[:asn_n_members]

        # add a "group_id" (if not already defined in the association)
        for (i, member) in enumerate(members):
            if member["expname"] is None:
                if self._models is None or self._models[i] is None:
                    raise Exception()
                member['_filename'] = None
            else:
                filename = os.path.join(asn_dir, member["expname"])
                member['_filename'] = filename
            if member.get("group_id") is None:
                try:
                    member["group_id"] = _file_to_group_id(filename)
                except (TypeError, AttributeError, KeyError):
                    member["group_id"] = 'exposure{0:04d}'.format(i + 1)

        # Pull the whole association table into meta.asn_table
        self.meta.asn_table = {}
        # TODO why is there a merge with an empty dict?
        properties.merge_tree(
            self.meta.asn_table._instance, asn_data
        )

        self._members = members
        # if we have no models, set them as None to signify they're not loaded
        if not self._models:
            self._models = [None] * len(members)

        if asn_file_path is not None:
            self.asn_table_name = os.path.basename(asn_file_path)
            # TODO why is pool only set if asn_file_path is not None?
            self.asn_pool_name = asn_data['asn_pool']

    def _load_member(self, index):
        member = self._members[index]
        model = datamodel_open(member['_filename'], memmap=self._memmap)

        # when model is opened, overwrite:
        # - exptype to meta.asn.exptype
        # - asn_table_name to meta.asn.table_name (if no AttributeError)
        # - asn_pool_name to meta.asn_pool_name (same if no AttributeError)
        # - group_id to meta.group_id
        model.meta.asn.exptype = member['exptype']
        model.meta.group_id = member['group_id']
        try:
            # TODO when do these fail? should they fail together?
            model.meta.asn.table_name = self.asn_table_name
            model.meta.asn.pool_name = self.asn_pool_name
        except AttributeError:
            pass
        return model

    def save(self,
             path=None,
             dir_path=None,
             save_model_func=None,
             **kwargs):
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
            def path(filename, idx=None):
                return filename
        elif not callable(path):
            path = make_file_with_index

        for idx, model in enumerate(self):
            if len(self) <= 1:
                idx = None
            if save_model_func is None:
                outpath, filename = os.path.split(
                    path(model.meta.filename, idx=idx)
                )
                if dir_path:
                    outpath = dir_path
                save_path = os.path.join(outpath, filename)
                try:
                    output_paths.append(
                        model.save(save_path, **kwargs)
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
        Assign a grouping ID by exposure, if not already assigned.

        If ``model.meta.group_id`` does not exist or it is `None`, then data
        from different detectors of the same exposure will be assigned the
        same group ID, which allows grouping by exposure in the ``tweakreg`` and
        ``skymatch`` steps. The following metadata is used when
        determining grouping:

        meta.observation.program_number
        meta.observation.observation_number
        meta.observation.visit_number
        meta.observation.visit_group
        meta.observation.sequence_id
        meta.observation.activity_id
        meta.observation.exposure_number

        If a model already has ``model.meta.group_id`` set, that value will be
        used for grouping.
        """
        group_dict = {}
        for (i, member) in enumerate(self._members):
            group_id = member['group_id']
            if group_id not in group_dict:
                group_dict[group_id] = []
            group_dict[group_id].append(self[i])

        return list(group_dict.values())

    @property
    def group_names(self):
        """
        Return list of names for the JwstDataModel groups by exposure.
        """
        group_ids = {}  # dictionary as an ordered set
        for member in self._members:
            if member['group_id'] not in group_ids:
                group_ids[member['group_id']] = None
        return list(group_ids.keys())

    def close(self):
        """Close all datamodels."""
        if not self._iscopy and self._models is not None:
            for model in self._models:
                if isinstance(model, JwstDataModel):
                    model.close()

    @property
    def crds_observatory(self):
        """
        Get the CRDS observatory for this container.  Used when selecting
        step/pipeline parameter files when the container is a pipeline input.

        Returns
        -------
        str
        """
        return "jwst"

    def get_crds_parameters(self):
        """
        Get CRDS parameters for this container.  Used when selecting
        step/pipeline parameter files when the container is a pipeline input.

        Returns
        -------
        dict
        """
        with self._open_first_science_exposure() as model:
            return model.get_crds_parameters()

    def _open_first_science_exposure(self):
        """
        Open first model with exptype SCIENCE, or the first model
        if none exists.

        Returns
        -------
        stdatamodels.JwstDataModel
        """
        for (i, member) in enumerate(self._members):
            if member['exptype'].upper() == "SCIENCE":
                return self[i]
        return self[0]


def _attrs_to_group_id(
        program_number,
        observation_number,
        visit_number,
        visit_group,
        sequence_id,
        activity_id,
        exposure_number,
    ):
    return (
        f"jw{program_number}{observation_number}{visit_number}"
        f"_{visit_group}{sequence_id}{activity_id}"
        f"_{exposure_number}"
    )


def _file_to_group_id(filename):
    """
    Compute a "group_id" without loading the file
    as a DataModel
    """
    # use astropy.io.fits directly to read header keywords
    # avoiding the DataModel overhead
    # TODO look up attribute to keyword in core schema
    with fits.open(filename) as ff:
        header = ff["PRIMARY"].header
        program_number = header["PROGRAM"]
        observation_number = header["OBSERVTN"]
        visit_number = header["VISIT"]
        visit_group = header["VISITGRP"]
        sequence_id = header["SEQ_ID"]
        activity_id = header["ACT_ID"]
        exposure_number = header["EXPOSURE"]

    return _attrs_to_group_id(
        program_number,
        observation_number,
        visit_number,
        visit_group,
        sequence_id,
        activity_id,
        exposure_number,
    )


def _model_to_group_id(model):
    """
    Compute a "group_id" from a model
    """
    return _attrs_to_group_id(
        model.meta.observation.program_number,
        model.meta.observation.observation_number,
        model.meta.observation.visit_number,
        model.meta.observation.visit_group,
        model.meta.observation.sequence_id,
        model.meta.observation.activity_id,
        model.meta.observation.exposure_number,
    )


def _models_to_association(models, meta=None, member_meta=None):
    """
    Create an association table from a list of models
    What type of association? What's required?

    Both lvl2 and lvl3 require:
        - asn_id [string] special naming convention (maybe default to None?)
        - asn_pool [string] special naming? (maybe default to 'undetermined' or None?)
        - products [list]
            - name (optional) [string] (is this actually optional? calwebb3_image3 at least catches errors when it's missing)
            - members [list] each item a dict
                - expname [string] filename
                - exptype [string] science, background, ...
    """
    if meta is None:
        meta = {}
    if member_meta is None:
        member_meta = [{}] * len(models)

    if len(member_meta) != len(models):
        raise ValueError(f"len(member_meta)[{len(member_meta)}] != len(models)[{len(models)}]")

    asn_table = {
        "asn_id": None,
        "asn_pool": None,
    }
    asn_table |= meta

    # for each model generate:
    # - expname (filename)
    # - exptype (probably all science...)
    # - group_id (see above)
    # use this to populate members filenames
    members = []
    for (i, model) in enumerate(models):
        member = {
            "exptype": "science",
        }

        member |= member_meta[i]

        if "expname" not in member:
            member["expname"] = model.meta.filename

        if "group_id" not in member:
            try:
                member["group_id"] = model.meta.group_id
            except AttributeError:
                member["group_id"] = None
            if member["group_id"] is None:
                try:
                    member["group_id"] = _model_to_group_id(model)
                except (TypeError, AttributeError):
                    member["group_id"] = 'exposure{0:04d}'.format(i + 1)

        members.append(member)

    # add members to table as first product
    asn_table["products"] = [{"members": members}]

    return asn_table


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
    path_head, path_tail = os.path.split(file_path)
    base, ext = os.path.splitext(path_tail)
    if idx is not None:
        base = base + str(idx)
    return os.path.join(path_head, base + ext)
