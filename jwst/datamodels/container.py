import copy
from collections import OrderedDict
from collections.abc import Sequence
import os.path as op
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

_ONE_MB = 1 << 20
RECOGNIZED_MEMBER_FIELDS = ['tweakreg_catalog', 'group_id']

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

        - None: initializes an empty `ModelContainer` instance, to which
          DataModels can be added via the ``append()`` method.

    asn_exptypes: str
        list of exposure types from the asn file to read
        into the ModelContainer, if None read all the given files.

    asn_n_members : int
        Open only the first N qualifying members.

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
    >>> c.append(m)

    Notes
    -----
        The optional paramters ``save_open`` and ``return_open`` can be
        provided to control how the `JwstDataModel` are used by the
        :py:class:`ModelContainer`. If ``save_open`` is set to `False`, each input
        `JwstDataModel` instance in ``init`` will be written out to disk and
        closed, then only the filename for the `JwstDataModel` will be used to
        initialize the :py:class:`ModelContainer` object.
        Subsequent access of each member will then open the `JwstDataModel` file to
        work with it. If ``return_open`` is also `False`, then the `JwstDataModel`
        will be closed when access to the `JwstDataModel` is completed. The use of
        these parameters can minimize the amount of memory used by this object
        during processing, with these parameters being used
        by :py:class:`~jwst.outlier_detection.OutlierDetectionStep`.

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

        self._models = []
        self._iscopy = iscopy
        self.asn_exptypes = asn_exptypes
        self.asn_n_members = asn_n_members
        self.asn_table = {}
        self.asn_table_name = None
        self.asn_pool_name = None

        self._memmap = kwargs.get("memmap", False)
        self._return_open = kwargs.get('return_open', True)
        self._save_open = kwargs.get('save_open', True)

        if init is None:
            # Don't populate the container with models
            pass
        elif isinstance(init, fits.HDUList):
            if self._save_open:
                model = [datamodel_open(init, memmap=self._memmap)]
            else:
                model = init._file.name
                init.close()
            self._models.append(model)
        elif isinstance(init, list):
            if all(isinstance(x, (str, fits.HDUList, JwstDataModel)) for x in init):
                if self._save_open:
                    init = [datamodel_open(m, memmap=self._memmap) for m in init]
            else:
                raise TypeError("list must contain items that can be opened "
                                "with jwst.datamodels.open()")
            self._models = init
        elif isinstance(init, self.__class__):
            instance = copy.deepcopy(init._instance)
            self._schema = init._schema
            self._shape = init._shape
            self._asdf = AsdfFile(instance)
            self._instance = instance
            self._ctx = self
            self._models = init._models
            self._iscopy = True
        elif is_association(init):
            self.from_asn(init)
        elif isinstance(init, str):
            init_from_asn = self.read_asn(init)
            self.from_asn(init_from_asn, asn_file_path=init)
        else:
            raise TypeError('Input {0!r} is not a list of JwstDataModels or '
                            'an ASN file'.format(init))

    def __len__(self):
        return len(self._models)

    def __getitem__(self, index):
        m = self._models[index]
        if not isinstance(m, JwstDataModel) and self._return_open:
            m = datamodel_open(m, memmap=self._memmap)
        return m

    def __setitem__(self, index, model):
        self._models[index] = model

    def __delitem__(self, index):
        del self._models[index]

    def __iter__(self):
        for model in self._models:
            if not isinstance(model, JwstDataModel) and self._return_open:
                model = datamodel_open(model, memmap=self._memmap)
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
        result._schema = self._schema
        result._ctx = result
        for m in self._models:
            if isinstance(m, JwstDataModel):
                result.append(m.copy())
            else:
                result.append(m)
        return result

    @staticmethod
    def read_asn(filepath):
        """
        Load fits files from a JWST association file.

        Parameters
        ----------
        filepath : str
            The path to an association file.
        """
        # Prevent circular import:
        from ..associations import AssociationNotValidError, load_asn

        filepath = op.abspath(op.expanduser(op.expandvars(filepath)))
        try:
            with open(filepath) as asn_file:
                asn_data = load_asn(asn_file)
        except AssociationNotValidError as e:
            raise IOError("Cannot read ASN file.") from e
        return asn_data

    def from_asn(self, asn_data, asn_file_path=None):
        """
        Load fits files from a JWST association file.

        Parameters
        ----------
        asn_data : ~jwst.associations.Association
            An association dictionary

        asn_file_path: str
            Filepath of the association, if known.
        """
        # match the asn_exptypes to the exptype in the association and retain
        # only those file that match, as a list, if asn_exptypes is set to none
        # grab all the files
        if self.asn_exptypes:
            infiles = []
            logger.debug('Filtering datasets based on allowed exptypes {}:'
                         .format(self.asn_exptypes))
            for member in asn_data['products'][0]['members']:
                if any([x for x in self.asn_exptypes if re.match(member['exptype'],
                                                                 x, re.IGNORECASE)]):
                    infiles.append(member)
                    logger.debug('Files accepted for processing {}:'.format(member['expname']))
        else:
            infiles = [member for member
                       in asn_data['products'][0]['members']]

        if asn_file_path:
            asn_dir = op.dirname(asn_file_path)
        else:
            asn_dir = ''

        # Only handle the specified number of members.
        if self.asn_n_members:
            sublist = infiles[:self.asn_n_members]
        else:
            sublist = infiles
        try:
            for member in sublist:
                filepath = op.join(asn_dir, member['expname'])
                update_model = any(attr in member for attr in RECOGNIZED_MEMBER_FIELDS)
                if update_model or self._save_open:
                    m = datamodel_open(filepath, memmap=self._memmap)
                    m.meta.asn.exptype = member['exptype']
                    for attr, val in member.items():
                        if attr in RECOGNIZED_MEMBER_FIELDS:
                            if attr == 'tweakreg_catalog':
                                if val.strip():
                                    val = op.join(asn_dir, val)
                                else:
                                    val = None

                            setattr(m.meta, attr, val)

                    if not self._save_open:
                        m.save(filepath, overwrite=True)
                        m.close()
                else:
                    m = filepath

                self._models.append(m)

        except IOError:
            self.close()
            raise

        # Pull the whole association table into meta.asn_table
        self.meta.asn_table = {}
        properties.merge_tree(
            self.meta.asn_table._instance, asn_data
        )

        if asn_file_path is not None:
            self.asn_table_name = op.basename(asn_file_path)
            self.asn_pool_name = asn_data['asn_pool']
            for model in self:
                try:
                    model.meta.asn.table_name = self.asn_table_name
                    model.meta.asn.pool_name = self.asn_pool_name
                except AttributeError:
                    pass

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
                outpath, filename = op.split(
                    path(model.meta.filename, idx=idx)
                )
                if dir_path:
                    outpath = dir_path
                save_path = op.join(outpath, filename)
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
        for i, model in enumerate(self._models):
            params = []
            if not self._save_open:
                model = datamodel_open(model, memmap=self._memmap)

            if (hasattr(model.meta, 'group_id') and
                        model.meta.group_id not in [None, '']):
                group_id = model.meta.group_id

            else:
                for param in unique_exposure_parameters:
                    params.append(getattr(model.meta.observation, param))
                try:
                    group_id = (
                        'jw' + '_'.join(
                            [
                                ''.join(params[:3]),
                                ''.join(params[3:6]),
                                params[6],
                            ]
                        )
                    )
                    model.meta.group_id = group_id
                except TypeError:
                    model.meta.group_id = 'exposure{0:04d}'.format(i + 1)

                group_id = model.meta.group_id
                if not self._save_open and not self._return_open:
                    model.close()
                    model = self._models[i]

            if group_id in group_dict:
                group_dict[group_id].append(model)
            else:
                group_dict[group_id] = [model]

        return group_dict.values()

    @property
    def group_names(self):
        """
        Return list of names for the JwstDataModel groups by exposure.
        """
        result = []
        for group in self.models_grouped:
            result.append(group[0].meta.group_id)
        return result

    def close(self):
        """Close all datamodels."""
        if not self._iscopy:
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
        # Eventually ModelContainer will also be used for Roman, but this
        # will work for now:
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
        for exposure in self.meta.asn_table.products[0].members:
            if exposure.exptype.upper() == "SCIENCE":
                first_exposure = exposure.expname
                break
        else:
            first_exposure = self.meta.asn_table.products[0].members[0].expname

        return datamodel_open(first_exposure)

    def ind_asn_type(self, asn_exptype):
        """
        Determine the indices of models corresponding to ``asn_exptype``.

        Parameters
        ----------
        asn_exptype : str
            Exposure type as defined in an association, e.g. "science".

        Returns
        -------
        ind : list
            Indices of models in ModelContainer._models matching ``asn_exptype``.
        """
        ind = []
        for i, model in enumerate(self._models):
            if model.meta.asn.exptype.lower() == asn_exptype:
                ind.append(i)
        return ind

    def set_buffer(self, buffer_size, overlap=None):
        """Set buffer size for scrolling section-by-section access.

        Parameters
        ----------
        buffer_size : float, None
            Define size of buffer in MB for each section.
            If `None`, a default buffer size of 1MB will be used.

        overlap : int, optional
            Define the number of rows of overlaps between sections.
            If `None`, no overlap will be used.
        """
        self.overlap = 0 if overlap is None else overlap
        self.grow = 0

        with datamodel_open(self._models[0]) as model:
            imrows, imcols = model.data.shape
            data_item_size = model.data.itemsize
            data_item_type = model.data.dtype
            model.close()
        del model
        min_buffer_size = imcols * data_item_size

        self.buffer_size = min_buffer_size if buffer_size is None else (buffer_size * _ONE_MB)

        section_nrows = min(imrows, int(self.buffer_size // min_buffer_size))

        if section_nrows == 0:
            self.buffer_size = min_buffer_size
            logger.warning("WARNING: Buffer size is too small to hold a single row."
                           f"Increasing buffer size to {self.buffer_size / _ONE_MB}MB")
            section_nrows = 1

        nbr = section_nrows - self.overlap
        nsec = (imrows - self.overlap) // nbr
        if (imrows - self.overlap) % nbr > 0:
            nsec += 1

        self.n_sections = nsec
        self.nbr = nbr
        self.section_nrows = section_nrows
        self.imrows = imrows
        self.imcols = imcols
        self.imtype = data_item_type

    def get_sections(self):
        """Iterator to return the sections from all members of the container."""

        for k in range(self.n_sections):
            e1 = k * self.nbr
            e2 = e1 + self.section_nrows

            if k == self.n_sections - 1:  # last section
                e2 = min(e2, self.imrows)
                e1 = min(e1, e2 - self.overlap - 1)

            data_list = np.empty((len(self._models), e2 - e1, self.imcols),
                                 dtype=self.imtype)
            wht_list = np.empty((len(self._models), e2 - e1, self.imcols),
                                dtype=self.imtype)
            for i, model in enumerate(self._models):
                model = datamodel_open(model, memmap=self._memmap)

                data_list[i, :, :] = model.data[e1:e2].copy()
                wht_list[i, :, :] = model.wht[e1:e2].copy()
                model.close()
                del model

            yield (data_list, wht_list, (e1, e2))


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
