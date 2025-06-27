from collections import OrderedDict
from collections.abc import Sequence
import copy
import os.path as op
from pathlib import Path
import re
import logging
from astropy.io import fits

from stdatamodels.jwst.datamodels.model_base import JwstDataModel
from stdatamodels.jwst.datamodels.util import open as datamodel_open
from stdatamodels.jwst.datamodels.util import is_association

from jwst.datamodels.utils import attrs_to_group_id

__doctest_skip__ = ["ModelContainer"]

__all__ = ["ModelContainer"]

RECOGNIZED_MEMBER_FIELDS = ["tweakreg_catalog", "group_id"]
EMPTY_ASN_TABLE = {
    "asn_id": None,
    "asn_pool": None,
    "products": [{"name": "", "members": [{"exptype": "", "expname": ""}]}],
}

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class ModelContainer(Sequence):
    """
    A list-like container for holding DataModels.

    This functions like a list for holding DataModel objects.  It can be
    iterated through like a list, DataModels within the container can be
    addressed by index, and the datamodels can be grouped into a list of
    lists for grouped looping, useful for NIRCam where grouping together
    all detectors of a given exposure is useful for some pipeline steps.

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
    """

    def __init__(self, init=None, asn_exptypes=None, asn_n_members=None, **kwargs):  # noqa: ARG002
        """
        Initialize the container.

        Parameters
        ----------
        init : file path, list of DataModels, or None
            If a file path, initialize from an association table.
            If a list, can be a list of DataModels of any type
            If None, initializes an empty `ModelContainer` instance, to which
            DataModels can be added via the ``append()`` method.

        asn_exptypes : str
            List of exposure types from the asn file to read
            into the ModelContainer, if None read all the given files.

        asn_n_members : int
            Open only the first N qualifying members.
        """
        self._models = []
        self.asn_exptypes = asn_exptypes
        self.asn_n_members = asn_n_members
        self.asn_table = copy.deepcopy(EMPTY_ASN_TABLE)
        self.asn_table_name = None
        self.asn_pool_name = None
        self.asn_file_path = None

        if init is None:
            # Don't populate the container with models
            pass
        elif isinstance(init, list):
            if all(isinstance(x, (str, fits.HDUList, JwstDataModel)) for x in init):
                for m in init:
                    self._models.append(datamodel_open(m))
                # set asn_table_name and product name to first datamodel stem
                # since they were not provided
                fname = self._models[0].meta.filename
                if fname is not None:
                    root = Path(fname).name.split(".")[0]
                    default_name = "_".join(root.split("_")[:-1])  # remove old suffix
                else:
                    default_name = ""
                self.asn_table_name = default_name
                self.asn_table["products"][0]["name"] = default_name
            else:
                raise TypeError(
                    "list must contain items that can be opened with jwst.datamodels.open()"
                )
        elif isinstance(init, self.__class__):
            for m in init:
                self._models.append(datamodel_open(m))
            self.asn_exptypes = init.asn_exptypes
            self.asn_n_members = init.asn_n_members
            self.asn_table = init.asn_table
            self.asn_table_name = init.asn_table_name
            self.asn_pool_name = init.asn_pool_name
            self.asn_file_path = init.asn_file_path
        elif is_association(init):
            self.from_asn(init)
        elif isinstance(init, (str, Path)):
            init_from_asn = self.read_asn(init)
            self.asn_file_path = init
            self.from_asn(init_from_asn)
        else:
            raise TypeError(f"Input {init} is not a list of JwstDataModels or an ASN file")

    def __len__(self):
        return len(self._models)

    def __getitem__(self, index):
        return self._models[index]

    def __setitem__(self, index, model):
        self._models[index] = model

    def __delitem__(self, index):
        del self._models[index]

    def __iter__(self):
        yield from self._models

    def insert(self, index, model):  # noqa: D102
        self._models.insert(index, model)

    def append(self, model):  # noqa: D102
        self._models.append(model)

    def extend(self, model):  # noqa: D102
        self._models.extend(model)

    def pop(self, index=-1):  # noqa: D102
        self._models.pop(index)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def copy(self, memo=None):
        """
        Make a deep copy of the container.

        Parameters
        ----------
        memo : dict
            Keeps track of elements that have already been copied to avoid infinite recursion.

        Returns
        -------
        ModelContainer
            A deep copy of the container and all the models in it.
        """
        result = self.__class__(init=None)
        for m in self._models:
            result.append(m.copy(memo=memo))
        return result

    @staticmethod
    def read_asn(filepath):
        """
        Load fits files from a JWST association file.

        Parameters
        ----------
        filepath : str
            The path to an association file.

        Returns
        -------
        dict
            An association dictionary
        """
        # Prevent circular import:
        from jwst.associations import AssociationNotValidError, load_asn

        filepath = Path(op.expandvars(filepath)).expanduser().resolve()
        try:
            with Path(filepath).open() as asn_file:
                asn_data = load_asn(asn_file)
        except AssociationNotValidError as e:
            raise OSError("Cannot read ASN file.") from e
        return asn_data

    def from_asn(self, asn_data):
        """
        Load fits files from a JWST association file.

        Parameters
        ----------
        asn_data : ~jwst.associations.Association
            An association dictionary
        """
        # match the asn_exptypes to the exptype in the association and retain
        # only those file that match, as a list, if asn_exptypes is set to none
        # grab all the files
        if self.asn_exptypes:
            infiles = []
            logger.debug(f"Filtering datasets based on allowed exptypes {self.asn_exptypes}:")
            for member in asn_data["products"][0]["members"]:
                if any(re.match(member["exptype"], x, re.IGNORECASE) for x in self.asn_exptypes):
                    infiles.append(member)
                    logger.debug("Files accepted for processing {}:".format(member["expname"]))
        else:
            infiles = list(asn_data["products"][0]["members"])

        if self.asn_file_path:
            asn_dir = Path(self.asn_file_path).parent
        else:
            asn_dir = Path()

        # Only handle the specified number of members.
        if self.asn_n_members:
            sublist = infiles[: self.asn_n_members]
        else:
            sublist = infiles
        try:
            for member in sublist:
                filepath = asn_dir / member["expname"]
                m = datamodel_open(filepath)
                m.meta.asn.exptype = member["exptype"]
                for attr, val in member.items():
                    if attr in RECOGNIZED_MEMBER_FIELDS:
                        if attr == "tweakreg_catalog":
                            if val.strip():
                                val = asn_dir / val
                            else:
                                val = None

                        setattr(m.meta, attr, val)
                self._models.append(m)

        except OSError:
            self.close()
            raise

        # Pull the whole association table into the asn_table attribute
        self.asn_table = copy.deepcopy(asn_data)

        if self.asn_file_path is not None:
            self.asn_table_name = Path(self.asn_file_path).name
            self.asn_pool_name = asn_data["asn_pool"]
            for model in self:
                try:
                    model.meta.asn.table_name = self.asn_table_name
                    model.meta.asn.pool_name = self.asn_pool_name
                except AttributeError:
                    pass

    def save(self, path=None, save_model_func=None, **kwargs):
        """
        Write out models in container to FITS or ASDF.

        Parameters
        ----------
        path : str or None
            - If None, the `meta.filename` is used for each model.
            - If a string, the string is used as a root and an index is
              appended, along with the '.fits' extension.

        save_model_func : func or None
            Alternate function to save each model instead of
            the models `save` method. Takes one argument, the model,
            and keyword argument `idx` for an index.

        **kwargs : dict
            Additional parameters to be passed to the `save` method of each
            model.

        Returns
        -------
        output_paths : [str[, ...]]
            List of output file paths of where the models were saved.
        """
        output_paths = []
        for idx, model in enumerate(self):
            if save_model_func is None:
                if path is None:
                    save_path = model.meta.filename
                else:
                    if len(self) <= 1:
                        idx = ""
                    if path.endswith(".fits"):
                        save_path = path.replace(".fits", f"{idx}.fits")
                    else:
                        save_path = f"{path}{idx}.fits"
                output_paths.append(model.save(save_path, **kwargs))
            else:
                output_paths.append(save_model_func(model, idx=idx))
        return output_paths

    @property
    def models_grouped(self):
        """
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

        Returns
        -------
        list
            A list of lists of datamodels grouped by exposure.
        """
        group_dict = OrderedDict()
        for i, model in enumerate(self._models):
            if hasattr(model.meta, "group_id") and model.meta.group_id not in [None, ""]:
                group_id = model.meta.group_id

            else:
                try:
                    group_id = attrs_to_group_id(model.meta.observation)
                except KeyError:
                    # If the required keys are not present, assign a default group ID
                    group_id = f"exposure{i + 1:04d}"

                model.meta.group_id = group_id

            if group_id in group_dict:
                group_dict[group_id].append(model)
            else:
                group_dict[group_id] = [model]

        return group_dict.values()

    @property
    def group_names(self):
        """
        List all the group names in the container.

        Returns
        -------
        list
            A list of group names.
        """
        result = []
        for group in self.models_grouped:
            result.append(group[0].meta.group_id)
        return result

    def close(self):
        """Close all datamodels."""
        for model in self._models:
            if isinstance(model, JwstDataModel):
                model.close()

    @property
    def crds_observatory(self):
        """
        Return the observatory name for CRDS queries.

        Returns
        -------
        str
            The observatory name for CRDS queries.
        """
        return "jwst"

    def get_crds_parameters(self):
        """
        Get CRDS parameters for this container.

        Notes
        -----
        stpipe requires ModelContainer to have a crds_observatory attribute in order
        to pass through step.run(), but it is never accessed.
        """
        msg = (
            "stpipe uses the get_crds_parameters method from the 0th model in the "
            "ModelContainer. This method is currently not used."
        )
        raise NotImplementedError(msg)

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
