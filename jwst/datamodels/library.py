import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
from asdf.tags.core import NDArrayType
from astropy.time import Time
from stdatamodels.jwst.datamodels import read_metadata
from stdatamodels.jwst.datamodels.util import open as datamodels_open
from stpipe.library import AbstractModelLibrary, BorrowError, NoGroupID

from jwst.associations import AssociationNotValidError, load_asn
from jwst.datamodels.utils import attrs_to_group_id

__all__ = ["ModelLibrary"]


class ModelLibrary(AbstractModelLibrary):
    """
    JWST implementation of the ModelLibrary.

    ModelLibrary is a container designed to allow
    efficient processing of datamodel instances created from an association.
    See the `stpipe library documentation <https://stpipe.readthedocs.io/en/latest/model_library.html>`_
    for more information.
    """

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

    @property
    def exptypes(self):
        """
        List exposure types for all members in the library.

        Returns
        -------
        list
            List of exposure types for all members in the library.
        """
        return [member["exptype"] for member in self._members]

    @property
    def asn_dir(self):
        """
        Return the directory from which the association was loaded.

        Returns
        -------
        str
            The directory from which the association was loaded.
        """
        return str(self._asn_dir)

    @property
    def on_disk(self):
        """
        Return the library's on_disk attribute.

        If True, the library is using temporary files to store the models when not in use.
        If False, the library is storing the models in memory.

        Returns
        -------
        bool
            Whether the library is on disk.
        """
        return self._on_disk

    def indices_for_exptype(self, exptype):
        """
        Determine the indices of models corresponding to ``exptype``.

        Parameters
        ----------
        exptype : str
            Exposure type as defined in an association, e.g. "science". case-insensitive

        Returns
        -------
        ind : list
            Indices of models in ModelLibrary with member exposure types matching ``exptype``.

        Notes
        -----
        Library does NOT need to be open (i.e., this can be called outside the `with` context)
        """
        return [
            i
            for i, member in enumerate(self._members)
            if member["exptype"].lower() == exptype.lower()
        ]

    def _model_to_filename(self, model):
        model_filename = model.meta.filename
        if model_filename is None:
            model_filename = "model.fits"
        return model_filename

    def _datamodels_open(self, filename, **kwargs):
        return datamodels_open(filename, **kwargs)

    @classmethod
    def _load_asn(cls, asn_path):
        try:
            with Path(asn_path).open() as asn_file:
                asn_data = load_asn(asn_file)
        except AssociationNotValidError as e:
            raise OSError("Cannot read ASN file.") from e
        return asn_data

    def _filename_to_group_id(self, filename):
        """
        Compute a "group_id" without loading the file as a DataModel.

        Parameters
        ----------
        filename : str
            The name of the file to read.

        Returns
        -------
        str
            The meta.group_id stored in the ASDF extension (if it exists)
            or a group_id calculated from the FITS headers.
        """
        # use read_metadata to read header keywords
        # avoiding the DataModel overhead
        meta = read_metadata(filename, flatten=False)["meta"]
        if "group_id" in meta.keys():
            return meta["group_id"]
        try:
            return attrs_to_group_id(meta["observation"])

        except KeyError as e:
            msg = f"Cannot find header keyword {e} in {filename}"
            raise NoGroupID(msg) from e

    def _model_to_exptype(self, model):
        return model.meta.asn.exptype

    def _model_to_group_id(self, model):
        """
        Compute a "group_id" from a model using the DataModel interface.

        Parameters
        ----------
        model : DataModel
            The model from which to to extract the group_id.

        Returns
        -------
        str
            The group_id string.

        Raises
        ------
        NoGroupID
            If the model does not have a meta.group_id, and one cannot be built from
            the model's meta.observation attributes.
        """
        if group_id := getattr(model.meta, "group_id", None):
            return group_id
        if model.meta.hasattr("observation"):
            try:
                return attrs_to_group_id(model.meta.observation)
            except KeyError as e:
                raise NoGroupID(f"Cannot build group_id from model.meta.observation: {e}") from e
        raise NoGroupID(f"{model} missing group_id: meta.observation was not found.")

    def _assign_member_to_model(self, model, member):
        model.meta.asn.exptype = member["exptype"]
        for attr in ("group_id", "tweakreg_catalog"):
            if attr in member:
                setattr(model.meta, attr, member[attr])
        if not hasattr(model.meta, "asn"):
            model.meta["asn"] = {}

        if "table_name" in self.asn.keys():
            model.meta.asn.table_name = self.asn["table_name"]

        if "asn_pool" in self.asn.keys():  # do not clobber existing values
            model.meta.asn.pool_name = self.asn["asn_pool"]

    def get_crds_parameters(self):
        """
        Override base method so lazy-loading of metadata can be used to get CRDS parameters.

        Get the "crds_parameters" from either:
            - the first "science" member (based on exptype)
            - the first model (if no "science" member is found)

        If no "science" members are found in the library a ``UserWarning``
        will be issued.

        Returns
        -------
        crds_parameters : dict
            The result of ``get_crds_parameters`` called on the selected
            model.
        """
        with self:
            idx = None
            for i, member in enumerate(self._members):
                if member["exptype"].lower() == "science":
                    idx = i
                    break
            if idx is None:
                warnings.warn(
                    "get_crds_parameters failed to find any science members. "
                    "The first model was used to determine the parameters",
                    UserWarning,
                    stacklevel=2,
                )
                idx = 0

            return self.read_metadata(idx)

    def read_metadata(self, idx, flatten=True):
        """
        Read metadata for a model at the given index.

        Parameters
        ----------
        idx : int
            Index of the model in the library.

        Returns
        -------
        dict
            The metadata dictionary for the model.
        """
        # if model was already borrowed, calling code should just read the metadata
        # from the model itself, so we raise an error here
        if idx in self._ledger:
            raise BorrowError("Attempt to read metadata from model that is already borrowed")

        # find model in _loaded_models, temp_filenames, or asn_dir
        if self._on_disk:
            if idx in self._temp_filenames:
                # if model has been modified, find its temp filename
                filename = self._temp_filenames[idx]
            else:
                # otherwise, find the filename in the asn_dir
                member = self._members[idx]
                filename = Path(self._asn_dir) / member["expname"]
        else:
            if idx in self._loaded_models:
                # if this model is in memory, retrieve parameters from it directly
                model = self._loaded_models[idx]
                return _read_meta_from_open_model(model, flatten=flatten)
            else:
                # otherwise, find the filename in the asn_dir
                member = self._members[idx]
                filename = Path(self._asn_dir) / member["expname"]

        meta = read_metadata(filename, flatten=flatten)
        meta = self._assign_member_to_meta(meta, self._members[idx], flatten)
        return meta

    def _assign_member_to_meta(self, meta, member, flatten):
        """
        Update meta dict with asn-related attributes, similar to _assign_member_to_model.

        Parameters
        ----------
        meta : dict
            The metadata dictionary to update.
        member : dict
            The member dictionary containing association attributes.
        flatten : bool
            If True, the metadata will be flattened to a single level.
            If False, the metadata will be nested.

        Returns
        -------
        dict
            The updated metadata dictionary with association attributes.
        """
        if flatten:
            meta["meta.asn.exptype"] = member["exptype"]
            for attr in ("group_id", "tweakreg_catalog"):
                if attr in member:
                    meta["meta." + attr] = member[attr]
            if "table_name" in self.asn.keys():
                meta["meta.asn.table_name"] = self.asn["table_name"]
            if "asn_pool" in self.asn.keys():
                meta["meta.asn.pool_name"] = self.asn["asn_pool"]
        else:
            meta["meta"].setdefault("asn", {})
            meta["meta"]["asn"]["exptype"] = member["exptype"]
            for attr in ("group_id", "tweakreg_catalog"):
                if attr in member:
                    meta["meta"][attr] = member[attr]
            if "table_name" in self.asn.keys():
                meta["meta"]["asn"]["table_name"] = self.asn["table_name"]
            if "asn_pool" in self.asn.keys():
                meta["meta"]["asn"]["pool_name"] = self.asn["asn_pool"]
        return meta


def _read_meta_from_open_model(model, flatten):
    """
    Read metadata from an open model.

    Parameters
    ----------
    model : DataModel
        The model from which to read metadata.
    flatten : bool
        If True, the metadata will be flattened to a single level.
        If False, the metadata will be nested.

    Returns
    -------
    dict
        The metadata dictionary for the model.
    """
    if flatten:
        return model.get_crds_parameters()

    def convert_val(val):
        if isinstance(val, datetime):
            return val.isoformat()
        elif isinstance(val, Time):
            return str(val)
        return val

    def recurse(tree):
        """
        Walk the tree and add metadata to a new dictionary.

        Do all conversions and exclusions model.get_crds_parameters does,
        but without flattening.

        Parameters
        ----------
        tree : dict
            An asdf-like nested dictionary structure, e.g. model._instance.

        Returns
        -------
        dict
            A new nested dictionary containing the metadata.
        """
        new_tree = {}
        for key, val in tree.items():
            if key == "wcs":
                continue
            val = convert_val(val)
            if isinstance(val, dict):
                new_tree[key] = recurse(val)
            elif isinstance(val, list):
                new_list = []
                for item in val:
                    if isinstance(item, dict):
                        new_list.append(recurse(item))
                    elif isinstance(item, list):
                        # handle nested lists
                        new_list.append(recurse({"_idx": item})["_idx"])  # placeholder key
                    elif isinstance(item, (np.ndarray, NDArrayType)):
                        continue
                    else:
                        new_list.append(convert_val(item))
                new_tree[key] = new_list
            elif isinstance(val, (np.ndarray, NDArrayType)):
                continue
            elif isinstance(val, (str, int, float, complex, bool)):
                new_tree[key] = val
            else:
                continue  # skip unsupported types
        return new_tree

    return recurse(model._instance)  # noqa: SLF001
