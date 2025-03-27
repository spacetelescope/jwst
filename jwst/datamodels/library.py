import io

import asdf
from pathlib import Path
from astropy.io import fits
from stdatamodels.jwst.datamodels.util import open as datamodels_open
from stpipe.library import AbstractModelLibrary, NoGroupID

from jwst.associations import AssociationNotValidError, load_asn

__all__ = ["ModelLibrary"]


class ModelLibrary(AbstractModelLibrary):
    """
    JWST implementation of the ModelLibrary.

    ModelLibrary is a container designed to allow
    efficient processing of datamodel instances created from an association.
    See the `stpipe library documentation <https://stpipe.readthedocs.io/en/latest/model_library.html`
    for details.
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
        # use astropy.io.fits directly to read header keywords
        # avoiding the DataModel overhead
        try:
            with fits.open(filename) as ff:
                if "ASDF" in ff:
                    asdf_yaml = asdf.util.load_yaml(
                        io.BytesIO(ff["ASDF"].data.tobytes()), tagged=True
                    )
                    if group_id := asdf_yaml.get("meta", {}).get("group_id"):
                        return group_id
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
        if hasattr(model.meta, "observation"):
            return _attrs_to_group_id(
                model.meta.observation.program_number,
                model.meta.observation.observation_number,
                model.meta.observation.visit_number,
                model.meta.observation.visit_group,
                model.meta.observation.sequence_id,
                model.meta.observation.activity_id,
                model.meta.observation.exposure_number,
            )
        raise NoGroupID(f"{model} missing group_id")

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


def _attrs_to_group_id(
    program_number,
    observation_number,
    visit_number,
    visit_group,
    sequence_id,
    activity_id,
    exposure_number,
):
    """
    Combine a number of file metadata values into a ``group_id`` string.

    Parameters
    ----------
    program_number : int
        Program number.
    observation_number : int
        Observation number.
    visit_number : int
        Visit number.
    visit_group : str
        Visit group.
    sequence_id : str
        Sequence ID.
    activity_id : str
        Activity ID.
    exposure_number : int
        Exposure number.

    Returns
    -------
    str
        The group_id string.
    """
    for val in (
        program_number,
        observation_number,
        visit_number,
        visit_group,
        sequence_id,
        activity_id,
        exposure_number,
    ):
        if val is None:
            raise NoGroupID(f"Missing required value for group_id: {val}")
    return (
        f"jw{program_number}{observation_number}{visit_number}"
        f"_{visit_group}{sequence_id}{activity_id}"
        f"_{exposure_number}"
    )
