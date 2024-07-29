import io

import asdf
from astropy.io import fits
from stdatamodels.jwst.datamodels.util import open as datamodels_open
from stpipe.library import AbstractModelLibrary, NoGroupID

from jwst.associations import AssociationNotValidError, load_asn

__all__ = ["ModelLibrary"]


class ModelLibrary(AbstractModelLibrary):
    @property
    def crds_observatory(self):
        return "jwst"

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
            with open(asn_path) as asn_file:
                asn_data = load_asn(asn_file)
        except AssociationNotValidError as e:
            raise OSError("Cannot read ASN file.") from e
        return asn_data

    def _filename_to_group_id(self, filename):
        """
        Compute a "group_id" without loading the file as a DataModel

        This function will return the meta.group_id stored in the ASDF
        extension (if it exists) or a group_id calculated from the
        FITS headers.
        """
        # use astropy.io.fits directly to read header keywords
        # avoiding the DataModel overhead
        # TODO look up attribute to keyword in core schema
        with fits.open(filename) as ff:
            if "ASDF" in ff:
                asdf_yaml = asdf.util.load_yaml(io.BytesIO(ff['ASDF'].data.tobytes()))
                if group_id := asdf_yaml.get('meta', {}).get('group_id'):
                    return group_id
            header = ff["PRIMARY"].header
            program_number = header["PROGRAM"]
            observation_number = header["OBSERVTN"]
            visit_number = header["VISIT"]
            visit_group = header["VISITGRP"]
            sequence_id = header["SEQ_ID"]
            activity_id = header["ACT_ID"]
            exposure_number = header["EXPOSURE"]

        # FIXME try except and NoGroupID...
        return _attrs_to_group_id(
            program_number,
            observation_number,
            visit_number,
            visit_group,
            sequence_id,
            activity_id,
            exposure_number,
        )

    def _model_to_group_id(self, model):
        """
        Compute a "group_id" from a model using the DataModel interface
        """
        if group_id := getattr(model.meta, "group_id", None):
            return group_id
        # FIXME try except and NoGroupID...
        return _attrs_to_group_id(
            model.meta.observation.program_number,
            model.meta.observation.observation_number,
            model.meta.observation.visit_number,
            model.meta.observation.visit_group,
            model.meta.observation.sequence_id,
            model.meta.observation.activity_id,
            model.meta.observation.exposure_number,
        )

    def _assign_member_to_model(self, model, member):
        for attr in ("group_id", "tweakreg_catalog", "exptype"):
            if attr in member:
                setattr(model.meta, attr, member[attr])
        if not hasattr(model.meta, "asn"):
            model.meta["asn"] = {}

        model.meta.asn.table_name = self.asn.get("table_name", "")
        model.meta.asn.pool_name = self.asn.get("asn_pool", "")


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
    Combine a number of file metadata values into a ``group_id`` string
    """
    return (
        f"jw{program_number}{observation_number}{visit_number}"
        f"_{visit_group}{sequence_id}{activity_id}"
        f"_{exposure_number}"
    )
