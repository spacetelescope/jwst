from astropy.io import fits
from stdatamodels import fits_support

from .rules import make_blender
from ._schemautil import parse_schema
from ._tablebuilder import TableBuilder, table_to_schema


__all__ = ["ModelBlender"]


class ModelBlender:
    """
    Class to "blend" metadata from several datamodels.

    The output "blended" model will contain:

        - metadata for a combined model
        - a table with metadata of each datamodel

    Input models can be added to the blender using `ModelBlender.accumulate`
    and the output/combined model updated using `ModelBlender.finalize_model`.

    All input/accumulated models must be of the same type.

    >>> blender = ModelBlender()
    >>> blender.accumulate(input_model_a)  # doctest: +SKIP
    >>> blender.accumulate(input_model_b)  # doctest: +SKIP
    >>> blender.finalize_model(combined_model)  # doctest: +SKIP
    """

    def __init__(self, blend_ignore_attrs=None):
        """
        Create a new `ModelBlender`.

        Parameters
        ----------
        blend_ignore_attrs : list or None
            A list of metadata attributes to ignore during blending.
            These attributes will not be set on the output/combined.
            These attributes must be strings containing the dotted
            path of each attribute (for example "meta.filename").
            (Note that "meta.wcs" will always be ignored).
        """
        self._model_type = None
        self._first_header_meta = None
        self._blenders = None
        self._table_builder = None
        self._blend_ignore_attrs = ["meta.wcs"]
        if blend_ignore_attrs is not None:
            self._blend_ignore_attrs.extend(blend_ignore_attrs)

    def accumulate(self, model):
        """
        Blend metadata for model.

        Process model adding its metadata to the blended
        metadata and the metadata table.

        Parameters
        ----------
        model : `jwst.datamodels.JwstDataModel`
            The datamodel to blend.
        """
        if self._first_header_meta is None:
            self._model_type = type(model)
            # search the schema for other metadata to "blend" and to add to the table
            attr_to_columns, attr_to_blend_rules, schema_ignores = parse_schema(model.schema)

            # update ignores list for items in schema that can't be blended
            self._blend_ignore_attrs.extend(schema_ignores)

            # capture the entire contents of the first model metadata
            self._first_header_meta = {}
            for attr, v in model.to_flat_dict(include_arrays=False).items():
                if not attr.startswith("meta"):
                    continue
                if any(attr.startswith(i) for i in self._blend_ignore_attrs):
                    continue
                self._first_header_meta[attr] = v

            # make "blenders" for the metadata with special rules
            self._blenders = {}
            for attr, rule in attr_to_blend_rules.items():
                if rule == "first":
                    continue
                if any(attr.startswith(i) for i in self._blend_ignore_attrs):
                    continue
                self._blenders[attr] = make_blender(rule)

            # make a table builder using the mapping from the schema
            self._table_builder = TableBuilder(attr_to_columns)
        else:
            if type(model) != self._model_type:  # noqa: E721
                raise ValueError(
                    f"model of type {type(model)} "
                    f"does not match previous type({self._model_type}). "
                    "ModelBlender only supports blending models of the same model type."
                )

        # convert the model to a flat header
        header = model.to_flat_dict(include_arrays=False)

        # add the header to the table
        self._table_builder.header_to_row(header)

        # and perform any special blending
        for attr, blender in self._blenders.items():
            if attr in header:
                blender.accumulate(header[attr])

    def _finalize_metadata(self):
        # start with the entire contents of the first model
        meta = self._first_header_meta.copy()
        for attr, blender in self._blenders.items():
            meta[attr] = blender.finalize()
        return meta

    def _finalize_table(self):
        return self._table_builder.build_table()

    def finalize_model(self, model):
        """
        Update model with the blend results.

        Add blended metadata and the accumulated metadata table to
        the provided datamodel. The update process involves:

            - setting the model metadata to the blended metadata values
            - adding an "hdrtab" attribute (containing the metadata table)
            - updating the model schema to save "hdrtab"

        The provided model will be updated in-place.

        Parameters
        ----------
        model : `jwst.datamodels.JwstDataModel`
            A datamodel that will have its metadata set
            to the blended metadata and have the metadata
            table assigned to the "hdrtab" attribute.
        """
        # update metadata of the output model based on the results
        # of the "blenders"
        for attr, val in self._finalize_metadata().items():
            try:
                model[attr] = val
            except KeyError:
                # Ignore keys that are in the asdf tree but not in the schema
                pass

        # patch the table into the output model and the schema
        table = self._finalize_table()
        schema = table_to_schema(table)
        model.add_schema_entry("hdrtab", schema)

        # because astropy will silently mangle boolean columns on write
        # we roundtrip the data through a BinTableHDU here to allow
        # stdatamodels to correct the data in a way that won't result in mangling
        model.hdrtab = fits_support.from_fits_hdu(fits.BinTableHDU.from_columns(table), schema)
