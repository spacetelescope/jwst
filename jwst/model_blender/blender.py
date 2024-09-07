from astropy.io import fits
from stdatamodels import fits_support

from .rules import make_blender
from .schemautil import parse_schema
from .tablebuilder import TableBuilder, table_to_schema


__all__ = ['ModelBlender']


class ModelBlender:
    def __init__(self, blend_ignore_attrs=None):
        self._first_header_meta = None
        self._blenders = None
        self._table_builder = None
        self._blend_ignore_attrs = ['meta.wcs']
        if blend_ignore_attrs is not None:
            self._blend_ignore_attrs.extend(blend_ignore_attrs)

    def accumulate(self, model):
        if self._first_header_meta is None:
            # search the schema for other metadata to "blend" and to add to the table
            attr_to_columns, attr_to_blend_rules, schema_ignores  = parse_schema(model.schema)

            # the previous model blending code accessed every attribute
            # that will be added to the table. This causes those attributes
            # to fill in defaults (if they don't have a defined value).
            # Access all those attributes to allow the new model blender to
            # behave in the same way.
            for attr in attr_to_columns:
                model[attr]

            # update ignores list for items in schema that can't be blended
            self._blend_ignore_attrs.extend(schema_ignores)

            # capture the entire contents of the first model metadata
            self._first_header_meta = {
                attr: v
                for attr, v in model.to_flat_dict(include_arrays=False).items()
                if attr.startswith('meta') and not any((attr.startswith(i) for i in self._blend_ignore_attrs))
            }

            # make "blenders" for the metadata with special rules
            self._blenders = {
                attr: make_blender(rule)
                for attr, rule in attr_to_blend_rules.items()
                # first is implied
                if rule != 'first' and not any((attr.startswith(i) for i in self._blend_ignore_attrs))
            }

            # make a table builder using the mapping from the schema
            self._table_builder = TableBuilder(attr_to_columns)

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
        model.add_schema_entry('hdrtab', schema)

        # because astropy will silently mangle boolean columns on write
        # we roundtrip the data through a BinTableHDU here to allow
        # stdatamodels to correct the data in a way that won't result in mangling
        model.hdrtab = fits_support.from_fits_hdu(fits.BinTableHDU.from_columns(table), schema)
