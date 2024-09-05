from .rules import make_blender
from .schemautil import parse_schema
from .tablebuilder import TableBuilder, table_to_schema


class ModelBlender:
    def __init__(self, blend_ignore_attrs=None):
        self._first_header_meta = None
        self._blenders = None
        self._table_builder = None
        # TODO for now hard-code the polynomial_info ignore
        self._blend_ignore_attrs = ['meta.wcs', 'meta.background.polynomial_info']
        if self._blend_ignore_attrs is not None:
            self._blend_ignore_attrs.extend(blend_ignore_attrs)

    def accumulate(self, model):
        if self._first_header_meta is None:
            # capture the entire contents of the first model metadata
            self._first_header_meta = {
                attr: v
                for attr, v in model.to_flat_dict().items()
                if attr.startswith('meta') and not any((attr.startswith(i) for i in self._blend_ignore_attrs))
            }

            # search the schema for other metadata to "blend" and to add to the table
            attr_to_columns, attr_to_blend_rules = parse_schema(model.schema)

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
        header = model.to_flat_dict()

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
        model.add_schema_entry('hdrtab', table_to_schema(table))
        model.hdrtab = table
