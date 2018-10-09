import os
from .model_base import DataModel

__all__ = ['AsnModel']


class AsnModel(DataModel):
    """
    A data model for association tables.
    """
    schema_url = "asn.schema.yaml"
    supported_formats = ['yaml', 'json', 'fits']

    def __init__(self, init=None, **kwargs):
        super(AsnModel, self).__init__(init=init, **kwargs)

        # apply logic to identify output product
        self.parse_table()

    def parse_table(self):
        self.output = None
        self.output_rootname = None
        self.inputs = []
        self.input_rootnames = []
        self.num_inputs = 0

        if not len(self.asn_table):
            return

        for i, etype in enumerate(self.asn_table.exptype):
            if 'prod' in etype.lower():
                self.output_rootname = self.asn_table.expname[i]
                for fmt in self.supported_formats:
                    fname = "{0}.{1}".format(self.output_rootname, fmt)
                    if os.path.exists(fname):
                        self.output = fname
                        break
            else:
                rootname = self.asn_table.expname[i]
                for fmt in self.supported_formats:
                    fname = "{0}.{1}".format(rootname, fmt)
                    if os.path.exists(fname):
                        self.inputs.append(fname)
                        break
                self.input_rootnames.append(rootname)
                self.num_inputs += 1
