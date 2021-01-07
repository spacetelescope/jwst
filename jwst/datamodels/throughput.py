from .reference import ReferenceFileModel


__all__ = ['ThroughputModel']


class ThroughputModel(ReferenceFileModel):
    """
    A data model for filter throughput.

    Parameters
    __________
    filter_table : numpy table
         Filter throughput table
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/throughput.schema"
