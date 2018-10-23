from .reference import ReferenceFileModel


__all__ = ['ThroughputModel']


class ThroughputModel(ReferenceFileModel):
    """
    A data model for filter throughput.

    Attributes
    __________
    filter_table : numpy table
         Filter throughput table
    """
    schema_url = "throughput.schema.yaml"
