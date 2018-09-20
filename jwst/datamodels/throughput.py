from .reference import ReferenceFileModel


__all__ = ['ThroughputModel']


class ThroughputModel(ReferenceFileModel):
    """
    A data model for filter throughput.
    """
    schema_url = "throughput.schema.yaml"
