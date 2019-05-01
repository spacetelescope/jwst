from .reference import ReferenceFileModel


__all__ = ['RSCDModel']


class RSCDModel(ReferenceFileModel):
    """
    A data model for the RSCD reference file.

    Parameters
    __________
    rscd_table : numpy table
        Reference file for RSCD correction
        A table with seven columns, three string-valued that identify which
        row to select, and four float columns containing coefficients.
    """
    schema_url = "rscd.schema.yaml"
