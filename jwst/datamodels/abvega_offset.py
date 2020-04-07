from .reference import ReferenceFileModel

__all__ = ['ABVegaOffsetModel']


class ABVegaOffsetModel(ReferenceFileModel):
    """
    A data model containing offsets to convert from AB to Vega
    magnitudes.

    Parameters
    ----------
    abvega_offset : `astropy.table.Table`
        An astropy table containing offsets to convert from AB to Vega
        magnitudes.  The ``abvega_offset`` column represents m_AB -
        m_Vega.

        There are three types of tables, depending on the instrument, each
        with different column selectors.  The columns names and data types
        are:

        * FGS

          - detector: str
          - abvega_offset: float32

        * NIRCam and NIRISS

          - filter: str
          - pupil: str
          - abvega_offset: float32

        * MIRI

          - filter: str
          - abvega_offset: float32
    """

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/abvega_offset.schema"
