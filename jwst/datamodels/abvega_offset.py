import warnings

from stdatamodels.validate import ValidationWarning
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

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/abvegaoffset.schema"

    def validate(self):
        super(ABVegaOffsetModel, self).validate()
        try:
            assert len(self.abvega_offset) > 0
            assert (self.meta.instrument.name in
                    ('FGS', 'MIRI', 'NIRCAM', 'NIRISS'))
            assert 'abvega_offset' in self.abvega_offset.colnames

            if self.meta.instrument.name == 'FGS':
                assert 'detector' in self.abvega_offset.colnames

            if self.meta.instrument.name == 'MIRI':
                assert 'filter' in self.abvega_offset.colnames

            if self.meta.instrument.name in ('NIRCAM', 'NIRISS'):
                assert 'filter' in self.abvega_offset.colnames
                assert 'pupil' in self.abvega_offset.colnames

        except AssertionError as errmsg:
            if self._strict_validation:
                raise AssertionError(errmsg)
            else:
                warnings.warn(str(errmsg), ValidationWarning)
