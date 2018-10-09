from .model_base import DataModel


__all__ = ['DrizProductModel']


class DrizProductModel(DataModel):
    """
    A data model for drizzle-generated products.
    """
    schema_url = "drizproduct.schema.yaml"

    @property
    def hdrtab(self):
        return self._extra_fits.HDRTAB.data

    @hdrtab.setter
    def hdrtab(self, v):
        self._extra_fits.HDRTAB.data = v
