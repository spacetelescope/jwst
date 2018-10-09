from .model_base import DataModel

__all__ = ['AmiLgModel']


class AmiLgModel(DataModel):
    """
    A data model for AMI LG analysis results.
    """
    schema_url = "amilg.schema.yaml"

    def get_primary_array_name(self):
        """
        Returns the name "primary" array for this model, which
        controls the size of other arrays that are implicitly created.
        This is intended to be overridden in the subclasses if the
        primary array's name is not "data".
        """
        return 'fit_image'
