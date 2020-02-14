from .reference import ReferenceFileModel


__all__ = ['PathlossModel']


class PathlossModel(ReferenceFileModel):
    """
    A data model for pathloss correction information.

    Parameters
    __________
    apertures.items.pointsource_data : numpy float32 array
         Point source pathloss

    apertures.items.pointsource_err : numpy float32 array
         Point source pathloss variance

    apertures.items.uniform_data : numpy float32 array
         Uniform source pathloss

    apertures.items.uniform_err : numpy float32 array
         Uniform source pathloss variance
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/pathloss.schema"
