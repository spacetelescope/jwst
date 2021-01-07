from .reference import ReferenceFileModel


__all__ = ['PsfMaskModel']


class PsfMaskModel(ReferenceFileModel):
    """
    A data model for coronagraphic 2D PSF mask reference files

    Parameters
    __________
    data : numpy float32 array
         The PSF mask
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/psfmask.schema"
