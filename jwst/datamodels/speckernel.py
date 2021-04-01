from .reference import ReferenceFileModel


__all__ = ['SpecKernelModel']


class SpecKernelModel(ReferenceFileModel):
    """
    A data model for 2D spectral kernels.

    Parameters
    __________
    wavelengths : numpy float32 array
         Wavelengths

    kernels : numpy float32 array
         Kernel values
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/speckernel.schema"
