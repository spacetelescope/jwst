from __future__ import absolute_import, unicode_literals, division, print_function

from astropy.modeling.core import Model
from astropy import units as u
from . import model_base

from .extension import BaseExtension
from jwst.transforms.jwextension import JWSTExtension
from gwcs.extension import GWCSExtension


jwst_extensions = [GWCSExtension(), JWSTExtension(), BaseExtension()]

__all__ = ['DistortionModel']


class DistortionModel(model_base.DataModel):
    """
    A model for a reference file of type "distortion".
    """
    schema_url = "distortion.schema.yaml"

    def __init__(self, init=None, model=None, input_units=None, output_units=None, **kwargs):

        super(DistortionModel, self).__init__(init=init, **kwargs)
        if init is None:
            if model is None:
                raise TypeError('"model" is a required parameter.')
            if output_units is None:
                raise TypeError("Both input_units and output_units should be specified.")
            if input_units is None:
                raise TypeError("Both input_units and output_units should be specified.")
        if model is not None:
            self.model = model
        if input_units is not None:
            self.meta.input_units = input_units
        if output_units is not None:
            self.meta.output_units = output_units
        self.meta.reftype = "distortion"

    def on_save(self, path=None):
        """
        This is a hook that is called just before saving the file.
        It can be used, for example, to update values in the metadata
        that are based on the content of the data.

        Override it in the subclass to make it do something, but don't
        forget to "chain up" to the base class, since it does things
        there, too.

        Parameters
        ----------
        path : str
            The path to the file that we're about to save to.
        """

        self.meta.model = self.model

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        assert isinstance(self.model, Model)
        assert isinstance(self.meta.input_units, (str, u.Unit))
        assert isinstance(self.meta.output_units, (str, u.Unit))
        assert self.meta.instrument.name in [NIRCAM, NIRSPEC, MIRI, TFI, FGS, NIRISS]

