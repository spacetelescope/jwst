from __future__ import absolute_import, unicode_literals, division, print_function

import warnings
from . import model_base
from astropy.modeling import models
import six

from .extension import BaseExtension
from jwst.transforms.jwextension import JWSTExtension
from gwcs.extension import GWCSExtension


jwst_extensions = [GWCSExtension(), JWSTExtension(), BaseExtension()]

__all__ = ['DistortionModel']


class DistortionModel(model_base.DataModel):
    """

    """
    schema_url = "distortion.schema.yaml"
    referencefile_schema_url ="referencefile.schema.yaml"

    """
    reffile={'author': 'ND', 'pedigree': 'ground', 'description': 'dist', 'useafter':'106-03-02'}
    meta = {"exposure": {"type": "MIR_IMAGE|MIR_FOCUS"},
            "instrument": {"name": "MIRI",
                           "filter": "F140X",
                           "channel": "12",
                           "band": "MEDIUM",
                           },
            "reffile": reffile
           }
    """
    def __init__(self, init=None, model=None, input_units="pixel", output_units='degree',
                 reffile=None, instrument=None, exp_type=None, filter=None, **kwargs):

        super(DistortionModel, self).__init__(init=init, **kwargs)
        if init is None:
            if model is None:
                raise TypeError('"model" is a required parameter.')
            if input_units is None or output_units is None:
                raise TypeError("Both input_units and output_units should be specified.")
        if init is None and reffile is None:
            warnings.warn("Reference file meta info in reffile was not provided.")
        if model is not None:
            self.model = model
        self.meta.input_units = input_units
        self.meta.output_units=output_units
        self.meta.reftype = "distortion"
        #self.meta.instrument = instrument
        #self.meta.exposure.type = exp_type
        #self.meta.instrument.filter = filter

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
        #if isinstance(path, six.string_types):
        #    self.meta.filename = os.path.basename(path)

        #self.meta.date = Time(datetime.datetime.now())
        #self.meta.date.format = 'isot'
        #self.meta.date = self.meta.date#.value
        #self.meta.model_type = self.__class__.__name__
        self.meta.model = self.model

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")
