import logging

from stdatamodels.jwst import datamodels as dm
from stdatamodels import filetype
from jwst.datamodels import ModelLibrary, ImageModel  # type: ignore[attr-defined]
from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.resample.resample_utils import load_custom_wcs

from . import resample
from ..stpipe import Step

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ResampleStep"]


# Force use of all DQ flagged data except for DO_NOT_USE and NON_SCIENCE
GOOD_BITS = '~DO_NOT_USE+NON_SCIENCE'


class ResampleStep(Step):
    """
    Resample input data onto a regular grid using the drizzle algorithm.

    .. note::
        When supplied via ``output_wcs``, a custom WCS overrides other custom
        WCS parameters such as ``output_shape`` (now computed from by
        ``output_wcs.bounding_box``), ``crpix``

    Parameters
    -----------
    input :  ~jwst.datamodels.JwstDataModel or ~jwst.associations.Association
        Single filename for either a single image or an association table.
    """

    class_alias = "resample"

    spec = """
        pixfrac = float(min=0.0, max=1.0, default=1.0)  # Pixel shrinkage factor
        kernel = option('square','gaussian','point','turbo','lanczos2','lanczos3',default='square')  # Flux distribution kernel
        fillval = string(default='NAN')  # Output value for pixels with no weight or flux
        weight_type = option('ivm', 'exptime', None, default='ivm')  # Input image weighting type
        output_shape = int_list(min=2, max=2, default=None)  # [x, y] order
        crpix = float_list(min=2, max=2, default=None)
        crval = float_list(min=2, max=2, default=None)
        rotation = float(default=None)  # Output image Y-axis PA relative to North
        pixel_scale_ratio = float(default=1.0)  # Ratio of input to output pixel scale
        pixel_scale = float(default=None)  # Absolute pixel scale in arcsec
        output_wcs = string(default='')  # Custom output WCS
        single = boolean(default=False)  # Resample each input to its own output grid
        blendheaders = boolean(default=True)  # Blend metadata from inputs into output
        in_memory = boolean(default=True)  # Keep images in memory
    """ # noqa: E501

    reference_file_types: list = []

    def process(self, input):

        if isinstance(input, str):
            ext = filetype.check(input)
            if ext in ("fits", "asdf"):
                input = dm.open(input)
        if isinstance(input, ModelLibrary):
            input_models = input
        elif isinstance(input, (str, dict, list)):
            input_models = ModelLibrary(input, on_disk=not self.in_memory)
        elif isinstance(input, ImageModel):
            input_models = ModelLibrary([input], on_disk=not self.in_memory)
            output = input.meta.filename
            self.blendheaders = False
        else:
            raise RuntimeError(f"Input {input} is not a 2D image.")

        try:
            output = input_models.asn["products"][0]["name"]
        except KeyError:
            # coron data goes through this path by the time it gets to
            # resampling.
            # TODO: figure out why and make sure asn_table is carried along
            output = None

        # Check that input models are 2D images
        with input_models:
            example_model = input_models.borrow(0)
            data_shape = example_model.data.shape
            input_models.shelve(example_model, 0, modify=False)
            if len(data_shape) != 2:
                # resample can only handle 2D images, not 3D cubes, etc
                raise RuntimeError(f"Input {example_model} is not a 2D image.")
            del example_model

            # Make sure all input models have consistent NaN and DO_NOT_USE values
            for model in input_models:
                match_nans_and_flags(model)
                input_models.shelve(model)
            del model

        # Setup drizzle-related parameters
        kwargs = self.get_drizpars()

        # Call the resampling routine
        if self.single:
            resamp = resample.ResampleImage(
                input_models,
                output=output,
                enable_var=False,
                compute_err="driz_err",
                **kwargs
            )
            result = resamp.resample_many_to_many(
                in_memory=self.in_memory
            )

        else:
            resamp = resample.ResampleImage(
                input_models,
                output=output,
                enable_var=True,
                compute_err="from_var",
                **kwargs
            )
            result = resamp.resample_many_to_one()

        return result

    @staticmethod
    def check_list_pars(vals, name, min_vals=None):
        """
        Validate step parameters that may take a 2-element list.

        Parameters
        ----------
        vals : list or None
            Values to validate.
        name : str
            Parameter name.
        min_vals : list, optional
            Minimum allowed values for the parameter. Must
            have 2 values.

        Returns
        -------
        values : list
            The validated list of values.

        Raises
        ------
        ValueError
            If the values do not have expected values.
        """
        if vals is None:
            return None
        if len(vals) != 2:
            raise ValueError(f"List '{name}' must have exactly two elements.")
        n = sum(x is None for x in vals)
        if n == 2:
            return None
        elif n == 0:
            if min_vals and sum(x >= y for x, y in zip(vals, min_vals)) != 2:
                raise ValueError(f"'{name}' values must be larger or equal to {list(min_vals)}")
            return list(vals)
        else:
            raise ValueError(f"Both '{name}' values must be either None or not None.")

    def get_drizpars(self):
        """
        Load all drizzle-related parameter values into kwargs list.
        """
        # Define the keys pulled from step parameters
        kwargs = dict(
            pixfrac=self.pixfrac,
            kernel=self.kernel,
            fillval=self.fillval,
            weight_type=self.weight_type,
            good_bits=GOOD_BITS,
            blendheaders=self.blendheaders,
        )

        # Custom output WCS parameters.
        output_shape = self.check_list_pars(
            self.output_shape,
            'output_shape',
            min_vals=[1, 1]
        )
        kwargs['output_wcs'] = load_custom_wcs(
            self.output_wcs,
            output_shape
        )

        wcs_pars = {
            'crpix': self.check_list_pars(self.crpix, 'crpix'),
            'crval': self.check_list_pars(self.crval, 'crval'),
            'rotation': self.rotation,
            'pixel_scale': self.pixel_scale,
            'pixel_scale_ratio': self.pixel_scale_ratio,
            'output_shape': None if output_shape is None else output_shape[::-1],
        }

        kwargs['wcs_pars'] = wcs_pars

        # Report values to processing log
        for k, v in kwargs.items():
            self.log.debug('   {}={}'.format(k, v))

        return kwargs
