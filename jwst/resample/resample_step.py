import logging

from stdatamodels.jwst import datamodels as dm
from stdatamodels import filetype
from jwst.datamodels import ModelLibrary, ImageModel  # type: ignore[attr-defined]
from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.resample.resample_utils import load_custom_wcs

from . import resample
from jwst.stpipe import Step

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ResampleStep"]


# Force use of all DQ flagged data except for DO_NOT_USE and NON_SCIENCE
GOOD_BITS = "~DO_NOT_USE+NON_SCIENCE"


class ResampleStep(Step):
    """
    Resample imaging data onto a regular grid using the drizzle algorithm.

    .. note::
        When supplied via ``output_wcs``, a custom WCS overrides other custom
        WCS parameters such as ``output_shape`` (now computed from by
        ``output_wcs.bounding_box``), ``crpix``
    """

    class_alias = "resample"

    spec = """
        pixfrac = float(min=0.0, max=1.0, default=1.0)  # Pixel shrinkage factor
        kernel = option('square','gaussian','point','turbo','lanczos2','lanczos3',default='square')  # Flux distribution kernel
        fillval = string(default='NAN')  # Output value for pixels with no weight or flux
        weight_type = option('ivm', 'exptime', None, default='ivm')  # Input image weighting type
        output_shape = int_list(min=2, max=2, default=None)  # [x, y] order
        crpix = float_list(min=2, max=2, default=None)  # 0-based image coordinates of the reference pixel
        crval = float_list(min=2, max=2, default=None)  # world coordinates of the reference pixel
        rotation = float(default=None)  # Output image Y-axis PA relative to North
        pixel_scale_ratio = float(default=1.0)  # Ratio of output to input pixel scale.
        pixel_scale = float(default=None)  # Absolute pixel scale in arcsec
        output_wcs = string(default='')  # Custom output WCS
        single = boolean(default=False)  # Resample each input to its own output grid
        blendheaders = boolean(default=True)  # Blend metadata from inputs into output
        in_memory = boolean(default=True)  # Keep images in memory
        enable_ctx = boolean(default=True)  # Compute and report the context array
        enable_err = boolean(default=True)  # Compute and report the err array
        report_var = boolean(default=True)  # Report the variance array
    """  # noqa: E501

    reference_file_types: list = []

    def process(self, input_data):
        """
        Run the resample step on the input data.

        Parameters
        ----------
        input_data : str, ImageModel, or any asn-type input loadable into ModelLibrary
            Filename pointing to an ImageModel or an association, or the ImageModel or
            association itself.

        Returns
        -------
        ModelLibrary or ImageModel
            The final output data. If the `single` parameter is set to True, then this
            is a single ImageModel; otherwise, it is a ModelLibrary.
        """
        if isinstance(input_data, str):
            ext = filetype.check(input_data)
            if ext in ("fits", "asdf"):
                input_data = dm.open(input_data)
        if isinstance(input_data, ModelLibrary):
            input_models = input_data
        elif isinstance(input_data, (str, dict, list)):
            input_models = ModelLibrary(input_data, on_disk=not self.in_memory)
        elif isinstance(input_data, ImageModel):
            input_models = ModelLibrary([input_data], on_disk=not self.in_memory)
            output = input_data.meta.filename
            self.blendheaders = False
        else:
            raise TypeError(f"Input {input_data} is not a 2D image.")

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
                input_models, output=output, enable_var=False, compute_err="driz_err", **kwargs
            )
            result = resamp.resample_many_to_many(in_memory=self.in_memory)

        else:
            if self.enable_err:
                # If error is enabled, we compute the error from the variance
                compute_err = "from_var"
                enable_var = True
                report_var = self.report_var
            else:
                # otherwise do not compute the error arrays at all
                enable_var = False
                compute_err = None
                report_var = False
            resamp = resample.ResampleImage(
                input_models,
                output=output,
                enable_ctx=self.enable_ctx,
                enable_var=enable_var,
                report_var=report_var,
                compute_err=compute_err,
                **kwargs,
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
            if min_vals and sum(x >= y for x, y in zip(vals, min_vals, strict=True)) != 2:
                raise ValueError(f"'{name}' values must be larger or equal to {list(min_vals)}")
            return list(vals)
        else:
            raise ValueError(f"Both '{name}' values must be either None or not None.")

    def get_drizpars(self):
        """
        Load all drizzle-related parameter values into kwargs list.

        Returns
        -------
        kwargs : dict
            Dictionary of drizzle parameters
        """
        # Define the keys pulled from step parameters
        kwargs = {
            "pixfrac": self.pixfrac,
            "kernel": self.kernel,
            "fillval": self.fillval,
            "weight_type": self.weight_type,
            "good_bits": GOOD_BITS,
            "blendheaders": self.blendheaders,
        }

        # Custom output WCS parameters.
        output_shape = self.check_list_pars(self.output_shape, "output_shape", min_vals=[1, 1])
        kwargs["output_wcs"] = load_custom_wcs(self.output_wcs, output_shape)

        wcs_pars = {
            "crpix": self.check_list_pars(self.crpix, "crpix"),
            "crval": self.check_list_pars(self.crval, "crval"),
            "rotation": self.rotation,
            "pixel_scale": self.pixel_scale,
            "pixel_scale_ratio": self.pixel_scale_ratio,
            "output_shape": None if output_shape is None else output_shape[::-1],
        }

        kwargs["wcs_pars"] = wcs_pars

        # Report values to processing log
        for k, v in kwargs.items():
            self.log.debug(f"   {k}={v}")

        return kwargs
