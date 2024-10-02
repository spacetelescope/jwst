import json
import logging
import os
import re

from jwst.datamodels import ModelLibrary, ImageModel
import gwcs
from stcal.resample import resampled_wcs_from_models, OutputTooLargeError
from stcal.resample.utils import load_custom_wcs

from . import resample
from ..associations.asn_from_list import asn_from_list
from ..stpipe import Step


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ResampleStep", "MissingFileName"]

# Force use of all DQ flagged data except for DO_NOT_USE and NON_SCIENCE
GOOD_BITS = '~DO_NOT_USE+NON_SCIENCE'
_OUPUT_EXT = ".fits"


class MissingFileName(ValueError):
    """ Raised when in_memory is False but no output file name has been
        provided.
    """


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
        allowed_memory = float(default=None)  # Fraction of memory to use for the combined image
        in_memory = boolean(default=True)  # Keep images in memory
    """

    reference_file_types: list = []

    def process(self, input):
        output_file_name = None

        if isinstance(input, (str, dict, list)):
            input = ModelLibrary(input, on_disk=not self.in_memory)
        elif isinstance(input, ImageModel):  # TODO: do we need to support this?
            input = ModelLibrary([input], on_disk=not self.in_memory)
            # TODO: I don't get the purpose of this:
            output_file_name = input.meta.filename  # <-- ?????
            self.blendheaders = False
        elif not isinstance(input, ModelLibrary):
            raise RuntimeError(f"Input {repr(input)} is not a 2D image.")

        input_models = resample.LibModelAccess(input)

        # try to figure output file name.
        # TODO: review whether this is the intended action - not sure this
        # code reproduces what's currently in the pipeline but also,
        # not sure that code makes sense.
        if output_file_name is not None:
            output_file_name = output_file_name.strip()

        if output_file_name and output_file_name.endswith(_OUPUT_EXT):
            self._output_dir = ''
            self._output_file_name = output_file_name
        else:
            self._output_dir = output_file_name
            self._output_file_name = None

        try:
            output_file_name = input_models.asn["products"][0]["name"]
        except KeyError:
            # coron data goes through this path by the time it gets to
            # resampling.
            # TODO: figure out why and make sure asn_table is carried along
            pass

        resampled_models = []

        # Setup drizzle-related parameters
        kwargs = self.get_drizpars()

        if self.single:
            output_suffix = "outlier_i2d"

            # define output WCS if needed using
            if kwargs.pop("output_wcs", None) is None:
                output_shape = (None if self.output_shape is None
                                else self.output_shape[::-1])
                kwargs["output_wcs"], *_ = resampled_wcs_from_models(
                    input_models,
                    pixel_scale_ratio=self.pixel_scale_ratio,
                    pixel_scale=self.pixel_scale,
                    output_shape=output_shape,
                    rotation=self.rotation,
                    crpix=self.crpix,
                    crval=self.crval,
                )

            group_ids = input_models.group_indices

            # Call the resampling routine for each group of images
            for group_id in group_ids:
                input_models.set_active_group(group_id)
                log.info(f"Resampling images in group {group_id}")

                try:
                    resampler = resample.ResampleImage(
                        input_models,
                        enable_ctx=False,
                        enable_var=False,
                        **kwargs,
                    )
                    model = resampler.run()

                except OutputTooLargeError as e:
                    log.error("Not enough available memory for resample.")
                    log.error(e.msg)
                    return input

                except Exception as e:
                    log.error(
                        "The following exception occured while resampling."
                    )
                    log.error(e.msg)
                    return input

                # output file name for the resampled model:
                while resampler.input_file_names:
                    ref_file_name = resampler.input_file_names.pop()
                    if ref_file_name:
                        output_file_name = self.resampled_file_name_from_input(
                            ref_file_name,
                            suffix=output_suffix,
                        )
                        model.meta.filename = output_file_name
                        break
                else:
                    ref_file_name = None
                    if not self.in_memory:
                        raise MissingFileName(
                            "Unable to determine output file name which is "
                            "required when in_memory=False."
                        )

                if self.in_memory:
                    resampled_models.append(model)
                else:
                    # save model to file and append its file name to the output
                    # list of resampled models:
                    model.save(
                        os.path.join(self._output_dir, model.meta.filename)
                    )
                    resampled_models.append(model.meta.filename)
                    log.info(
                        f"Resampled image model saved to {model.meta.filename}"
                    )

                resampled_models.append(model.meta.filename)
                del model

        else:
            if not self.in_memory and output_file_name is None:
                raise MissingFileName(
                    "Unable to determine output file name which is "
                    "required when in_memory=False."
                )

            if output_file_name and not output_file_name.endswith(_OUPUT_EXT):
                sep = '' if output_file_name[-1] == '_' else '_'
                output_file_name = output_file_name + f"{sep}i2d{_OUPUT_EXT}"

            try:
                resampler = resample.ResampleImage(
                    input_models,
                    enable_ctx=True,
                    enable_var=True,
                    **kwargs,
                )
                model = resampler.run()
            except OutputTooLargeError as e:
                log.error("Not enough available memory for resample.")
                log.error(e.msg)
                return input
            except Exception as e:
                log.error("The following exception occured while resampling.")
                log.error(e.msg)
                return input

            model.meta.filename = output_file_name

            if self.in_memory:
                resampled_models.append(model)
            else:
                model.save(model.meta.filename)
                resampled_models.append(model.meta.filename)
                log.info(
                    f"Resampled image model saved to {model.meta.filename}"
                )

            del model

        # make a ModelLibrary obj and save it to asn if requested:
        if self.in_memory:
            return ModelLibrary(resampled_models, on_disk=False)
        else:
            # build ModelLibrary as an association from the output file names
            asn = asn_from_list(resampled_models, product_name=output_file_name)
            # serializes the asn and converts to dict
            asn_dict = json.loads(asn.dump()[1])
            return ModelLibrary(asn_dict, on_disk=True)

    def resampled_file_name_from_input(self, input_file_name, suffix):
        """ Form output file name from input image name """
        indx = input_file_name.rfind('.')
        output_type = input_file_name[indx:]
        output_root = '_'.join(
            input_file_name.replace(output_type, '').split('_')[:-1]
        )
        output_file_name = f'{output_root}_{suffix}{output_type}'
        return output_file_name

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
            wht_type=self.weight_type,
            good_bits=GOOD_BITS,
            blendheaders=self.blendheaders,
            allowed_memory=self.allowed_memory,
        )

        # Custom output WCS parameters.
        wcs_pars = {
            'output_shape': self.check_list_pars(
                self.output_shape,
                'output_shape',
                min_vals=[1, 1]
            ),
            'crpix': self.check_list_pars(self.crpix, 'crpix'),
            'crval': self.check_list_pars(self.crval, 'crval'),
            'rotation': self.rotation,
            'pixel_scale': self.pixel_scale,
            'pixel_scale_ratio': self.pixel_scale_ratio,
        }
        kwargs["wcs_pars"] = wcs_pars
        if isinstance(self.output_wcs, str):
            kwargs["output_wcs"] = load_custom_wcs(
                self.output_wcs,
                wcs_pars["output_shape"]
            )
        elif isinstance(self.output_wcs, gwcs.WCS):
            if self.output_shape is not None:
                self.output_wcs.array_shape = self.output_shape[::-1]
                self.output_wcs.pixel_shape = self.output_shape
            kwargs["output_wcs"] = self.output_wcs

        # Report values to processing log
        for k, v in kwargs.items():
            self.log.debug('   {}={}'.format(k, v))

        return kwargs

    def update_fits_wcs(self, model):
        """
        Update FITS WCS keywords of the resampled image.
        """
        # Delete any SIP-related keywords first
        pattern = r"^(cd[12]_[12]|[ab]p?_\d_\d|[ab]p?_order)$"
        regex = re.compile(pattern)

        keys = list(model.meta.wcsinfo.instance.keys())
        for key in keys:
            if regex.match(key):
                del model.meta.wcsinfo.instance[key]

        # Write new PC-matrix-based WCS based on GWCS model
        transform = model.meta.wcs.forward_transform
        model.meta.wcsinfo.crpix1 = -transform[0].offset.value + 1
        model.meta.wcsinfo.crpix2 = -transform[1].offset.value + 1
        model.meta.wcsinfo.cdelt1 = transform[3].factor.value
        model.meta.wcsinfo.cdelt2 = transform[4].factor.value
        model.meta.wcsinfo.ra_ref = transform[6].lon.value
        model.meta.wcsinfo.dec_ref = transform[6].lat.value
        model.meta.wcsinfo.crval1 = model.meta.wcsinfo.ra_ref
        model.meta.wcsinfo.crval2 = model.meta.wcsinfo.dec_ref
        model.meta.wcsinfo.pc1_1 = transform[2].matrix.value[0][0]
        model.meta.wcsinfo.pc1_2 = transform[2].matrix.value[0][1]
        model.meta.wcsinfo.pc2_1 = transform[2].matrix.value[1][0]
        model.meta.wcsinfo.pc2_2 = transform[2].matrix.value[1][1]
        model.meta.wcsinfo.ctype1 = "RA---TAN"
        model.meta.wcsinfo.ctype2 = "DEC--TAN"

        # Remove no longer relevant WCS keywords
        rm_keys = ["v2_ref", "v3_ref", "ra_ref", "dec_ref", "roll_ref",
                   "v3yangle", "vparity"]
        for key in rm_keys:
            if key in model.meta.wcsinfo.instance:
                del model.meta.wcsinfo.instance[key]
