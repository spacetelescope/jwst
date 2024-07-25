import logging
import re
from copy import deepcopy

import asdf

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer

from . import resample
from ..stpipe import Step
from ..assign_wcs import util

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
        allowed_memory = float(default=None)  # Fraction of memory to use for the combined image
        in_memory = boolean(default=True)  # Keep images in memory
    """

    reference_file_types = []

    def process(self, input):

        input = datamodels.open(input)

        if isinstance(input, ModelContainer):
            input_models = input
            try:
                output = input_models.meta.asn_table.products[0].name
            except AttributeError:
                # coron data goes through this path by the time it gets to
                # resampling.
                # TODO: figure out why and make sure asn_table is carried along
                output = None
        else:
            input_models = ModelContainer([input])
            input_models.asn_pool_name = input.meta.asn.pool_name
            input_models.asn_table_name = input.meta.asn.table_name
            output = input.meta.filename
            self.blendheaders = False

        # Check that input models are 2D images
        if len(input_models[0].data.shape) != 2:
            # resample can only handle 2D images, not 3D cubes, etc
            raise RuntimeError("Input {} is not a 2D image.".format(input_models[0]))

        # Setup drizzle-related parameters
        kwargs = self.get_drizpars()

        # Call the resampling routine
        resamp = resample.ResampleData(input_models, output=output, **kwargs)
        result = resamp.do_drizzle()

        for model in result:
            model.meta.cal_step.resample = 'COMPLETE'
            self.update_fits_wcs(model)
            util.update_s_region_imaging(model)
            model.meta.asn.pool_name = input_models.asn_pool_name
            model.meta.asn.table_name = input_models.asn_table_name

            # if pixel_scale exists, it will override pixel_scale_ratio.
            # calculate the actual value of pixel_scale_ratio based on pixel_scale
            # because source_catalog uses this value from the header.
            if self.pixel_scale is None:
                model.meta.resample.pixel_scale_ratio = self.pixel_scale_ratio
            else:
                model.meta.resample.pixel_scale_ratio = resamp.pscale_ratio
            model.meta.resample.pixfrac = kwargs['pixfrac']

        if len(result) == 1:
            result = result[0]

        input_models.close()
        return result

    @staticmethod
    def _check_list_pars(vals, name, min_vals=None):
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

    @staticmethod
    def _load_custom_wcs(asdf_wcs_file, output_shape):
        if not asdf_wcs_file:
            return None

        with asdf.open(asdf_wcs_file) as af:
            wcs = deepcopy(af.tree["wcs"])
            wcs.pixel_area = af.tree.get("pixel_area", None)
            wcs.array_shape = af.tree.get("pixel_shape", None)
            wcs.array_shape = af.tree.get("array_shape", None)

        if output_shape is not None:
            wcs.array_shape = output_shape[::-1]
            wcs.pixel_shape = output_shape
        elif wcs.pixel_shape is not None:
            wcs.array_shape = wcs.pixel_shape[::-1]
        elif wcs.array_shape is not None:
            wcs.pixel_shape = wcs.array_shape[::-1]
        elif wcs.bounding_box is not None:
            wcs.array_shape = tuple(
                int(axs[1] + 0.5)
                for axs in wcs.bounding_box.bounding_box(order="C")
            )
        else:
            raise ValueError(
                "Step argument 'output_shape' is required when custom WCS "
                "does not have neither of 'array_shape', 'pixel_shape', or "
                "'bounding_box' attributes set."
            )

        return wcs

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
            single=self.single,
            blendheaders=self.blendheaders,
            allowed_memory = self.allowed_memory,
            in_memory = self.in_memory
        )

        # Custom output WCS parameters.
        kwargs['output_shape'] = self._check_list_pars(
            self.output_shape,
            'output_shape',
            min_vals=[1, 1]
        )
        kwargs['output_wcs'] = self._load_custom_wcs(
            self.output_wcs,
            kwargs['output_shape']
        )
        kwargs['crpix'] = self._check_list_pars(self.crpix, 'crpix')
        kwargs['crval'] = self._check_list_pars(self.crval, 'crval')
        kwargs['rotation'] = self.rotation
        kwargs['pscale'] = self.pixel_scale
        kwargs['pscale_ratio'] = self.pixel_scale_ratio

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
        rm_keys = ['v2_ref', 'v3_ref', 'ra_ref', 'dec_ref', 'roll_ref',
                   'v3yangle', 'vparity']
        for key in rm_keys:
            if key in model.meta.wcsinfo.instance:
                del model.meta.wcsinfo.instance[key]
