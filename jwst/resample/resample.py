from copy import deepcopy
import logging
import json
import math
import os
import re
import warnings

import numpy as np

from drizzle.resample import Drizzle
from spherical_geometry.polygon import SphericalPolygon

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels.dqflags import pixel

from stcal.resample import (
    Resample,
    OutputTooLargeError,
)
from stcal.resample.utils import (
    compute_wcs_pixel_area,
    is_imaging_wcs,
    resample_range,
)

from jwst.datamodels import ModelLibrary
from jwst.associations.asn_from_list import asn_from_list

from jwst.model_blender.blender import ModelBlender
from jwst.resample import resample_utils
from jwst.assign_wcs import util as assign_wcs_util


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = [
    "ResampleImage",
    "OutputTooLargeError",
    "is_imaging_wcs",
]

_SUPPORTED_CUSTOM_WCS_PARS = [
    'pixel_scale_ratio',
    'pixel_scale',
    'output_shape',
    'crpix',
    'crval',
    'rotation',
]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ResampleImage(Resample):
    dq_flag_name_map = pixel

    def __init__(self, input_models, pixfrac=1.0, kernel="square",
                 fillval="NAN", wht_type="ivm", good_bits=0,
                 blendheaders=True, output_wcs=None, wcs_pars=None,
                 output=None, enable_ctx=True, enable_var=True,
                 compute_err=None, allowed_memory=None, asn_id=None,
                 in_memory=True):
        """
        Parameters
        ----------
        input_models : LibModelAccessBase
            A `LibModelAccessBase`-based object allowing iterating over
            all contained models of interest.

        pixfrac : float, optional
            The fraction of a pixel that the pixel flux is confined to. The
            default value of 1 has the pixel flux evenly spread across the
            image. A value of 0.5 confines it to half a pixel in the linear
            dimension, so the flux is confined to a quarter of the pixel area
            when the square kernel is used.

        kernel: {"square", "gaussian", "point", "turbo", "lanczos2", "lanczos3"}, optional
            The name of the kernel used to combine the input. The choice of
            kernel controls the distribution of flux over the kernel.
            The square kernel is the default.

            .. warning::
               The "gaussian" and "lanczos2/3" kernels **DO NOT**
               conserve flux.

        fillval: float, None, str, optional
            The value of output pixels that did not have contributions from
            input images' pixels. When ``fillval`` is either `None` or
            ``"INDEF"`` and ``out_img`` is provided, the values of ``out_img``
            will not be modified. When ``fillval`` is either `None` or
            ``"INDEF"`` and ``out_img`` is **not provided**, the values of
            ``out_img`` will be initialized to `numpy.nan`. If ``fillval``
            is a string that can be converted to a number, then the output
            pixels with no contributions from input images will be set to this
            ``fillval`` value.

        wht_type : {"exptime", "ivm"}, optional
            The weighting type for adding models' data. For ``wht_type="ivm"``
            (the default), the weighting will be determined per-pixel using
            the inverse of the read noise (VAR_RNOISE) array stored in each
            input image. If the ``VAR_RNOISE`` array does not exist,
            the variance is set to 1 for all pixels (i.e., equal weighting).
            If ``weight_type="exptime"``, the weight will be set equal
            to the measurement time (``TMEASURE``) when available and to
            the exposure time (``EFFEXPTM``) otherwise.

        good_bits : int, str, None, optional
            An integer bit mask, `None`, a Python list of bit flags, a comma-,
            or ``'|'``-separated, ``'+'``-separated string list of integer
            bit flags or mnemonic flag names that indicate what bits in models'
            DQ bitfield array should be *ignored* (i.e., zeroed).

            When co-adding models using :py:meth:`add_model`, any pixels with
            a non-zero DQ values are assigned a weight of zero and therefore
            they do not contribute to the output (resampled) data.
            ``good_bits`` provides a mean to ignore some of the DQ bitflags.

            When ``good_bits`` is an integer, it must be
            the sum of all the DQ bit values from the input model's
            DQ array that should be considered "good" (or ignored). For
            example, if pixels in the DQ array can be
            combinations of 1, 2, 4, and 8 flags and one wants to consider DQ
            "defects" having flags 2 and 4 as being acceptable, then
            ``good_bits`` should be set to 2+4=6. Then a pixel with DQ values
            2,4, or 6 will be considered a good pixel, while a pixel with
            DQ value, e.g., 1+2=3, 4+8=12, etc. will be flagged as
            a "bad" pixel.

            Alternatively, when ``good_bits`` is a string, it can be a
            comma-separated or '+' separated list of integer bit flags that
            should be summed to obtain the final "good" bits. For example,
            both "4,8" and "4+8" are equivalent to integer ``good_bits=12``.

            Finally, instead of integers, ``good_bits`` can be a string of
            comma-separated mnemonics. For example, for JWST, all the following
            specifications are equivalent:

            `"12" == "4+8" == "4, 8" == "JUMP_DET, DROPOUT"`

            In order to "translate" mnemonic code to integer bit flags,
            ``Resample.dq_flag_name_map`` attribute must be set to either
            a dictionary (with keys being mnemonc codes and the values being
            integer flags) or a `~astropy.nddata.BitFlagNameMap`.

            In order to reverse the meaning of the flags
            from indicating values of the "good" DQ flags
            to indicating the "bad" DQ flags, prepend '~' to the string
            value. For example, in order to exclude pixels with
            DQ flags 4 and 8 for computations and to consider
            as "good" all other pixels (regardless of their DQ flag),
            use a value of ``~4+8``, or ``~4,8``. A string value of
            ``~0`` would be equivalent to a setting of ``None``.

            Default value (0) will make *all* pixels with non-zero DQ
            values be considered "bad" pixels, and the corresponding data
            pixels will be assigned zero weight and thus these pixels
            will not contribute to the output resampled data array.

            Set `good_bits` to `None` to turn off the use of model's DQ
            array.

            For more details, see documentation for
            `astropy.nddata.bitmask.extend_bit_flag_map`.

        blendheaders : bool, optional
            Indicates whether to blend metadata from all input models and
            store the combined result to the output model.

        output_wcs : dict, WCS object, None, optional
            Specifies output WCS either directly as a WCS or a dictionary
            with keys ``'wcs'`` (WCS object) and ``'pixel_scale'``
            (pixel scale in arcseconds). ``'pixel_scale'``, when provided,
            will be used for computation of drizzle scaling factor. When it is
            not provided, output pixel scale will be *estimated* from the
            provided WCS object. ``output_wcs`` object is required when
            ``output_model`` is `None`. ``output_wcs`` is ignored when
            ``output_model`` is provided.

        wcs_pars : dict, None, optional
            A dictionary of custom WCS parameters used to define an
            output WCS from input models' outlines. This argument is ignored
            when ``output_wcs`` is specified.

            List of supported parameters (keywords in the dictionary):

                - ``pixel_scale_ratio`` : float

                    Desired pixel scale ratio defined as the ratio of the
                    desired output pixel scale to the first input model's pixel
                    scale computed from this model's WCS at the fiducial point
                    (taken as the ``ref_ra`` and ``ref_dec`` from the
                    ``wcsinfo`` meta attribute of the first input image).
                    Ignored when ``pixel_scale`` is specified. Default value
                    is ``1.0``.

                - ``pixel_scale`` : float, None

                    Desired pixel scale (in degrees) of the output WCS. When
                    provided, overrides ``pixel_scale_ratio``. Default value
                    is `None`.

                - ``output_shape`` : tuple of two integers (int, int), None

                    Shape of the image (data array) using ``np.ndarray``
                    convention (``ny`` first and ``nx`` second). This value
                    will be assigned to ``pixel_shape`` and ``array_shape``
                    properties of the returned WCS object. Default value is
                    `None`.

                - ``rotation`` : float, None

                    Position angle of output image's Y-axis relative to North.
                    A value of ``0.0`` would orient the final output image to
                    be North up. The default of `None` specifies that the
                    images will not be rotated, but will instead be resampled
                    in the default orientation for the camera with the x and y
                    axes of the resampled image corresponding approximately
                    to the detector axes. Ignored when ``transform`` is
                    provided. Default value is `None`.

                - ``crpix`` : tuple of float, None

                    Position of the reference pixel in the resampled image
                    array. If ``crpix`` is not specified, it will be set to
                    the center of the bounding box of the returned WCS object.
                    Default value is `None`.

                - ``crval`` : tuple of float, None

                    Right ascension and declination of the reference pixel.
                    Automatically computed if not provided. Default value is
                    `None`.

        output : str, None, optional
            Filename for the output model.

        accumulate : bool, optional
            Indicates whether resampled models should be added to the
            provided ``output_model`` data or if new arrays should be
            created.

        enable_ctx : bool, optional
            Indicates whether to create a context image. If ``disable_ctx``
            is set to `True`, parameters ``out_ctx``, ``begin_ctx_id``, and
            ``max_ctx_id`` will be ignored.

        enable_var : bool, optional
            Indicates whether to resample variance arrays.

        compute_err : {"from_var", "driz_err"}, None, optional
            - ``"from_var"``: compute output model's error array from
              all (Poisson, flat, readout) resampled variance arrays.
              Setting ``compute_err`` to ``"from_var"`` will assume
              ``enable_var`` was set to `True` regardless of actual
              value of the parameter ``enable_var``.

            - ``"driz_err"``: compute output model's error array by drizzling
              together all input models' error arrays.

            Error array will be assigned to ``'err'`` key of the output model.

            .. note::
                At this time, output error array is not equivalent to
                error propagation results.

        allowed_memory : float, None
            Fraction of memory allowed to be used for resampling. If
            ``allowed_memory`` is `None` then no check for available memory
            will be performed.

        in_memory : bool, optional

        asn_id : str, None, optional

        """
        self.input_models = input_models
        self.output_jwst_model = None

        self.output_dir = None
        self.output_filename = output
        if output is not None and '.fits' not in str(output):
            self.output_dir = output
            self.output_filename = None
        self.intermediate_suffix = 'outlier_i2d'

        self.blendheaders = blendheaders
        if blendheaders:
            self._blender = ModelBlender(
                blend_ignore_attrs=[
                    'meta.photometry.pixelarea_steradians',
                    'meta.photometry.pixelarea_arcsecsq',
                    'meta.filename',
                ]
            )

        self.in_memory = in_memory
        self.asn_id = asn_id

        # check wcs_pars has supported keywords:
        if wcs_pars is None:
            wcs_pars = {}
        elif wcs_pars:
            unsup = []
            unsup = set(wcs_pars.keys()).difference(_SUPPORTED_CUSTOM_WCS_PARS)
            if unsup:
                raise KeyError(
                    "Unsupported custom WCS parameters: "
                    f"{','.join(map(repr, unsup))}."
                )

        # determine output WCS:
        shape = wcs_pars.get("output_shape")

        if output_wcs is None:
            if (pscale := wcs_pars.get("pixel_scale")) is not None:
                pscale /= 3600.0
            wcs, _, ps, ps_ratio = resample_utils.resampled_wcs_from_models(
                input_models,
                pixel_scale_ratio=wcs_pars.get("pixel_scale_ratio", 1.0),
                pixel_scale=pscale,
                output_shape=shape,
                rotation=wcs_pars.get("rotation"),
                crpix=wcs_pars.get("crpix"),
                crval=wcs_pars.get("crval"),
            )

            output_wcs = {
                "wcs": wcs,
                "pixel_scale": 3600.0 * ps,
                "pixel_scale_ratio": ps_ratio,
            }

        else:
            if shape is None:
                if output_wcs.array_shape is None:
                    raise ValueError(
                        "Custom WCS objects must have the 'array_shape' "
                        "attribute set (defined)."
                    )
            else:
                output_wcs = deepcopy(output_wcs)
                output_wcs.array_shape = shape

        super().__init__(
            n_input_models=len(input_models),
            pixfrac=pixfrac,
            kernel=kernel,
            fillval=fillval,
            wht_type=wht_type,
            good_bits=good_bits,
            output_wcs=output_wcs,
            output_model=None,
            accumulate=False,
            enable_ctx=enable_ctx,
            enable_var=enable_var,
            compute_err=compute_err,
            allowed_memory=allowed_memory,
        )

    def _input_model_to_dict(self, model):
        # wcs = model.meta.wcs

        model_dict = {
            # arrays:
            "data": model.data,
            "dq": model.dq,

            # meta:
            "filename": model.meta.filename,
            "group_id": model.meta.group_id,
            "s_region": model.meta.wcsinfo.s_region,
            "wcs": model.meta.wcs,
            "wcsinfo": model.meta.wcsinfo,

            "exposure_time": model.meta.exposure.exposure_time,
            "start_time": model.meta.exposure.start_time,
            "end_time": model.meta.exposure.end_time,
            "duration": model.meta.exposure.duration,
            "measurement_time": model.meta.exposure.measurement_time,
            "effective_exposure_time": model.meta.exposure.effective_exposure_time,
            "elapsed_exposure_time": model.meta.exposure.elapsed_exposure_time,

            "pixelarea_steradians": model.meta.photometry.pixelarea_steradians,
            "pixelarea_arcsecsq": model.meta.photometry.pixelarea_arcsecsq,

            "level": model.meta.background.level,  # sky level
            "subtracted": model.meta.background.subtracted,

            # spectroscopy-specific:
            "instrument_name": model.meta.instrument.name,
            "exposure_type": model.meta.exposure.type,
        }

        if self._enable_var:
            model_dict["var_flat"] = model.var_flat
            model_dict["var_rnoise"] = model.var_rnoise
            model_dict["var_poisson"] = model.var_poisson

        elif (self.weight_type is not None and
                self.weight_type.startswith('ivm')):
            model_dict["var_rnoise"] = model.var_rnoise

        if self._compute_err == "driz_err":
            model_dict["err"] = model.err

        return model_dict

    def _create_output_jwst_model(self, ref_input_model=None):
        """ Create a new blank model and update it's meta with info from ``ref_input_model``. """
        output_model = datamodels.ImageModel(None)  # tuple(self.output_wcs.array_shape))

        # update meta data and wcs
        if ref_input_model is not None:
            output_model.update(ref_input_model)
        output_model.meta.wcs = self.output_wcs
        return output_model

    def _update_output_model(self, model, info_dict):
        model.data = info_dict["data"]
        model.wht = info_dict["wht"]
        if self._enable_ctx:
            model.con = info_dict["con"]
        if self._compute_err:
            model.err = info_dict["err"]
        if self._enable_var:
            model.var_rnoise = info_dict["var_rnoise"]
            model.var_flat = info_dict["var_flat"]
            model.var_poisson = info_dict["var_poisson"]

        model.meta.wcs = info_dict["wcs"]
        model.meta.photometry.pixelarea_steradians = info_dict["pixelarea_steradians"]
        model.meta.photometry.pixelarea_arcsecsq = info_dict["pixelarea_arcsecsq"]

        model.meta.resample.pointings = info_dict["pointings"]
        # model.meta.resample.n_coadds = info_dict["n_coadds"]

        model.meta.resample.pixel_scale_ratio = info_dict["pixel_scale_ratio"]
        model.meta.resample.pixfrac = info_dict["pixfrac"]
        model.meta.resample.kernel = info_dict["kernel"]
        model.meta.resample.fillval = info_dict["fillval"]
        model.meta.resample.weight_type = info_dict["weight_type"]

        model.meta.exposure.exposure_time = info_dict["exposure_time"]
        model.meta.exposure.start_time = info_dict["start_time"]
        model.meta.exposure.end_time = info_dict["end_time"]
        model.meta.exposure.duration = info_dict["duration"]
        model.meta.exposure.measurement_time = info_dict["measurement_time"]
        model.meta.exposure.effective_exposure_time = info_dict["effective_exposure_time"]
        model.meta.exposure.elapsed_exposure_time = info_dict["elapsed_exposure_time"]

        model.meta.cal_step.resample = 'COMPLETE'

    def add_model(self, model):
        """ Resamples model image and either variance data (if ``enable_var``
        was `True`) or error data (if ``enable_err`` was `True`) and adds
        them using appropriate weighting to the corresponding
        arrays of the output model. It also updates resampled data weight,
        the context array (if ``enable_ctx`` is `True`), relevant output
        model's values such as "n_coadds".

        Whenever ``model`` has a unique group ID that was never processed
        before, the "pointings" value of the output model is incremented and
        the "group_id" attribute is updated. Also, time counters are updated
        with new values from the input ``model`` by calling
        :py:meth:`~Resample.update_time`.

        Parameters
        ----------
        model : dict
            A dictionary containing data arrays and other meta attributes
            and values of actual models used by pipelines.

        """
        super().add_model(self._input_model_to_dict(model))
        if self.output_jwst_model is None:
            self.output_jwst_model = self._create_output_jwst_model(
                ref_input_model=model
            )
        if self.blendheaders:
            self._blender.accumulate(model)

    def finalize(self, free_memory=True):
        """ Finalizes all computations and frees temporary objects.

        ``finalize`` calls :py:meth:`~Resample.finalize_resample_variance` and
        :py:meth:`~Resample.finalize_time_info`.

        .. warning::
          If ``enable_var=True`` and :py:meth:`~Resample.finalize` is called
          with ``free_memory=True`` then intermediate arrays holding variance
          weights will be lost and so continuing adding new models after
          a call to :py:meth:`~Resample.finalize` will result in incorrect
          variance.

        """
        if self.blendheaders:
            self._blender.finalize_model(self.output_jwst_model)
        super().finalize(free_memory=True)

        self._update_output_model(
            self.output_jwst_model,
            self.output_model,
        )

        if is_imaging_wcs(self.output_jwst_model.meta.wcs):
            # only for an imaging WCS:
            self.update_fits_wcsinfo(self.output_jwst_model)
            assign_wcs_util.update_s_region_imaging(self.output_jwst_model)
        else:
            assign_wcs_util.update_s_region_spectral(self.output_jwst_model)

        self.output_jwst_model.meta.cal_step.resample = 'COMPLETE'

    def reset_arrays(self, reset_output=True, n_input_models=None):
        """ Initialize/reset `Drizzle` objects, `ModelBlender`, output model
        and arrays, and time counters. Output WCS and shape are not modified
        from `Resample` object initialization. This method needs to be called
        before calling :py:meth:`add_model` for the first time if
        :py:meth:`finalize` was previously called.

        Parameters
        ----------
        reset_output : bool, optional
            When `True` a new output model will be created. Otherwise new
            models will be resampled and added to existing output data arrays.

        n_input_models : int, None, optional
            Number of input models expected to be resampled. When provided,
            this is used to estimate memory requirements and optimize memory
            allocation for the context array.

        """
        super().reset_arrays(
            reset_output=reset_output,
            n_input_models=n_input_models
        )
        if self.blendheaders:
            self._blender = ModelBlender(
                blend_ignore_attrs=[
                    'meta.photometry.pixelarea_steradians',
                    'meta.photometry.pixelarea_arcsecsq',
                    'meta.filename',
                ]
            )
        self.output_jwst_model = None

    def _create_output_model(self, ref_input_model=None):
        """ Create a new blank model and update it's meta with info from
        ``ref_input_model``.
        """
        output_model = datamodels.ImageModel(None)

        # update meta data and wcs
        if ref_input_model is not None:
            output_model.update(ref_input_model)
        output_model.meta.wcs = self._output_wcs

        return output_model

    def resample_group(self, indices, compute_error=False):
        """ Resample multiple input images that belong to a single
        ``group_id`` as specified by ``indices``.

        Parameters
        ----------
        indices : list

        compute_error : bool, optional
            If True, an approximate error image will be resampled
            alongside the science image.
        """
        if self.output_jwst_model is not None:
            self.reset_arrays(reset_output=True)

        output_model_filename = ''

        log.info(f"{len(indices)} exposures to drizzle together")
        for index in indices:
            model = self.input_models.borrow(index)
            if self.output_jwst_model is None:
                # Determine output file type from input exposure filenames
                # Use this for defining the output filename
                indx = model.meta.filename.rfind('.')
                output_type = model.meta.filename[indx:]
                output_root = '_'.join(model.meta.filename.replace(
                    output_type,
                    ''
                ).split('_')[:-1])
                output_model_filename = (
                    f'{output_root}_'
                    f'{self.intermediate_suffix}{output_type}'
                )

            if isinstance(model, datamodels.SlitModel):
                # must call this explicitly to populate area extension
                # although the existence of this extension may not be necessary
                model.area = model.area

            self.add_model(model)
            self.input_models.shelve(model, index, modify=False)
            del model

        self.finalize()
        copy_asn_info_from_library(self.input_models, self.output_jwst_model)
        self.output_jwst_model.meta.filename = output_model_filename
        return self.output_jwst_model

    def resample_many_to_many(self):
        """Resample many inputs to many outputs where outputs have a common frame.

        Coadd only different detectors of the same exposure, i.e. map NRCA5 and
        NRCB5 onto the same output image, as they image different areas of the
        sky.

        Used for outlier detection
        """
        output_models = []
        for group_id, indices in self.input_models.group_indices.items():

            output_model = self.resample_group(self.input_models, indices)

            if not self.in_memory:
                # Write out model to disk, then return filename
                output_name = output_model.meta.filename
                if self.output_dir is not None:
                    output_name = os.path.join(self.output_dir, output_name)
                output_model.save(output_name)
                log.info(f"Saved model in {output_name}")
                output_models.append(output_name)
            else:
                output_models.append(output_model)

        if self.in_memory:
            # build ModelLibrary as a list of in-memory models
            return ModelLibrary(output_models, on_disk=False)
        else:
            # build ModelLibrary as an association from the output files
            # this saves memory if there are multiple groups
            asn = asn_from_list(output_models, product_name='outlier_i2d')
            asn_dict = json.loads(asn.dump()[1])  # serializes the asn and converts to dict
            return ModelLibrary(asn_dict, on_disk=True)

    def resample_many_to_one(self):
        """Resample and coadd many inputs to a single output.

        Used for stage 3 resampling
        """
        # if self.output_jwst_model is not None:
        #     self.reset_arrays(reset_output=True)

        log.info("Resampling science and variance data")

        with self.input_models:
            for model in self.input_models:
                self.add_model(model)
                self.input_models.shelve(model)

        self.finalize()
        self.output_jwst_model.meta.filename = self.output_filename
        copy_asn_info_from_library(self.input_models, self.output_jwst_model)

        return self.output_jwst_model

    @staticmethod
    def update_fits_wcsinfo(model):
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


class ResampleData:
    """
    This is the controlling routine for the resampling process.

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input model, such as pixfrac,
         weight type, exposure time (if relevant), and kernel, and merges
         them with any user-provided values.
      2. Creates output WCS based on input images and define mapping function
         between all input arrays and the output array.
      3. Updates output data model with output arrays from drizzle, including
         a record of metadata from all input models.
    """

    def __init__(self, input_models, output=None, single=False, blendheaders=True,
                 pixfrac=1.0, kernel="square", fillval="NAN", wht_type="ivm",
                 good_bits=0, pscale_ratio=1.0, pscale=None, **kwargs):
        """
        Parameters
        ----------
        input_models : library of objects
            library of data models, one for each input image

        output : str
            filename for output

        kwargs : dict
            Other parameters.

            .. note::
                ``output_shape`` is in the ``x, y`` order.

            .. note::
                ``in_memory`` controls whether or not the resampled
                array from ``resample_many_to_many()``
                should be kept in memory or written out to disk and
                deleted from memory. Default value is `True` to keep
                all products in memory.
        """
        self.output_dir = None
        self.output_filename = output
        if output is not None and '.fits' not in str(output):
            self.output_dir = output
            self.output_filename = None
        self.intermediate_suffix = 'outlier_i2d'

        self.pscale_ratio = pscale_ratio
        self.single = single
        self.blendheaders = blendheaders
        self.pixfrac = pixfrac
        self.kernel = kernel
        self.fillval = fillval
        self.weight_type = wht_type
        self.good_bits = good_bits
        self.in_memory = kwargs.get('in_memory', True)
        self.input_pixscale0 = None  # computed pixel scale of the first image (deg)
        self._recalc_pscale_ratio = pscale is not None

        log.info(f"Driz parameter kernel: {self.kernel}")
        log.info(f"Driz parameter pixfrac: {self.pixfrac}")
        log.info(f"Driz parameter fillval: {self.fillval}")
        log.info(f"Driz parameter weight_type: {self.weight_type}")

        # output_wcs = kwargs["wcs_pars"].get('output_wcs', None)
        # output_shape = kwargs["wcs_pars"].get('output_shape', None)
        # crpix = kwargs["wcs_pars"].get('crpix', None)
        # crval = kwargs["wcs_pars"].get('crval', None)
        # rotation = kwargs["wcs_pars"].get('rotation', None)

        output_wcs = kwargs.get('output_wcs', None)
        output_shape = kwargs.get('output_shape', None)
        crpix = kwargs.get('crpix', None)
        crval = kwargs.get('crval', None)
        rotation = kwargs.get('rotation', None)

        self.asn_id = kwargs.get('asn_id', None)

        if pscale is not None:
            log.info(f'Output pixel scale: {pscale} arcsec.')
            pscale /= 3600.0
        else:
            log.info(f'Output pixel scale ratio: {pscale_ratio}')

        if output_wcs:
            # Use user-supplied reference WCS for the resampled image:
            self.output_wcs = output_wcs
            if output_shape is not None:
                self.output_wcs.array_shape = output_shape[::-1]

            if output_wcs.pixel_area is None:
                self.output_pix_area = compute_wcs_pixel_area(self.output_wcs)
                if self.output_pix_area is None:
                    raise ValueError(
                        "Unable to compute output pixel area from 'output_wcs'."
                    )
            else:
                self.output_pix_area = output_wcs.pixel_area

        else:
            # Define output WCS based on all inputs, including a reference WCS:
            self.output_wcs = resample_utils.make_output_wcs(
                input_models,
                ref_wcs=output_wcs,
                pscale_ratio=self.pscale_ratio,
                pscale=pscale,
                rotation=rotation,
                shape=None if output_shape is None else output_shape[::-1],
                crpix=crpix,
                crval=crval
            )
            # Estimate output pixel area in Sr. NOTE: in principle we could
            # use the same algorithm as for when output_wcs is provided by the
            # user.
            tr = self.output_wcs.pipeline[0].transform
            self.output_pix_area = (
                np.deg2rad(tr['cdelt1'].factor.value) *
                np.deg2rad(tr['cdelt2'].factor.value)
            )

        self.output_array_shape = tuple(self.output_wcs.array_shape)
        log.debug(f"Output mosaic size: {self.output_array_shape}")

        if pscale is None:
            pscale = np.rad2deg(math.sqrt(self.output_pix_area))
            log.info(f'Computed output pixel scale: {3600.0 * pscale} arcsec.')

        self.pscale = pscale  # in deg

        log.debug(f"Output mosaic size: {tuple(self.output_wcs.pixel_shape)}")


    def _create_output_model(self, ref_input_model=None):
        """ Create a new blank model and update it's meta with info from ``ref_input_model``. """
        output_model = datamodels.ImageModel(None)  # tuple(self.output_wcs.array_shape))

        # update meta data and wcs
        if ref_input_model is not None:
            output_model.update(ref_input_model)
        output_model.meta.wcs = self.output_wcs
        if self.output_pix_area is not None:
            output_model.meta.photometry.pixelarea_steradians = self.output_pix_area
            output_model.meta.photometry.pixelarea_arcsecsq = (
                self.output_pix_area * np.rad2deg(3600)**2
            )
        return output_model

    def do_drizzle(self, input_models):
        """Pick the correct drizzling mode based on self.single
        """
        if self.single:
            return self.resample_many_to_many(input_models)
        else:
            return self.resample_many_to_one(input_models)

    def _get_intensity_scale(self, img):
        """
        Compute an intensity scale from the input and output pixel area.

        For imaging data, the scaling is used to account for differences
        between the nominal pixel area and the average pixel area for
        the input data.

        For spectral data, the scaling is used to account for flux
        conservation with non-unity pixel scale ratios, when the
        data units are flux density.

        Parameters
        ----------
        img : DataModel
            The input data model.

        Returns
        -------
        iscale : float
            The scale to apply to the input data before drizzling.
        """
        input_pixflux_area = img.meta.photometry.pixelarea_steradians
        if input_pixflux_area:
            if 'SPECTRAL' in img.meta.wcs.output_frame.axes_type:
                # Use the nominal area as is
                input_pixel_area = input_pixflux_area

                # If input image is in flux density units, correct the
                # flux for the user-specified change to the spatial dimension
                if resample_utils.is_flux_density(img.meta.bunit_data):
                    input_pixel_area *= self.pscale_ratio
            else:
                img.meta.wcs.array_shape = img.data.shape
                input_pixel_area = compute_wcs_pixel_area(
                    img.meta.wcs,
                    img.data.shape,
                )
                if input_pixel_area is None:
                    raise ValueError(
                        "Unable to compute input pixel area from WCS of input "
                        f"image {repr(img.meta.filename)}."
                    )
                if self.input_pixscale0 is None:
                    self.input_pixscale0 = np.rad2deg(
                        math.sqrt(input_pixel_area)
                    )
                    if self._recalc_pscale_ratio:
                        self.pscale_ratio = self.pscale / self.input_pixscale0
            iscale = math.sqrt(input_pixflux_area / input_pixel_area)
        else:
            iscale = 1.0
        return iscale

    def resample_group(self, input_models, indices, compute_error=False):
        """Apply resample_many_to_many for one group

        Parameters
        ----------
        input_models : ModelLibrary

        indices : list

        compute_error : bool, optional
            If True, an approximate error image will be resampled
            alongside the science image.
        """
        output_model = None

        # Initialize the output with the wcs
        driz = Drizzle(
            out_shape=self.output_array_shape,
            kernel=self.kernel,
            fillval=self.fillval,
            disable_ctx=True,
        )
        # Also make a temporary model to hold error data
        if compute_error:
            driz_error = Drizzle(
                out_shape=self.output_array_shape,
                kernel=self.kernel,
                fillval=self.fillval,
                disable_ctx=True,
            )
        log.info(f"{len(indices)} exposures to drizzle together")
        for index in indices:
            img = input_models.borrow(index)
            if output_model is None:
                output_model = self._create_output_model(
                    ref_input_model=img
                )
                # Determine output file type from input exposure filenames
                # Use this for defining the output filename
                indx = img.meta.filename.rfind('.')
                output_type = img.meta.filename[indx:]
                output_root = '_'.join(img.meta.filename.replace(
                    output_type,
                    ''
                ).split('_')[:-1])
                output_model.meta.filename = (
                    f'{output_root}_'
                    f'{self.intermediate_suffix}{output_type}'
                )
                copy_asn_info_from_library(input_models, output_model)

            if isinstance(img, datamodels.SlitModel):
                # must call this explicitly to populate area extension
                # although the existence of this extension may not be necessary
                img.area = img.area

            iscale = self._get_intensity_scale(img)
            log.debug(f'Using intensity scale iscale={iscale}')

            inwht = resample_utils.build_driz_weight(
                img,
                weight_type=self.weight_type,
                good_bits=self.good_bits
            )

            # apply sky subtraction
            blevel = img.meta.background.level
            if not img.meta.background.subtracted and blevel is not None:
                data = img.data - blevel
            else:
                data = img.data

            xmin, xmax, ymin, ymax = resample_range(
                data.shape,
                img.meta.wcs.bounding_box
            )
            pixmap = resample_utils.calc_gwcs_pixmap(
                img.meta.wcs,
                self.output_wcs,
                img.data.shape,
            )

            driz.add_image(
                data=data,
                exptime=img.meta.exposure.exposure_time,  # GWCSDrizzle.add_image param default was 1.0
                pixmap=pixmap,
                scale=iscale,
                weight_map=inwht,
                wht_scale=1.0,  # hard-coded for JWST count-rate data
                pixfrac=self.pixfrac,
                in_units="cps",  # GWCSDrizzle.add_image param default
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
            )
            del data

            # make an approximate error image by drizzling it
            # in the same way the image is handled
            if compute_error:
                driz_error.add_image(
                    data=img.err,
                    exptime=img.meta.exposure.exposure_time,  # GWCSDrizzle.add_image param default
                    pixmap=pixmap,
                    scale=iscale,
                    weight_map=inwht,
                    wht_scale=1.0,  # hard-coded for JWST count-rate data
                    pixfrac=self.pixfrac,
                    in_units="cps",  # GWCSDrizzle.add_image param default
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                )
            input_models.shelve(img, index, modify=False)
            del img

        output_model.data = driz.out_img
        output_model.wht = driz.out_wht
        del driz
        # copy the drizzled error into the output model
        if compute_error:
            output_model.err = driz_error.out_img
            del driz_error

        return output_model

    def resample_many_to_many(self, input_models):
        """Resample many inputs to many outputs where outputs have a common frame.

        Coadd only different detectors of the same exposure, i.e. map NRCA5 and
        NRCB5 onto the same output image, as they image different areas of the
        sky.

        Used for outlier detection
        """
        output_models = []
        for group_id, indices in input_models.group_indices.items():

            output_model = self.resample_group(input_models, indices)

            if not self.in_memory:
                # Write out model to disk, then return filename
                output_name = output_model.meta.filename
                if self.output_dir is not None:
                    output_name = os.path.join(self.output_dir, output_name)
                output_model.save(output_name)
                log.info(f"Saved model in {output_name}")
                output_models.append(output_name)
            else:
                output_models.append(output_model)

        if self.in_memory:
            # build ModelLibrary as a list of in-memory models
            return ModelLibrary(output_models, on_disk=False)
        else:
            # build ModelLibrary as an association from the output files
            # this saves memory if there are multiple groups
            asn = asn_from_list(output_models, product_name='outlier_i2d')
            asn_dict = json.loads(asn.dump()[1])  # serializes the asn and converts to dict
            return ModelLibrary(asn_dict, on_disk=True)

    def resample_many_to_one(self, input_models):
        """Resample and coadd many inputs to a single output.

        Used for stage 3 resampling
        """
        output_model = None

        if self.blendheaders:
            blender = ModelBlender(
                blend_ignore_attrs=[
                    'meta.photometry.pixelarea_steradians',
                    'meta.photometry.pixelarea_arcsecsq',
                    'meta.filename',
                ]
            )

        driz = Drizzle(
            out_shape=self.output_array_shape,
            kernel=self.kernel,
            fillval=self.fillval,
            max_ctx_id=len(input_models),
            disable_ctx=False,
        )
        self._init_variance_arrays()
        self._init_exptime_counters()

        log.info("Resampling science and variance data")

        leading_group_idx = [v[0] for v in input_models.group_indices.values()]

        with input_models:
            for idx, img in enumerate(input_models):
                if output_model is None:
                    output_model = self._create_output_model(
                        ref_input_model=img
                    )
                    # Determine output file type from input exposure filenames
                    # Use this for defining the output filename
                    output_model.meta.filename = self.output_filename
                    output_model.meta.resample.weight_type = self.weight_type
                    output_model.meta.resample.pointings = len(input_models.group_names)

                    # copy over asn information
                    copy_asn_info_from_library(input_models, output_model)

                if idx in leading_group_idx:
                    self._update_exptime(img)

                if self.blendheaders:
                    blender.accumulate(img)

                iscale = self._get_intensity_scale(img)
                log.debug(f'Using intensity scale iscale={iscale}')

                inwht = resample_utils.build_driz_weight(
                    img,
                    weight_type=self.weight_type,
                    good_bits=self.good_bits,
                )

                # apply sky subtraction
                blevel = img.meta.background.level
                if not img.meta.background.subtracted and blevel is not None:
                    data = img.data - blevel
                else:
                    data = img.data.copy()

                in_image_limits = resample_range(
                    data.shape,
                    img.meta.wcs.bounding_box
                )
                xmin, xmax, ymin, ymax = in_image_limits

                pixmap = resample_utils.calc_gwcs_pixmap(
                    img.meta.wcs,
                    output_model.meta.wcs,
                    data.shape,
                )

                driz.add_image(
                    data=data,
                    exptime=img.meta.exposure.exposure_time,  # GWCSDrizzle.add_image param default
                    pixmap=pixmap,
                    scale=iscale,
                    weight_map=inwht,
                    wht_scale=1.0,  # hard-coded for JWST count-rate data
                    pixfrac=self.pixfrac,
                    in_units="cps",  # GWCSDrizzle.add_image param default
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                )
                # Resample variance arrays in input_models to output_model
                self._resample_variance_arrays(
                    model=img,
                    iscale=iscale,
                    inwht=inwht,
                    pixmap=pixmap,
                    in_image_limits=in_image_limits,
                    output_shape=self.output_array_shape,
                )

                del data, inwht

                input_models.shelve(img)

        # Since the context array is dynamic, it must be re-assigned
        # back to the product's `con` attribute.
        output_model.data = driz.out_img
        output_model.wht = driz.out_wht
        if driz.out_ctx is not None:
            output_model.con = driz.out_ctx

        del driz

        # compute final variances:
        self._compute_resample_variance_totals(output_model)

        var_components = [
            output_model.var_rnoise,
            output_model.var_poisson,
            output_model.var_flat
        ]
        output_model.err = np.sqrt(np.nansum(var_components, axis=0))

        # nansum returns zero for input that is all NaN -
        # set those values to NaN instead
        all_nan = np.all(np.isnan(var_components), axis=0)
        output_model.err[all_nan] = np.nan

        if self.blendheaders:
            blender.finalize_model(output_model)

        self._get_exptime_totals(output_model)

        return ModelLibrary([output_model,], on_disk=False)

    def _init_variance_arrays(self):
        shape = self.output_array_shape
        self._weighted_rn_var = np.full(shape, np.nan, dtype=np.float32)
        self._weighted_pn_var = np.full(shape, np.nan, dtype=np.float32)
        self._weighted_flat_var = np.full(shape, np.nan, dtype=np.float32)
        self._total_weight_rn_var = np.zeros(shape, dtype=np.float32)
        self._total_weight_pn_var = np.zeros(shape, dtype=np.float32)
        self._total_weight_flat_var = np.zeros(shape, dtype=np.float32)

    def _resample_variance_arrays(self, model, iscale, inwht, pixmap,
                                  in_image_limits, output_shape):
        xmin, xmax, ymin, ymax = in_image_limits

        # Do the read noise variance first, so it can be
        # used for weights if needed
        rn_var = self._resample_one_variance_array(
            "var_rnoise",
            input_model=model,
            iscale=iscale,
            inwht=inwht,
            pixmap=pixmap,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )

        # Find valid weighting values in the variance
        if rn_var is not None:
            mask = (rn_var > 0) & np.isfinite(rn_var)
        else:
            mask = np.full_like(rn_var, False)

        # Set the weight for the image from the weight type
        weight = np.ones(output_shape)
        if self.weight_type == "ivm" and rn_var is not None:
            weight[mask] = rn_var[mask] ** -1
        elif self.weight_type == "exptime":
            if resample_utils.check_for_tmeasure(model):
                weight[:] = model.meta.exposure.measurement_time
            else:
                weight[:] = model.meta.exposure.exposure_time

        # Weight and add the readnoise variance
        # Note: floating point overflow is an issue if variance weights
        # are used - it can't be squared before multiplication
        if rn_var is not None:
            mask = (rn_var >= 0) & np.isfinite(rn_var) & (weight > 0)
            self._weighted_rn_var[mask] = np.nansum(
                [
                    self._weighted_rn_var[mask],
                    rn_var[mask] * weight[mask] * weight[mask]
                ],
                axis=0
            )
            self._total_weight_rn_var[mask] += weight[mask]

        # Now do poisson and flat variance, updating only valid new values
        # (zero is a valid value; negative, inf, or NaN are not)
        pn_var = self._resample_one_variance_array(
            "var_poisson",
            input_model=model,
            iscale=iscale,
            inwht=inwht,
            pixmap=pixmap,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )
        if pn_var is not None:
            mask = (pn_var >= 0) & np.isfinite(pn_var) & (weight > 0)
            self._weighted_pn_var[mask] = np.nansum(
                [
                    self._weighted_pn_var[mask],
                    pn_var[mask] * weight[mask] * weight[mask]
                ],
                axis=0
            )
            self._total_weight_pn_var[mask] += weight[mask]

        flat_var = self._resample_one_variance_array(
            "var_flat",
            input_model=model,
            iscale=iscale,
            inwht=inwht,
            pixmap=pixmap,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )
        if flat_var is not None:
            mask = (flat_var >= 0) & np.isfinite(flat_var) & (weight > 0)
            self._weighted_flat_var[mask] = np.nansum(
                [
                    self._weighted_flat_var[mask],
                    flat_var[mask] * weight[mask] * weight[mask]
                ],
                axis=0
            )
            self._total_weight_flat_var[mask] += weight[mask]

    def _compute_resample_variance_totals(self, output_model):
        # Divide by the total weights, squared, and set in the output model.
        # Zero weight and missing values are NaN in the output.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "invalid value*", RuntimeWarning)
            warnings.filterwarnings("ignore", "divide by zero*", RuntimeWarning)

            output_variance = (
                self._weighted_rn_var / self._total_weight_rn_var /
                self._total_weight_rn_var
            )
            setattr(output_model, "var_rnoise", output_variance)

            output_variance = (
                self._weighted_pn_var / self._total_weight_pn_var /
                self._total_weight_pn_var
            )
            setattr(output_model, "var_poisson", output_variance)

            output_variance = (
                self._weighted_flat_var / self._total_weight_flat_var /
                self._total_weight_flat_var
            )
            setattr(output_model, "var_flat", output_variance)

        del (
            self._weighted_rn_var,
            self._weighted_pn_var,
            self._weighted_flat_var,
            self._total_weight_rn_var,
            self._total_weight_pn_var,
            self._total_weight_flat_var,
        )

    def _resample_one_variance_array(self, name, input_model, iscale,
                                     inwht, pixmap,
                                     xmin=None, xmax=None, ymin=None, ymax=None):
        """Resample one variance image from an input model.

        The error image is passed to drizzle instead of the variance, to
        better match kernel overlap and user weights to the data, in the
        pixel averaging process. The drizzled error image is squared before
        returning.
        """
        variance = getattr(input_model, name)
        if variance is None or variance.size == 0:
            log.debug(
                f"No data for '{name}' for model "
                f"{repr(input_model.meta.filename)}. Skipping ..."
            )
            return

        elif variance.shape != input_model.data.shape:
            log.warning(
                f"Data shape mismatch for '{name}' for model "
                f"{repr(input_model.meta.filename)}. Skipping ..."
            )
            return

        output_shape = self.output_array_shape

        # Resample the error array. Fill "unpopulated" pixels with NaNs.
        driz = Drizzle(
            out_shape=output_shape,
            kernel=self.kernel,
            fillval=np.nan,
            disable_ctx=True
        )

        log.debug(f"Pixmap shape: {pixmap[:,:,0].shape}")
        log.debug(f"Input Sci shape: {variance.shape}")
        log.debug(f"Output Sci shape: {output_shape}")

        # Call 'drizzle' to perform image combination
        log.info(f"Drizzling {variance.shape} --> {output_shape}")

        driz.add_image(
            data=np.sqrt(variance),
            exptime=input_model.meta.exposure.exposure_time,
            pixmap=pixmap,
            scale=iscale,
            weight_map=inwht,
            wht_scale=1.0,  # hard-coded for JWST count-rate data
            pixfrac=self.pixfrac,
            in_units="cps",
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )

        return driz.out_img ** 2

    def _init_exptime_counters(self):
        self._total_exposure_time = 0.
        self._exposure_times = {'start': [], 'end': []}
        self._duration = 0.0
        self._total_measurement_time = 0.0
        self._measurement_time_failures = []

    def _update_exptime(self, model):
        self._total_exposure_time += model.meta.exposure.exposure_time
        if not resample_utils.check_for_tmeasure(model):
            self._measurement_time_failures.append(1)
        else:
            self._total_measurement_time += model.meta.exposure.measurement_time
            self._measurement_time_failures.append(0)
        self._exposure_times['start'].append(model.meta.exposure.start_time)
        self._exposure_times['end'].append(model.meta.exposure.end_time)
        self._duration += model.meta.exposure.duration

    def _get_exptime_totals(self, output_model):
        # Update some basic exposure time values based on output_model
        output_model.meta.exposure.exposure_time = self._total_exposure_time
        if not any(self._measurement_time_failures):
            output_model.meta.exposure.measurement_time = self._total_measurement_time
        output_model.meta.exposure.start_time = min(self._exposure_times['start'])
        output_model.meta.exposure.end_time = max(self._exposure_times['end'])

        # Update other exposure time keywords:
        # XPOSURE (identical to the total effective exposure time, EFFEXPTM)
        xposure = self._total_exposure_time
        output_model.meta.exposure.effective_exposure_time = xposure
        # DURATION (identical to TELAPSE, elapsed time)
        output_model.meta.exposure.duration = self._duration
        output_model.meta.exposure.elapsed_exposure_time = self._duration

        del (
            self._total_exposure_time,
            self._exposure_times,
            self._duration,
            self._total_measurement_time,
            self._measurement_time_failures,
        )

    def update_exposure_times(self, output_model, input_models):
        """Modify exposure time metadata in-place"""
        self._init_exptime_counters()
        with input_models:
            for _, indices in input_models.group_indices.items():
                model = input_models.borrow(indices[0])
                self._update_exptime(model)
        self._get_exptime_totals(output_model)


def _get_boundary_points(xmin, xmax, ymin, ymax, dx=None, dy=None, shrink=0):
    """
    xmin, xmax, ymin, ymax - integer coordinates of pixel boundaries
    step - distance between points along an edge
    shrink - number of pixels by which to reduce `shape`

    Returns a list of points and the area of the rectangle
    """
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1

    if dx is None:
        dx = nx
    if dy is None:
        dy = ny

    if nx - 2 * shrink < 1 or ny - 2 * shrink < 1:
        raise ValueError("Image size is too small.")

    sx = max(1, int(np.ceil(nx / dx)))
    sy = max(1, int(np.ceil(ny / dy)))

    xmin += shrink
    xmax -= shrink
    ymin += shrink
    ymax -= shrink

    size = 2 * sx + 2 * sy
    x = np.empty(size)
    y = np.empty(size)

    b = np.s_[0:sx]  # bottom edge
    r = np.s_[sx:sx + sy]  # right edge
    t = np.s_[sx + sy:2 * sx + sy]  # top edge
    l = np.s_[2 * sx + sy:2 * sx + 2 * sy]  # left

    x[b] = np.linspace(xmin, xmax, sx, False)
    y[b] = ymin
    x[r] = xmax
    y[r] = np.linspace(ymin, ymax, sy, False)
    x[t] = np.linspace(xmax, xmin, sx, False)
    y[t] = ymax
    x[l] = xmin
    y[l] = np.linspace(ymax, ymin, sy, False)

    area = (xmax - xmin) * (ymax - ymin)
    center = (0.5 * (xmin + xmax), 0.5 * (ymin + ymax))

    return x, y, area, center, b, r, t, l


def compute_image_pixel_area(wcs):
    """ Computes pixel area in steradians.
    """
    if wcs.array_shape is None:
        raise ValueError("WCS must have array_shape attribute set.")

    valid_polygon = False
    spatial_idx = np.where(np.array(wcs.output_frame.axes_type) == 'SPATIAL')[0]

    ny, nx = wcs.array_shape
    ((xmin, xmax), (ymin, ymax)) = wcs.bounding_box

    xmin = max(0, int(xmin + 0.5))
    xmax = min(nx - 1, int(xmax - 0.5))
    ymin = max(0, int(ymin + 0.5))
    ymax = min(ny - 1, int(ymax - 0.5))
    if xmin > xmax:
        (xmin, xmax) = (xmax, xmin)
    if ymin > ymax:
        (ymin, ymax) = (ymax, ymin)

    k = 0
    dxy = [1, -1, -1, 1]
    ra, dec, center = np.nan, np.nan, (np.nan, np.nan)
    while xmin < xmax and ymin < ymax:
        try:
            x, y, image_area, center, b, r, t, l = _get_boundary_points(
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
                dx=min((xmax - xmin) // 4, 15),
                dy=min((ymax - ymin) // 4, 15)
            )
        except ValueError:
            return None

        world = wcs(x, y)
        ra = world[spatial_idx[0]]
        dec = world[spatial_idx[1]]

        limits = [ymin, xmax, ymax, xmin]

        for j in range(4):
            sl = [b, r, t, l][k]
            if not (np.all(np.isfinite(ra[sl])) and
                    np.all(np.isfinite(dec[sl]))):
                limits[k] += dxy[k]
                k = (k + 1) % 4
                break
            k = (k + 1) % 4
        else:
            valid_polygon = True
            break

        ymin, xmax, ymax, xmin = limits

    if not valid_polygon:
        return None

    world = wcs(*center)
    wcenter = (world[spatial_idx[0]], world[spatial_idx[1]])

    sky_area = SphericalPolygon.from_radec(ra, dec, center=wcenter).area()
    if sky_area > 2 * np.pi:
        log.warning(
            "Unexpectedly large computed sky area for an image. "
            "Setting area to: 4*Pi - area"
        )
        sky_area = 4 * np.pi - sky_area
    pix_area = sky_area / image_area

    return pix_area


def copy_asn_info_from_library(library, output_model):
    """
    Transfer association information from the input library to the output model.

    Parameters
    ----------
    library : ModelLibrary
        The input library of data models.

    output_model : DataModel
        The output data model to which the association information will be copied.
    """
    if not hasattr(library, "asn"):
        # No ASN table, occurs when input comes from ModelContainer in spectroscopic modes
        # in this case do nothing; the asn info will be passed along later
        # by code inside ResampleSpecStep
        return
    if (asn_pool := library.asn.get("asn_pool", None)) is not None:
        output_model.meta.asn.pool_name = asn_pool
    if (
        asn_table_name := library.asn.get("table_name", None)
    ) is not None:
        output_model.meta.asn.table_name = asn_table_name
