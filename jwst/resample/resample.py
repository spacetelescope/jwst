import logging
import json
from pathlib import Path
import re

import numpy as np

from spherical_geometry.polygon import SphericalPolygon

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels.dqflags import pixel

from stcal.resample import Resample
from stcal.resample.utils import is_imaging_wcs

from jwst.datamodels import ModelLibrary
from jwst.associations.asn_from_list import asn_from_list

from jwst.model_blender.blender import ModelBlender
from jwst.resample import resample_utils
from jwst.assign_wcs import util as assign_wcs_util


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = [
    "input_jwst_model_to_dict",
    "is_imaging_wcs",
    "ResampleImage",
]

_SUPPORTED_CUSTOM_WCS_PARS = [
    "pixel_scale_ratio",
    "pixel_scale",
    "output_shape",
    "crpix",
    "crval",
    "rotation",
]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ResampleImage(Resample):
    """Resample imaging data."""

    dq_flag_name_map = pixel

    def __init__(
        self,
        input_models,
        pixfrac=1.0,
        kernel="square",
        fillval="NAN",
        weight_type="ivm",
        good_bits=0,
        blendheaders=True,
        output_wcs=None,
        wcs_pars=None,
        output=None,
        enable_ctx=True,
        enable_var=True,
        report_var=True,
        compute_err=None,
        asn_id=None,
    ):
        """
        Initialize the ResampleImage object.

        Parameters
        ----------
        input_models : ModelLibrary
            A `ModelLibrary`-based object allowing iterating over
            all contained models of interest.

        pixfrac : float, optional
            The fraction of a pixel that the pixel flux is confined to. The
            default value of 1 has the pixel flux evenly spread across the
            image. A value of 0.5 confines it to half a pixel in the linear
            dimension, so the flux is confined to a quarter of the pixel area
            when the square kernel is used.

        kernel : {"square", "gaussian", "point", "turbo", "lanczos2", "lanczos3"}, optional
            The name of the kernel used to combine the input. The choice of
            kernel controls the distribution of flux over the kernel.
            The square kernel is the default.

            .. warning::
               The "gaussian" and "lanczos2/3" kernels **DO NOT**
               conserve flux.

        fillval : float, None, str, optional
            The value of output pixels that did not have contributions from
            input images' pixels. When ``fillval`` is either `None` or
            ``"INDEF"`` and ``out_img`` is provided, the values of ``out_img``
            will not be modified. When ``fillval`` is either `None` or
            ``"INDEF"`` and ``out_img`` is **not provided**, the values of
            ``out_img`` will be initialized to `numpy.nan`. If ``fillval``
            is a string that can be converted to a number, then the output
            pixels with no contributions from input images will be set to this
            ``fillval`` value.

        weight_type : {"exptime", "ivm"}, optional
            The weighting type for adding models' data. For
            ``weight_type="ivm"`` (the default), the weighting will be
            determined per-pixel using the inverse of the read noise
            (VAR_RNOISE) array stored in each input image.
            If the ``VAR_RNOISE`` array does not exist,
            the variance is set to 1 for all pixels (i.e., equal weighting).
            If ``weight_type="exptime"``, the weight will be set equal
            to the measurement time when available and to
            the exposure time otherwise.

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

        output_wcs : dict, None, optional
            Specifies output WCS as a dictionary
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

        enable_ctx : bool, optional
            Indicates whether to create a context image. If ``disable_ctx``
            is set to `True`, parameters ``out_ctx``, ``begin_ctx_id``, and
            ``max_ctx_id`` will be ignored.

        enable_var : bool, optional
            Indicates whether to resample variance arrays.

        report_var : bool, optional
            Indicates whether to report variance arrays in the output model.
            In order to get an error array when compute_err=from_var, enable_var
            must be True, but sometimes it's useful not to save var_rnoise,
            var_flat, and var_poisson arrays to decrease output file size.

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

        asn_id : str, None, optional
            The association id. The id is what appears in
            the :ref:`asn-jwst-naming`.
        """
        self.input_models = input_models
        self.output_jwst_model = None
        self._report_var = report_var

        self.output_dir = None
        self.output_filename = output
        if output is not None and ".fits" not in str(output):
            self.output_dir = output
            self.output_filename = None
        self.intermediate_suffix = "outlier_i2d"

        self.blendheaders = blendheaders
        if blendheaders:
            self._blender = ModelBlender(
                blend_ignore_attrs=[
                    "meta.photometry.pixelarea_steradians",
                    "meta.photometry.pixelarea_arcsecsq",
                    "meta.filename",
                ]
            )

        self.asn_id = asn_id

        # check wcs_pars has supported keywords:
        if wcs_pars is None:
            wcs_pars = {}
        elif wcs_pars:
            unsup = []
            unsup = set(wcs_pars.keys()).difference(_SUPPORTED_CUSTOM_WCS_PARS)
            if unsup:
                raise KeyError(f"Unsupported custom WCS parameters: {','.join(map(repr, unsup))}.")

        if output_wcs is None:
            # determine output WCS:
            shape = wcs_pars.get("output_shape")
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
            if wcs_pars:
                log.warning("Ignoring 'wcs_pars' since 'output_wcs' is not None.")
            if output_wcs["wcs"].array_shape is None:
                raise ValueError(
                    "Custom WCS objects must have the 'array_shape' attribute set (defined)."
                )

        super().__init__(
            n_input_models=len(input_models),
            pixfrac=pixfrac,
            kernel=kernel,
            fillval=fillval,
            weight_type=weight_type,
            good_bits=good_bits,
            output_wcs=output_wcs,
            enable_ctx=enable_ctx,
            enable_var=enable_var,
            compute_err=compute_err,
        )

    def input_model_to_dict(self, model, weight_type, enable_var, compute_err):
        """
        Convert a data model to a dictionary of keywords and values expected by `stcal.resample`.

        Parameters
        ----------
        model : DataModel
            A JWST data model.
        weight_type : str
            The weighting type for adding models' data.
        enable_var : bool
            Indicates whether to resample variance arrays.
        compute_err : str
            The method to compute the output model's error array.

        Returns
        -------
        dict
            A dictionary of keywords and values expected by `stcal.resample`.
        """
        return input_jwst_model_to_dict(
            model=model, weight_type=weight_type, enable_var=enable_var, compute_err=compute_err
        )

    def create_output_jwst_model(self, ref_input_model=None):
        """
        Create a new blank model and update its meta with info from ``ref_input_model``.

        Parameters
        ----------
        ref_input_model : `~jwst.datamodels.JwstDataModel`, optional
            The reference input model from which to copy meta data.

        Returns
        -------
        ImageModel
            A new blank model with updated meta data.
        """
        output_model = datamodels.ImageModel(None)  # tuple(self.output_wcs.array_shape))

        # update meta data and wcs
        if ref_input_model is not None:
            output_model.update(ref_input_model)
        output_model.meta.wcs = self.output_wcs
        return output_model

    def update_output_model(self, model, info_dict):
        """
        Add meta information to the output model.

        Parameters
        ----------
        model : ImageModel
            The output model to be updated.
        info_dict : dict
            A dictionary containing information about the resampling process.
        """
        model.data = info_dict["data"]
        model.wht = info_dict["wht"]
        if self._enable_ctx:
            model.con = info_dict["con"]
        if self._compute_err:
            model.err = info_dict["err"]
        elif model.meta.hasattr("bunit_err"):
            # bunit_err metadata is mapped to the err extension, so it must be removed
            # in order to fully remove the err extension.
            del model.meta.bunit_err
        if self._enable_var and self._report_var:
            model.var_rnoise = info_dict["var_rnoise"]
            model.var_flat = info_dict["var_flat"]
            model.var_poisson = info_dict["var_poisson"]

        model.meta.wcs = info_dict["wcs"]
        model.meta.photometry.pixelarea_steradians = info_dict["pixelarea_steradians"]
        model.meta.photometry.pixelarea_arcsecsq = info_dict["pixelarea_arcsecsq"]

        model.meta.resample.pointings = info_dict["pointings"]
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
        model.meta.exposure.effective_exposure_time = info_dict["exposure_time"]
        model.meta.exposure.elapsed_exposure_time = info_dict["elapsed_exposure_time"]

    def add_model(self, model):
        """
        Add a single input model to the resampling.

        Parameters
        ----------
        model : ImageModel
            A JWST data model to be resampled.
        """
        super().add_model(
            self.input_model_to_dict(
                model,
                weight_type=self.weight_type,
                enable_var=self._enable_var,
                compute_err=self._compute_err,
            )
        )
        if self.output_jwst_model is None:
            self.output_jwst_model = self.create_output_jwst_model(ref_input_model=model)
        if self.blendheaders:
            self._blender.accumulate(model)

    def finalize(self):
        """Perform final computations and set output model values and metadata."""
        if self.blendheaders:
            self._blender.finalize_model(self.output_jwst_model)
        super().finalize()

        self.update_output_model(
            self.output_jwst_model,
            self.output_model,
        )

        if is_imaging_wcs(self.output_jwst_model.meta.wcs):
            # only for an imaging WCS:
            self.update_fits_wcsinfo(self.output_jwst_model)
            assign_wcs_util.update_s_region_imaging(self.output_jwst_model)

        self.output_jwst_model.meta.cal_step.resample = "COMPLETE"

    def reset_arrays(self, n_input_models=None):
        """
        Initialize/reset between finalize() and add_model() calls.

        Resets or re-initializes `Drizzle` objects, `ModelBlender`, output model
        and arrays, and time counters. Output WCS and shape are not modified
        from `Resample` object initialization. This method needs to be called
        before calling :py:meth:`add_model` for the first time after
        :py:meth:`finalize` was previously called.

        Parameters
        ----------
        n_input_models : int, None, optional
            Number of input models expected to be resampled. When provided,
            this is used to estimate memory requirements and optimize memory
            allocation for the context array.
        """
        super().reset_arrays(n_input_models=n_input_models)
        if self.blendheaders:
            self._blender = ModelBlender(
                blend_ignore_attrs=[
                    "meta.photometry.pixelarea_steradians",
                    "meta.photometry.pixelarea_arcsecsq",
                    "meta.filename",
                ]
            )
        self.output_jwst_model = None

    def resample_group(self, indices):
        """
        Resample multiple input images belonging to a single ``group_id``.

        If ``output_jwst_model``
        was created by a previous call to this method, ``output_jwst_model``
        as well as other arrays (weights, context, etc.) will be cleared.
        Upon completion, this method calls :py:meth:`finalize` to compute
        final values for various attributes of the resampled model
        (e.g., exposure start and end times, etc.)

        Parameters
        ----------
        indices : list
            Indices of models in ``input_models`` model library (used
            to initialize this object) that have the same ``group_id``
            and need to be resampled together.

        Returns
        -------
        output_jwst_model
            Resampled model with populated data, weights, error arrays and
            other attributes.
        """
        if self.output_jwst_model is not None:
            self.reset_arrays(n_input_models=len(indices))

        output_model_filename = ""

        log.info(f"{len(indices)} exposures to drizzle together")
        first = True
        with self.input_models:
            for index in indices:
                model = self.input_models.borrow(index)
                model_modified = False
                if self.output_jwst_model is None:
                    # Determine output file type from input exposure filenames
                    # Use this for defining the output filename
                    indx = model.meta.filename.rfind(".")
                    output_type = model.meta.filename[indx:]
                    output_root = "_".join(
                        model.meta.filename.replace(output_type, "").split("_")[:-1]
                    )
                    output_model_filename = f"{output_root}_{self.intermediate_suffix}{output_type}"

                if isinstance(model, datamodels.SlitModel):
                    # must call this explicitly to populate area extension
                    # although the existence of this extension may not be
                    # necessary
                    model.area = model.area
                    model_modified = True

                self.add_model(model)
                if first:
                    self.output_jwst_model.meta.bunit_data = model.meta.bunit_data
                    first = False
                self.input_models.shelve(model, index, modify=model_modified)
                del model

        self.finalize()
        copy_asn_info_from_library(self.input_models, self.output_jwst_model)
        self.output_jwst_model.meta.filename = output_model_filename
        return self.output_jwst_model

    def resample_many_to_many(self, in_memory=True):
        """
        Resample many inputs to many outputs where outputs have a common frame.

        Coadd only different detectors of the same exposure, i.e. map NRCA5 and
        NRCB5 onto the same output image, as they image different areas of the
        sky.

        Used for outlier detection.

        Parameters
        ----------
        in_memory : bool, optional
            Indicates whether to return a `ModelLibrary` with resampled models
            loaded in memory or whether to serialize resampled models to
            files on disk and return a `ModelLibrary` with only the associacion
            info. See https://stpipe.readthedocs.io/en/latest/model_library.html#on-disk-mode
            for more details.

        Returns
        -------
        ModelLibrary
            A library of resampled models.
        """
        output_models = []

        for _group_id, indices in self.input_models.group_indices.items():
            output_model = self.resample_group(indices)

            if not in_memory:
                # Write out model to disk, then return filename
                output_name = output_model.meta.filename
                if self.output_dir is not None:
                    output_name = str(Path(self.output_dir) / output_name)
                output_model.save(output_name)
                log.info(f"Saved model in {output_name}")
                output_models.append(output_name)
            else:
                output_models.append(output_model)

        if in_memory:
            # build ModelLibrary as a list of in-memory models
            return ModelLibrary(output_models, on_disk=False)
        else:
            # build ModelLibrary as an association from the output files
            # this saves memory if there are multiple groups
            asn = asn_from_list(output_models, product_name="outlier_i2d", asn_id="abcdefg")
            asn_dict = json.loads(asn.dump()[1])  # serializes the asn and converts to dict
            return ModelLibrary(asn_dict, on_disk=True)

    def resample_many_to_one(self):
        """
        Resample and coadd many inputs to a single output.

        Used for stage 3 resampling.

        Returns
        -------
        ImageModel
            The resampled and coadded image.
        """
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

        Parameters
        ----------
        model : ImageModel
            The resampled image
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
        rm_keys = ["v2_ref", "v3_ref", "ra_ref", "dec_ref", "roll_ref", "v3yangle", "vparity"]
        for key in rm_keys:
            if key in model.meta.wcsinfo.instance:
                del model.meta.wcsinfo.instance[key]


def input_jwst_model_to_dict(model, weight_type, enable_var, compute_err):
    """
    Convert a data model to a dictionary of keywords and values expected by `stcal.resample`.

    Parameters
    ----------
    model : DataModel
        A JWST data model.
    weight_type : str
        The weighting type for adding models' data.
    enable_var : bool
        Indicates whether to resample variance arrays.
    compute_err : str
        The method to compute the output model's error array.

    Returns
    -------
    dict
        A dictionary of keywords and values expected by `stcal.resample`.
    """
    model_dict = {
        # arrays:
        "data": model.data,
        "dq": model.dq,
        # meta:
        "filename": model.meta.filename,
        "group_id": model.meta.group_id,
        "wcs": model.meta.wcs,
        "wcsinfo": model.meta.wcsinfo,
        "bunit_data": model.meta.bunit_data,
        "exposure_time": model.meta.exposure.exposure_time,
        "start_time": model.meta.exposure.start_time,
        "end_time": model.meta.exposure.end_time,
        "duration": model.meta.exposure.duration,
        "measurement_time": model.meta.exposure.measurement_time,
        "pixelarea_steradians": model.meta.photometry.pixelarea_steradians,
        "pixelarea_arcsecsq": model.meta.photometry.pixelarea_arcsecsq,
        "level": model.meta.background.level,  # sky level
        "subtracted": model.meta.background.subtracted,
        # spectroscopy-specific:
        "instrument_name": model.meta.instrument.name,
        "exposure_type": model.meta.exposure.type,
    }

    if enable_var:
        model_dict["var_flat"] = model.var_flat
        model_dict["var_rnoise"] = model.var_rnoise
        model_dict["var_poisson"] = model.var_poisson

    elif weight_type is not None and weight_type.startswith("ivm"):
        model_dict["var_rnoise"] = model.var_rnoise

    if compute_err == "driz_err":
        model_dict["err"] = model.err

    return model_dict


def _get_boundary_points(xmin, xmax, ymin, ymax, dx=None, dy=None, shrink=0):
    """
    Compute list of boundary points for a rectangle.

    Parameters
    ----------
    xmin, xmax, ymin, ymax : int
        Coordinates of pixel boundaries.
    dx, dy : int
        Distance between points along an edge in the X and Y directions, respectively.
    shrink : int
        Number of pixels by which to reduce `shape`

    Returns
    -------
    x, y : numpy.ndarray
        Arrays of X and Y coordinates of the boundary points.
    area : float
        Area of the rectangle.
    center : tuple
        Center of the rectangle.
    b, r, t, l : slice
        Slices for the bottom, right, top, and left edges, respectively.
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
    r = np.s_[sx : sx + sy]  # right edge
    t = np.s_[sx + sy : 2 * sx + sy]  # top edge
    l = np.s_[2 * sx + sy : 2 * sx + 2 * sy]  # left

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
    """
    Compute pixel area in steradians from a WCS.

    Parameters
    ----------
    wcs : gwcs.WCS
        A WCS object.

    Returns
    -------
    float
        Pixel area in steradians.
    """
    if wcs.array_shape is None:
        raise ValueError("WCS must have array_shape attribute set.")

    valid_polygon = False
    spatial_idx = np.where(np.array(wcs.output_frame.axes_type) == "SPATIAL")[0]

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
                dy=min((ymax - ymin) // 4, 15),
            )
        except ValueError:
            return None

        world = wcs(x, y)
        ra = world[spatial_idx[0]]
        dec = world[spatial_idx[1]]

        limits = [ymin, xmax, ymax, xmin]

        for _ in range(4):
            sl = [b, r, t, l][k]
            if not (np.all(np.isfinite(ra[sl])) and np.all(np.isfinite(dec[sl]))):
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
            "Unexpectedly large computed sky area for an image. Setting area to: 4*Pi - area"
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
    if (asn_table_name := library.asn.get("table_name", None)) is not None:
        output_model.meta.asn.table_name = asn_table_name
