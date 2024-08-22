import logging
import warnings
import json

import numpy as np
from astropy.io import fits

from stdatamodels.jwst import datamodels
from stcal.resample import ResampleModelIO, ResampleCoAdd, ResampleSingle
from stcal.resample.utils import get_tmeasure
from drizzle.resample import Drizzle
from stdatamodels.jwst.datamodels.dqflags import pixel
from stdatamodels.properties import ObjectNode

from jwst.datamodels import ModelLibrary
from jwst.associations.asn_from_list import asn_from_list

from . import gwcs_drizzle
from jwst.model_blender.blender import ModelBlender
from jwst.resample import resample_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = [
    "OutputTooLargeError",
    "ResampleJWSTModelIO",
    "ResampleJWSTSingle",
    "ResampleJWSTCoAdd",
    "ResampleData",
]


class OutputTooLargeError(RuntimeError):
    """Raised when the output is too large for in-memory instantiation"""


class ResampleJWSTModelIO(ResampleModelIO):
    attributes_to_meta = {
        "filename": "meta.filename",
        "group_id": "meta.group_id",
        "s_region": "meta.wcsinfo.s_region",
        "wcsinfo": "meta.wcsinfo",
        "wcs": "meta.wcs",
        "exptime": "meta.exptime",
        "exposure_time": "meta.exposure.exposure_time",
        "start_time": "meta.exposure.start_time",
        "end_time": "meta.exposure.end_time",
        "duration": "meta.exposure.duration",
        "measurement_time": "meta.exposure.measurement_time",
        "effective_exposure_time": "meta.exposure.effective_exposure_time",
        "elapsed_exposure_time": "meta.exposure.elapsed_exposure_time",

        "pixelarea_steradians": "meta.photometry.pixelarea_steradians",
        "pixelarea_arcsecsq": "meta.photometry.pixelarea_arcsecsq",

        "level": "meta.background.level",
        "subtracted": "meta.background.subtracted",

        "weight_type": "meta.resample.weight_type",
        "pointings": "meta.resample.pointings",
        "ncoadds": "meta.resample.ncoadds",
    }

    def get_model_attr_value(self, model, attribute_name):
        m = model
        meta_name = ResampleJWSTModelIO.attributes_to_meta.get(
            attribute_name,
            attribute_name
        )
        fields = meta_name.strip().split(".")
        while fields:
            m = getattr(m, fields.pop(0))
        if isinstance(m, ObjectNode):
            return m.instance
        return m

    def set_model_attr_value(self, model, attribute_name, value):
        m = model
        meta_name = ResampleJWSTModelIO.attributes_to_meta.get(
            attribute_name,
            attribute_name
        )
        fields = meta_name.strip().split(".")
        while len(fields) > 1:
            m = getattr(m, fields.pop(0))
        setattr(m, fields.pop(), value)

    def open_model(cls, file_name):
        return datamodels.open(file_name)

    def get_model_array(self, model, array_name, **kwargs):
        if isinstance(model, str):
            model = self.open_model(file_name=model)
        if "default" in kwargs:
            return getattr(model, array_name, kwargs["default"])
        else:
            return getattr(model, array_name)

    def set_model_array(self, model, array_name, data):
        """ model must be an open model - not a file name """
        setattr(model, array_name, data)

    def get_model_meta(self, model, attributes):
        meta = {}
        if isinstance(model, str):
            if 's_region' in attributes:
                attributes.pop(attributes.index('s_region'))
                with fits.open(model) as h:
                    meta['s_region'] = h[('sci', 1)].header['s_region']
            if attributes:
                model = self.open_model(model)

        for f in attributes:
            meta[f] = self.get_model_attr_value(model, attribute_name=f)

        return meta

    def set_model_meta(self, model, attributes):
        """ model must be an open model - not a file name """
        for k, v in attributes.items():
            self.set_model_attr_value(model, attribute_name=k, value=v)

    def close_model(self, model):
        self.save_model(model)
        # model.close()

    def save_model(self, model):
        if model.meta.filename:
            model.write(model.meta.filename, overwrite=True)

    def write_model(self, model, file_name, **kwargs):
        overwrite = kwargs.get("overwrite", False)
        model.write(file_name, overwrite=overwrite)

    def new_model(self, image_shape=None, file_name=None, copy_meta_from=None):
        """ Return a new model for the resampled output """
        model = datamodels.ImageModel(image_shape)
        model.meta.filename = file_name
        if copy_meta_from is not None:
            model.update(copy_meta_from)
        return model


class ResampleJWSTCoAdd(ResampleJWSTModelIO, ResampleCoAdd):
    # resample_array_names = [
    #     {'attr': 'data', 'variance', 'exptime']
    dq_flag_name_map = pixel

    def __init__(self, *args, blendheaders=True, **kwargs):
        super().__init__(*args, **kwargs)
        self._blendheaders = blendheaders

    # FIXME: this method will be moved completely to stcal once we have a method
    #        that can create output wcs from s_region.
    def _compute_output_wcs(self, **wcs_pars):
        """
        returns a distortion-free WCS object and its pixel scale.
        this code should be moved to stcal

        """
        # Define output WCS based on all inputs, including a reference WCS:
        output_shape = wcs_pars.get("output_shape", None)
        output_wcs = resample_utils.make_output_wcs(
            self._input_models,
            ref_wcs=None,
            pscale_ratio=wcs_pars.get("pixel_scale_ratio", 1.0),
            pscale=wcs_pars.get("pixel_scale", None),
            rotation=wcs_pars.get("rotation", 0.0),
            shape=None if output_shape is None else output_shape[::-1],
            crpix=wcs_pars.get("crpix", None),
            crval=wcs_pars.get("crval", None),
        )

        # Estimate output pixel area in Sr. NOTE: in principle we could
        # use the same algorithm as for when output_wcs is provided by the
        # user.
        tr = output_wcs.pipeline[0].transform
        output_pix_area = (
            np.deg2rad(tr['cdelt1'].factor.value) *
            np.deg2rad(tr['cdelt2'].factor.value)
        )
        return output_wcs, output_pix_area

    # TODO: Not sure about funct. signature and also I don't like it needs
    # to open input files again. Should we store meta of all inputs?
    # Should blendmeta.blendmodels be redesigned to blend one meta at a time?
    def blend_output_metadata(self):
        """ Create new output metadata based on blending all input metadata. """

        if not self._blendheaders:
            return

        ignore_list = [
            'meta.photometry.pixelarea_steradians',
            'meta.photometry.pixelarea_arcsecsq',
        ]
        return

    # FIXME: blendmodels must be redesigned to work with model library but
    #        most importantly, see if it can be done one at a time when the
    #        'run()' method is run in order to avoid unnecessary opening/closing
    #        of data models.

        log.info(f'Blending metadata for {self._output_filename}')
        blendmeta.blendmodels(
            self._output_model,
            inputs=self._input_models,
            output=self._output_filename,
            ignore=ignore_list
        )

    def final_post_processing(self):
        # update meta for the output model:
        self._output_model.meta.cal_step.resample = 'COMPLETE'
        _update_fits_wcsinfo(self._output_model)
        util.update_s_region_imaging(self._output_model)
        self._output_model.meta.asn.pool_name = self._input_models.asn.get("pool_name", None)
        self._output_model.meta.asn.table_name = self._input_models.asn.get("table_name", None)
        self._output_model.meta.resample.pixel_scale_ratio = self._pixel_scale_ratio
        self._output_model.meta.resample.pixfrac = self.pixfrac
        self.blend_output_metadata()

    def run(self):
        output_model = super().run()
        ml = ModelLibrary([output_model])
        return ml


class ResampleJWSTSingle(ResampleJWSTModelIO, ResampleSingle):
    dq_flag_name_map = pixel

    def run(self):
        output_models = super().run()
        ml = ModelLibrary(output_models)
        return ml


def _update_fits_wcsinfo(model):
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
    rm_keys = [
        'v2_ref',
        'v3_ref',
        'ra_ref',
        'dec_ref',
        'roll_ref',
        'v3yangle',
        'vparity',
    ]
    for key in rm_keys:
        if key in model.meta.wcsinfo.instance:
            del model.meta.wcsinfo.instance[key]


####################################################
#  Code below was left for spectral data for now   #
####################################################

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
                output_pix_area = compute_image_pixel_area(self.output_wcs)
                if output_pix_area is None:
                    raise ValueError(
                        "Unable to compute output pixel area from 'output_wcs'."
                    )
            else:
                output_pix_area = output_wcs.pixel_area

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
            output_pix_area = (
                np.deg2rad(tr['cdelt1'].factor.value) *
                np.deg2rad(tr['cdelt2'].factor.value)
            )

        if pscale is None:
            pscale = np.rad2deg(np.sqrt(output_pix_area))
            log.info(f'Computed output pixel scale: {3600.0 * pscale} arcsec.')

        self.pscale = pscale  # in deg

        log.debug('Output mosaic size: {}'.format(self.output_wcs.array_shape))

        allowed_memory = kwargs['allowed_memory']
        if allowed_memory is None:
            allowed_memory = os.environ.get('DMODEL_ALLOWED_MEMORY', allowed_memory)
        if allowed_memory:
            allowed_memory = float(allowed_memory)
            # make a small image model to get the dtype
            dtype = datamodels.ImageModel((1, 1)).data.dtype

            # get the available memory
            available_memory = psutil.virtual_memory().available + psutil.swap_memory().total

            # compute the output array size
            required_memory = np.prod(self.output_wcs.array_shape) * dtype.itemsize

            # compare used to available
            used_fraction = required_memory / available_memory
            can_allocate = used_fraction <= allowed_memory
        else:
            can_allocate = True

        if not can_allocate:
            raise OutputTooLargeError(
                f'Combined ImageModel size {self.output_wcs.array_shape} '
                f'requires {bytes2human(required_memory)}. '
                f'Model cannot be instantiated.'
            )
        self.blank_output = datamodels.ImageModel(tuple(self.output_wcs.array_shape))

        # update meta data and wcs
        with input_models:
            example_model = input_models.borrow(0)
            self.blank_output.update(example_model)
            input_models.shelve(example_model, 0, modify=False)
            del example_model
        self.blank_output.meta.wcs = self.output_wcs
        self.blank_output.meta.photometry.pixelarea_steradians = output_pix_area
        self.blank_output.meta.photometry.pixelarea_arcsecsq = (
            output_pix_area * np.rad2deg(3600)**2
        )


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
                input_pixel_area = compute_image_pixel_area(img.meta.wcs)
                if input_pixel_area is None:
                    raise ValueError(
                        "Unable to compute input pixel area from WCS of input "
                        f"image {repr(img.meta.filename)}."
                    )
                if self.input_pixscale0 is None:
                    self.input_pixscale0 = np.rad2deg(
                        np.sqrt(input_pixel_area)
                    )
                    if self._recalc_pscale_ratio:
                        self.pscale_ratio = self.pscale / self.input_pixscale0
            iscale = np.sqrt(input_pixflux_area / input_pixel_area)
        else:
            iscale = 1.0
        return iscale

    def resample_group(self, input_models, indices):
        """Apply resample_many_to_many for one group

        Parameters
        ----------
        input_models : ModelLibrary

        indices : list
        """
        output_model = self.blank_output.copy()

        copy_asn_info_from_library(input_models, output_model)

        with input_models:
            example_image = input_models.borrow(indices[0])

            # Determine output file type from input exposure filenames
            # Use this for defining the output filename
            indx = example_image.meta.filename.rfind('.')
            output_type = example_image.meta.filename[indx:]
            output_root = '_'.join(example_image.meta.filename.replace(
                output_type, '').split('_')[:-1])
            output_model.meta.filename = (
                f'{output_root}_'
                f'{self.intermediate_suffix}{output_type}')
            input_models.shelve(example_image, indices[0], modify=False)
            del example_image

            # Initialize the output with the wcs
            driz = gwcs_drizzle.GWCSDrizzle(output_model, pixfrac=self.pixfrac,
                                            kernel=self.kernel, fillval=self.fillval)

            log.info(f"{len(indices)} exposures to drizzle together")
            for index in indices:
                img = input_models.borrow(index)
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

                xmin, xmax, ymin, ymax = resample_utils._resample_range(
                    data.shape,
                    img.meta.wcs.bounding_box
                )

                driz.add_image(
                    data,
                    img.meta.wcs,
                    iscale=iscale,
                    inwht=inwht,
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax
                )
                del data
                input_models.shelve(img, index, modify=False)
                del img

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

        if not self.in_memory:
            # build ModelLibrary as an association from the output files
            # this saves memory if there are multiple groups
            asn = asn_from_list(output_models, product_name='outlier_i2d')
            asn_dict = json.loads(asn.dump()[1]) # serializes the asn and converts to dict
            return ModelLibrary(asn_dict, on_disk=True)
        # otherwise just build it as a list of in-memory models
        return ModelLibrary(output_models, on_disk=False)

    def resample_many_to_one(self, input_models):
        """Resample and coadd many inputs to a single output.

        Used for stage 3 resampling
        """
        output_model = self.blank_output.copy()
        output_model.meta.filename = self.output_filename
        output_model.meta.resample.weight_type = self.weight_type
        output_model.meta.resample.pointings = len(input_models.group_names)

        # copy over asn information
        copy_asn_info_from_library(input_models, output_model)

        if self.blendheaders:
            blender = ModelBlender(
                blend_ignore_attrs=[
                    'meta.photometry.pixelarea_steradians',
                    'meta.photometry.pixelarea_arcsecsq',
                    'meta.filename',
                ]
            )

        # Initialize the output with the wcs
        driz = gwcs_drizzle.GWCSDrizzle(output_model, pixfrac=self.pixfrac,
                                        kernel=self.kernel, fillval=self.fillval)

        log.info("Resampling science data")
        with input_models:
            for img in input_models:
                if self.blendheaders:
                    blender.accumulate(img)
                iscale = self._get_intensity_scale(img)
                log.debug(f'Using intensity scale iscale={iscale}')
                img.meta.iscale = iscale

                inwht = resample_utils.build_driz_weight(img,
                                                        weight_type=self.weight_type,
                                                        good_bits=self.good_bits)
                # apply sky subtraction
                blevel = img.meta.background.level
                if not img.meta.background.subtracted and blevel is not None:
                    data = img.data - blevel
                else:
                    data = img.data.copy()

                xmin, xmax, ymin, ymax = resample_utils._resample_range(
                    data.shape,
                    img.meta.wcs.bounding_box
                )

                driz.add_image(
                    data,
                    img.meta.wcs,
                    iscale=iscale,
                    inwht=inwht,
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax
                )
                del data, inwht
                input_models.shelve(img)

        if self.blendheaders:
            blender.finalize_model(output_model)

        # Resample variance arrays in input_models to output_model
        self.resample_variance_arrays(output_model, input_models)
        var_components = [
            output_model.var_rnoise,
            output_model.var_poisson,
            output_model.var_flat
        ]
        output_model.err = np.sqrt(np.nansum(var_components,axis=0))

        # nansum returns zero for input that is all NaN -
        # set those values to NaN instead
        all_nan = np.all(np.isnan(var_components), axis=0)
        output_model.err[all_nan] = np.nan

        self.update_exposure_times(output_model, input_models)

        return ModelLibrary([output_model,], on_disk=False)


    def resample_variance_arrays(self, output_model, input_models):
        """Resample variance arrays from input_models to the output_model.

        Variance images from each input model are resampled individually and
        added to a weighted sum. If weight_type is 'ivm', the inverse of the
        resampled read noise variance is used as the weight for all the variance
        components. If weight_type is 'exptime', the exposure time is used.

        The output_model is modified in place.
        """
        log.info("Resampling variance components")
        weighted_rn_var = np.full_like(output_model.data, np.nan)
        weighted_pn_var = np.full_like(output_model.data, np.nan)
        weighted_flat_var = np.full_like(output_model.data, np.nan)
        total_weight_rn_var = np.zeros_like(output_model.data)
        total_weight_pn_var = np.zeros_like(output_model.data)
        total_weight_flat_var = np.zeros_like(output_model.data)
        with input_models:
            for i, model in enumerate(input_models):
                # Do the read noise variance first, so it can be
                # used for weights if needed
                rn_var = self._resample_one_variance_array(
                    "var_rnoise", model, output_model)

                # Find valid weighting values in the variance
                if rn_var is not None:
                    mask = (rn_var > 0) & np.isfinite(rn_var)
                else:
                    mask = np.full_like(rn_var, False)

                # Set the weight for the image from the weight type
                weight = np.ones(output_model.data.shape)
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
                    weighted_rn_var[mask] = np.nansum(
                        [weighted_rn_var[mask],
                        rn_var[mask] * weight[mask] * weight[mask]],
                        axis=0
                    )
                    total_weight_rn_var[mask] += weight[mask]

                # Now do poisson and flat variance, updating only valid new values
                # (zero is a valid value; negative, inf, or NaN are not)
                pn_var = self._resample_one_variance_array(
                    "var_poisson", model, output_model)
                if pn_var is not None:
                    mask = (pn_var >= 0) & np.isfinite(pn_var) & (weight > 0)
                    weighted_pn_var[mask] = np.nansum(
                        [weighted_pn_var[mask],
                        pn_var[mask] * weight[mask] * weight[mask]],
                        axis=0
                    )
                    total_weight_pn_var[mask] += weight[mask]

                flat_var = self._resample_one_variance_array(
                    "var_flat", model, output_model)
                if flat_var is not None:
                    mask = (flat_var >= 0) & np.isfinite(flat_var) & (weight > 0)
                    weighted_flat_var[mask] = np.nansum(
                        [weighted_flat_var[mask],
                        flat_var[mask] * weight[mask] * weight[mask]],
                        axis=0
                    )
                    total_weight_flat_var[mask] += weight[mask]

                del model.meta.iscale
                del weight
                input_models.shelve(model, i)

            # We now have a sum of the weighted resampled variances.
            # Divide by the total weights, squared, and set in the output model.
            # Zero weight and missing values are NaN in the output.
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "invalid value*", RuntimeWarning)
                warnings.filterwarnings("ignore", "divide by zero*", RuntimeWarning)

                output_variance = (weighted_rn_var
                                / total_weight_rn_var / total_weight_rn_var)
                setattr(output_model, "var_rnoise", output_variance)

                output_variance = (weighted_pn_var
                                / total_weight_pn_var / total_weight_pn_var)
                setattr(output_model, "var_poisson", output_variance)

                output_variance = (weighted_flat_var
                                / total_weight_flat_var / total_weight_flat_var)
                setattr(output_model, "var_flat", output_variance)

            del weighted_rn_var, weighted_pn_var, weighted_flat_var
            del total_weight_rn_var, total_weight_pn_var, total_weight_flat_var


    def _resample_one_variance_array(self, name, input_model, output_model):
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

        # Make input weight map
        inwht = resample_utils.build_driz_weight(
            input_model,
            weight_type=self.weight_type,  # weights match science
            good_bits=self.good_bits
        )

        resampled_error = np.zeros_like(output_model.data)
        outwht = np.zeros_like(output_model.data)
        outcon = np.zeros_like(output_model.con)

        xmin, xmax, ymin, ymax = resample_utils._resample_range(
            variance.shape,
            input_model.meta.wcs.bounding_box
        )

        iscale = input_model.meta.iscale

        # Resample the error array. Fill "unpopulated" pixels with NaNs.
        self.drizzle_arrays(
            np.sqrt(variance),
            inwht,
            input_model.meta.wcs,
            output_model.meta.wcs,
            resampled_error,
            outwht,
            outcon,
            iscale=iscale,
            pixfrac=self.pixfrac,
            kernel=self.kernel,
            fillval=np.nan,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax
        )
        return resampled_error ** 2

    def update_exposure_times(self, output_model, input_models):
        """Modify exposure time metadata in-place"""
        total_exposure_time = 0.
        exposure_times = {'start': [], 'end': []}
        duration = 0.0
        total_measurement_time = 0.0
        measurement_time_failures = []
        with input_models:
            for _, indices in input_models.group_indices.items():
                model = input_models.borrow(indices[0])
                total_exposure_time += model.meta.exposure.exposure_time
                if not resample_utils.check_for_tmeasure(model):
                    measurement_time_failures.append(1)
                else:
                    total_measurement_time += model.meta.exposure.measurement_time
                    measurement_time_failures.append(0)
                exposure_times['start'].append(model.meta.exposure.start_time)
                exposure_times['end'].append(model.meta.exposure.end_time)
                duration += model.meta.exposure.duration
                input_models.shelve(model, indices[0], modify=False)

        # Update some basic exposure time values based on output_model
        output_model.meta.exposure.exposure_time = total_exposure_time
        if not any(measurement_time_failures):
            output_model.meta.exposure.measurement_time = total_measurement_time
        output_model.meta.exposure.start_time = min(exposure_times['start'])
        output_model.meta.exposure.end_time = max(exposure_times['end'])

        # Update other exposure time keywords:
        # XPOSURE (identical to the total effective exposure time, EFFEXPTM)
        xposure = total_exposure_time
        output_model.meta.exposure.effective_exposure_time = xposure
        # DURATION (identical to TELAPSE, elapsed time)
        output_model.meta.exposure.duration = duration
        output_model.meta.exposure.elapsed_exposure_time = duration

    @staticmethod
    def drizzle_arrays(insci, inwht, input_wcs, output_wcs, outsci, outwht,
                       outcon, uniqid=1, xmin=0, xmax=0, ymin=0, ymax=0,
                       iscale=1.0, pixfrac=1.0, kernel='square',
                       fillval="NAN", wtscale=1.0):
        """
        Low level routine for performing 'drizzle' operation on one image.

        The interface is compatible with STScI code. All images are Python
        ndarrays, instead of filenames. File handling (input and output) is
        performed by the calling routine.

        Parameters
        ----------

        insci : 2d array
            A 2d numpy array containing the input image to be drizzled.

        inwht : 2d array
            A 2d numpy array containing the pixel by pixel weighting.
            Must have the same dimensions as insci. If none is supplied,
            the weighting is set to one.

        input_wcs : gwcs.WCS object
            The world coordinate system of the input image.

        output_wcs : gwcs.WCS object
            The world coordinate system of the output image.

        outsci : 2d array
            A 2d numpy array containing the output image produced by
            drizzling. On the first call it should be set to zero.
            Subsequent calls it will hold the intermediate results.  This
            is modified in-place.

        outwht : 2d array
            A 2d numpy array containing the output counts. On the first
            call it should be set to zero. On subsequent calls it will
            hold the intermediate results.  This is modified in-place.

        outcon : 2d or 3d array, optional
            A 2d or 3d numpy array holding a bitmap of which image was an input
            for each output pixel. Should be integer zero on first call.
            Subsequent calls hold intermediate results.  This is modified
            in-place.

        uniqid : int, optional
            The id number of the input image. Should be one the first time
            this function is called and incremented by one on each subsequent
            call.

        xmin : int, optional
            This and the following three parameters set a bounding rectangle
            on the input image. Only pixels on the input image inside this
            rectangle will have their flux added to the output image. Xmin
            sets the minimum value of the x dimension. The x dimension is the
            dimension that varies quickest on the image. All four parameters
            are zero based, counting starts at zero.

        xmax : int, optional
            Sets the maximum value of the x dimension on the bounding box
            of the input image. If ``xmax = 0``, no maximum will
            be set in the x dimension (all pixels in a row of the input image
            will be resampled).

        ymin : int, optional
            Sets the minimum value in the y dimension on the bounding box. The
            y dimension varies less rapidly than the x and represents the line
            index on the input image.

        ymax : int, optional
            Sets the maximum value in the y dimension. If ``ymax = 0``,
            no maximum will be set in the y dimension (all pixels in a column
            of the input image will be resampled).

        iscale : float, optional
            A scale factor to be applied to pixel intensities of the
            input image before resampling.

        pixfrac : float, optional
            The fraction of a pixel that the pixel flux is confined to. The
            default value of 1 has the pixel flux evenly spread across the image.
            A value of 0.5 confines it to half a pixel in the linear dimension,
            so the flux is confined to a quarter of the pixel area when the square
            kernel is used.

        kernel: str, optional
            The name of the kernel used to combine the input. The choice of
            kernel controls the distribution of flux over the kernel. The kernel
            names are: "square", "gaussian", "point", "turbo", "lanczos2",
            and "lanczos3". The square kernel is the default.

        fillval: str, optional
            The value a pixel is set to in the output if the input image does
            not overlap it. The default value of NAN sets NaN values.

        Returns
        -------
        A tuple with three values: a version string, the number of pixels
        on the input image that do not overlap the output image, and the
        number of complete lines on the input image that do not overlap the
        output input image.

        """

        # Insure that the fillval parameter gets properly interpreted for use with tdriz
        if util.is_blank(str(fillval)):
            fillval = 'NAN'
        else:
            fillval = str(fillval)

        if insci.dtype > np.float32:
            insci = insci.astype(np.float32)

        # Add input weight image if it was not passed in
        if inwht is None:
            inwht = np.ones_like(insci)

        # Compute what plane of the context image this input would
        # correspond to:
        planeid = int((uniqid - 1) / 32)

        # Check if the context image has this many planes
        if outcon.ndim == 3:
            nplanes = outcon.shape[0]
        elif outcon.ndim == 2:
            nplanes = 1
        else:
            nplanes = 0

        if nplanes <= planeid:
            raise IndexError("Not enough planes in drizzle context image")

        # Alias context image to the requested plane if 3d
        if outcon.ndim == 3:
            outcon = outcon[planeid]

        # Compute the mapping between the input and output pixel coordinates
        # for use in drizzle.cdrizzle.tdriz
        pixmap = resample_utils.calc_gwcs_pixmap(input_wcs, output_wcs, insci.shape)

        log.debug(f"Pixmap shape: {pixmap[:,:,0].shape}")
        log.debug(f"Input Sci shape: {insci.shape}")
        log.debug(f"Output Sci shape: {outsci.shape}")

        log.info(f"Drizzling {insci.shape} --> {outsci.shape}")

        _vers, _nmiss, _nskip = cdrizzle.tdriz(
            insci, inwht, pixmap,
            outsci, outwht, outcon,
            uniqid=uniqid,
            xmin=xmin, xmax=xmax,
            ymin=ymin, ymax=ymax,
            scale=iscale,
            pixfrac=pixfrac,
            kernel=kernel,
            in_units="cps",
            expscale=1.0,
            wtscale=wtscale,
            fillstr=fillval
        )


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
