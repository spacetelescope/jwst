"""Primary code for performing outlier detection on JWST observations."""

from functools import partial
import numpy as np

from astropy.stats import sigma_clip
from scipy import ndimage
from drizzle.cdrizzle import tblot

from .. import datamodels
from ..resample import resample
from ..resample.resample_utils import build_driz_weight, calc_gwcs_pixmap
from ..stpipe.step import Step

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

CRBIT = np.uint32(
    datamodels.dqflags.pixel['DO_NOT_USE'] + datamodels.dqflags.pixel['OUTLIER']
)


__all__ = ["OutlierDetection", "flag_cr", "abs_deriv"]


class OutlierDetection:
    """Main class for performing outlier detection.

    This is the controlling routine for the outlier detection process.
    It loads and sets the various input data and parameters needed by
    the various functions and then controls the operation of this process
    through all the steps used for the detection.

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input model and merges
         them with any user-provided values
      2. Resamples all input images into grouped observation mosaics.
      3. Creates a median image from all grouped observation mosaics.
      4. Blot median image to match each original input image.
      5. Perform statistical comparison between blotted image and original
         image to identify outliers.
      6. Updates input data model DQ arrays with mask of detected outliers.

    """

    default_suffix = 'i2d'

    def __init__(self, input_models, reffiles=None, **pars):
        """
        Initialize the class with input ModelContainers.

        Parameters
        ----------
        input_models : list of DataModels, str
            list of data models as ModelContainer or ASN file,
            one data model for each input image

        pars : dict, optional
            Optional user-specified parameters to modify how outlier_detection
            will operate.  Valid parameters include:
            - resample_suffix

        """
        self.inputs = input_models
        self.reffiles = reffiles

        self.outlierpars = {}
        if 'outlierpars' in reffiles:
            self._get_outlier_pars()
        self.outlierpars.update(pars)
        # Insure that self.input_models always refers to a ModelContainer
        # representation of the inputs

        # Define how file names are created
        self.make_output_path = pars.get(
            'make_output_path',
            partial(Step._make_output_path, None)
        )

    def _convert_inputs(self):
        """Convert input into datamodel required for processing.

        This method converts `self.inputs` into a version of
        `self.input_models` suitable for processing by the class.

        This base class works on imaging data, and relies on use of the
        ModelContainer class as the format needed for processing. However,
        the input may not always be a ModelContainer object, so this method
        will convert the input to a ModelContainer object for processing.
        Additionally, sub-classes may redefine this to set up the input as
        whatever format the sub-class needs for processing.

        """
        bits = self.outlierpars['good_bits']
        if isinstance(self.inputs, datamodels.ModelContainer):
            self.input_models = self.inputs
            self.converted = False
        else:
            self.input_models = datamodels.ModelContainer()
            num_inputs = self.inputs.data.shape[0]
            log.debug("Converting CubeModel to ModelContainer with {} images".
                      format(num_inputs))
            for i in range(self.inputs.data.shape[0]):
                image = datamodels.ImageModel(data=self.inputs.data[i],
                                              err=self.inputs.err[i],
                                              dq=self.inputs.dq[i])
                image.meta = self.inputs.meta
                image.wht = build_driz_weight(image,
                                              weight_type='exptime',
                                              good_bits=bits)
                self.input_models.append(image)
            self.converted = True

    def _get_outlier_pars(self):
        """Extract outlier detection parameters from reference file."""
        # start by interpreting input data models to define selection criteria
        input_dm = self.input_models[0]
        filtname = input_dm.meta.instrument.filter
        if hasattr(self.input_models, 'group_names'):
            num_groups = len(self.input_models.group_names)
        else:
            num_groups = 1

        ref_model = datamodels.OutlierParsModel(self.reffiles['outlierpars'])

        # look for row that applies to this set of input data models
        # NOTE:
        #  This logic could be replaced by a method added to the DrizParsModel
        #  object to select the correct row based on a set of selection
        #  parameters
        row = None
        outlierpars = ref_model.outlierpars_table

        # flag to support wild-card rows in outlierpars table
        filter_match = False
        for n, filt, num in zip(range(1, outlierpars.numimages.shape[0] + 1),
                                outlierpars.filter, outlierpars.numimages):
            # only remember this row if no exact match has already been made
            # for the filter. This allows the wild-card row to be anywhere in
            # the table; since it may be placed at beginning or end of table.

            if filt == "ANY" and not filter_match and num_groups >= num:
                row = n
            # always go for an exact match if present, though...
            if filtname == filt and num_groups >= num:
                row = n
                filter_match = True

        # With presence of wild-card rows, code should never trigger this logic
        if row is None:
            log.error("No row found in %s that matches input data.",
                      self.reffiles)
            raise ValueError

        # read in values from that row for each parameter
        for kw in list(self.outlierpars.keys()):
            self.outlierpars[kw] = \
                                ref_model['outlierpars_table.{0}'.format(kw)]

    def build_suffix(self, **pars):
        """Build suffix.

        Class-specific method for defining the resample_suffix attribute
        using a suffix specific to the sub-class.

        """
        # Parse any user-provided filename suffix for resampled products
        self.resample_suffix = '_outlier_{}.fits'.format(
            pars.get('resample_suffix', self.default_suffix))
        if 'resample_suffix' in pars:
            del pars['resample_suffix']
        log.debug("Defined output product suffix as: {}".format(
            self.resample_suffix))

    def do_detection(self):
        """Flag outlier pixels in DQ of input images."""
        self._convert_inputs()
        self.build_suffix(**self.outlierpars)

        pars = self.outlierpars
        save_intermediate_results = pars['save_intermediate_results']
        if pars['resample_data']:
            # Start by creating resampled/mosaic images for
            # each group of exposures
            sdriz = resample.ResampleData(self.input_models, single=True,
                                          blendheaders=False, **pars)
            sdriz.do_drizzle()
            drizzled_models = sdriz.output_models
            for model in drizzled_models:
                if save_intermediate_results:
                    log.info("Writing out resampled exposures...")
                    model_output_path = self.make_output_path(
                        basepath=model.meta.filename,
                        suffix='outlier_i2d')
                    model.save(model_output_path)
        else:
            drizzled_models = self.input_models
            for i in range(len(self.input_models)):
                drizzled_models[i].wht = build_driz_weight(
                    self.input_models[i],
                    weight_type='exptime',
                    good_bits=pars['good_bits'])

        # Initialize intermediate products used in the outlier detection
        median_model = datamodels.ImageModel(
                                        init=drizzled_models[0].data.shape)
        median_model.update(drizzled_models[0])
        median_model.meta.wcs = drizzled_models[0].meta.wcs

        # Perform median combination on set of drizzled mosaics
        median_model.data = self.create_median(drizzled_models)

        if save_intermediate_results:
            median_output_path = self.make_output_path(
                basepath=self.input_models[0].meta.filename,
                suffix='median'
            )
            log.info("Writing out MEDIAN image to: {}".format(
                median_output_path
            ))
            median_model.save(median_output_path)

        if pars['resample_data']:
            # Blot the median image back to recreate each input image specified
            # in the original input list/ASN/ModelContainer
            blot_models = self.blot_median(median_model)
            if save_intermediate_results:
                log.info("Writing out BLOT images...")
                for model in blot_models:
                    model_path = self.make_output_path(
                        basename=model.meta.filename,
                        suffix='blot'
                    )
                    model.save(model_path)
        else:
            # Median image will serve as blot image
            blot_models = datamodels.ModelContainer()
            for i in range(len(self.input_models)):
                blot_models.append(median_model)

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        self.detect_outliers(blot_models)

        # clean-up (just to be explicit about being finished with
        # these results)
        del median_model, blot_models

    def create_median(self, resampled_models):
        """Create a median image from the singly resampled images.

        NOTES
        -----
        This version is simplified from astrodrizzle's version in the
        following ways:
        - type of combination: fixed to 'median'
        - 'minmed' not implemented as an option
        """
        resampled_sci = [i.data.copy() for i in resampled_models]
        resampled_weight = [i.wht for i in resampled_models]

        maskpt = self.outlierpars.get('maskpt', 0.7)

        # Create a mask for each input image, masking out areas where there is
        # no data or the data has very low weight
        badmasks = []
        for weight in resampled_weight:
            # Create boolean masks for weight being zero or NaN
            mask_zero_weight = np.equal(weight, 0.)
            mask_nans = np.isnan(weight)
            # Combine the masks
            weight_masked = np.ma.array(weight, mask=np.logical_or(
                mask_zero_weight, mask_nans))
            # Sigma-clip the unmasked data
            weight_masked = sigma_clip(weight_masked, sigma=3, maxiters=5)
            mean_weight = np.mean(weight_masked)
            # Mask pixels where weight falls below maskpt percent
            weight_threshold = mean_weight * maskpt
            badmask = np.less(weight, weight_threshold)
            log.debug("Percentage of pixels with low weight: {}".format(
                np.sum(badmask) / len(weight.flat) * 100))
            badmasks.append(badmask)

        # Fill resampled_sci array with nan's where mask values are True
        #pdb.set_trace()
        for f1, f2 in zip(resampled_sci, badmasks):
            for elem1, elem2 in zip(f1, f2):
                elem1[elem2] = np.nan

        # For a of stack of images with "bad" data replaced with Nan
        # use np.nanmedian to compute the median.
        log.info("Generating median from {} images".format(len(resampled_sci)))
        median_image = np.nanmedian(resampled_sci, axis = 0)
        del resampled_sci

        return median_image

    def blot_median(self, median_model):
        """Blot resampled median image back to the detector images."""
        interp = self.outlierpars.get('interp', 'poly5')
        sinscl = self.outlierpars.get('sinscl', 1.0)

        # Initialize container for output blot images
        blot_models = datamodels.ModelContainer()

        log.info("Blotting median...")

        for model in self.input_models:
            blotted_median = model.copy()
            blot_root = '_'.join(model.meta.filename.replace(
                '.fits', '').split('_')[:-1])
            blotted_median.meta.filename = '{}_blot.fits'.format(blot_root)

            # clean out extra data not related to blot result
            blotted_median.err = None
            blotted_median.dq = None
            # apply blot to re-create model.data from median image
            blotted_median.data = gwcs_blot(median_model, model, interp=interp,
                                                     sinscl=sinscl)
            blot_models.append(blotted_median)

        return blot_models

    def detect_outliers(self, blot_models):
        """Flag DQ array for cosmic rays in input images.

        The science frame in each ImageModel in input_models is compared to
        the corresponding blotted median image in blot_models.  The result is
        an updated DQ array in each ImageModel in input_models.

        Parameters
        ----------
        input_models: JWST ModelContainer object
            data model container holding science ImageModels, modified in place

        blot_models : JWST ModelContainer object
            data model container holding ImageModels of the median output frame
            blotted back to the wcs and frame of the ImageModels in
            input_models

        Returns
        -------
        None
            The dq array in each input model is modified in place

        """

        for image, blot in zip(self.input_models, blot_models):
            image_copy = image.copy()
            flag_cr(image_copy, blot, **self.outlierpars)
            image.dq = image_copy.dq.copy()

        if self.converted:
            # Make sure actual input gets updated with new results
            for i in range(len(self.input_models)):
                self.inputs.dq[i, :, :] = self.input_models[i].dq


def flag_cr(sci_image, blot_image, **pars):
    """Masks outliers in science image by updating DQ in-place

    Mask blemishes in dithered data by comparing a science image
    with a model image and the derivative of the model image.

    Parameters
    ----------
    sci_image : ImageModel
        the science data

    blot_image : ImageModel
        the blotted median image of the dithered science frames

    pars : dict
        the user parameters for Outlier Detection

    Default parameters:

    grow     = 1               # Radius to mask [default=1 for 3x3]
    ctegrow  = 0               # Length of CTE correction to be applied
    snr      = "5.0 4.0"       # Signal-to-noise ratio
    scale    = "1.2 0.7"       # scaling factor applied to the derivative
    backg    = 0               # Background value

    """
    grow = pars.get('grow', 1)
    ctegrow = pars.get('ctegrow', 0)  # not provided by outlierpars
    backg = pars.get('backg', 0)
    snr1, snr2 = [float(val) for val in pars.get('snr', '5.0 4.0').split()]
    scl1, scl2 = [float(val) for val in pars.get('scale', '1.2 0.7').split()]

    # Get background level of science data if it has not been subtracted, so it
    # can be added into the level of the blotted data, which has been
    # background-subtracted
    if (sci_image.meta.background.subtracted is False and
        sci_image.meta.background.level is not None):
        subtracted_background = sci_image.meta.background.level
        log.debug(f"Adding background level {subtracted_background} to blotted image")
    else:
        # No subtracted background.  Allow user-set value, which defaults to 0
        subtracted_background = backg

    exptime = sci_image.meta.exposure.exposure_time

    sci_data = sci_image.data * exptime
    blot_data = (blot_image.data + subtracted_background) * exptime
    blot_deriv = abs_deriv(blot_data)

    err_data = np.nan_to_num(sci_image.err)

    # Define output cosmic ray mask to populate
    cr_mask = np.zeros(sci_image.shape, dtype=np.uint8)

    #
    #
    #    COMPUTATION PART I
    #
    #
    # Model the noise and create a CR mask
    diff_noise = np.abs(sci_data - blot_data)
    ta = np.sqrt(np.abs(blot_data) + err_data ** 2)
    t2 = scl1 * blot_deriv + snr1 * ta
    tmp1 = np.logical_not(np.greater(diff_noise, t2))

    # Convolve mask with 3x3 kernel
    kernel = np.ones((3, 3), dtype=np.uint8)
    tmp2 = np.zeros(tmp1.shape, dtype=np.int32)
    ndimage.convolve(tmp1, kernel, output=tmp2, mode='nearest', cval=0)
    #
    #
    #    COMPUTATION PART II
    #
    #
    # Create a second CR Mask
    xt2 = scl2 * blot_deriv + snr2 * ta

    np.logical_not(np.greater(diff_noise, xt2) & np.less(tmp2, 9), cr_mask)

    #
    #
    #    COMPUTATION PART III
    #
    #
    # Flag additional cte 'radial' and 'tail' pixels surrounding CR
    # pixels as CRs

    # In both the 'radial' and 'length' kernels below, 0=good and
    # 1=bad, so that upon convolving the kernels with cr_mask, the
    # convolution output will have low->bad and high->good from which
    # 2 new arrays are created having 0->bad and 1->good. These 2 new
    # arrays are then AND'ed to create a new cr_mask.

    # recast cr_mask to int for manipulations below; will recast to
    # Bool at end
    cr_mask_orig_bool = cr_mask.copy()
    cr_mask = cr_mask_orig_bool.astype(np.int8)

    # make radial convolution kernel and convolve it with original cr_mask
    cr_grow_kernel = np.ones((grow, grow))
    cr_grow_kernel_conv = cr_mask.copy()
    ndimage.convolve(cr_mask, cr_grow_kernel, output=cr_grow_kernel_conv)

    # make tail convolution kernel and (shortly) convolve it with
    # original cr_mask
    cr_ctegrow_kernel = np.zeros((2 * ctegrow + 1, 2 * ctegrow + 1))
    cr_ctegrow_kernel_conv = cr_mask.copy()

    # which pixels are masked by tail kernel depends on readout direction
    # We could put useful info in here for CTE masking if needed.  Code
    # remains below.  For now, we set to zero, which turns off CTE masking.
    ctedir = 0
    if (ctedir == 1):
        cr_ctegrow_kernel[0:ctegrow, ctegrow] = 1
    if (ctedir == -1):
        cr_ctegrow_kernel[ctegrow + 1:2 * ctegrow + 1, ctegrow] = 1
    if (ctedir == 0):
        pass

    # finally do the tail convolution
    ndimage.convolve(cr_mask, cr_ctegrow_kernel, output=cr_ctegrow_kernel_conv)

    # select high pixels from both convolution outputs; then 'and' them to
    # create new cr_mask
    where_cr_grow_kernel_conv = np.where(cr_grow_kernel_conv < grow * grow,
                                         0, 1)
    where_cr_ctegrow_kernel_conv = np.where(cr_ctegrow_kernel_conv < ctegrow,
                                            0, 1)

    # combine masks and cast back to Bool
    np.logical_and(where_cr_ctegrow_kernel_conv,
                   where_cr_grow_kernel_conv, cr_mask)
    cr_mask = cr_mask.astype(bool)

    count_sci = np.count_nonzero(sci_image.dq)
    count_cr = np.count_nonzero(cr_mask)
    log.debug("Pixels in input DQ: {}".format(count_sci))
    log.debug("Pixels in cr_mask:  {}".format(count_cr))

    # Update the DQ array in the input image.
    sci_image.dq = np.bitwise_or(sci_image.dq, np.invert(cr_mask) * CRBIT)


def abs_deriv(array):
    """Take the absolute derivate of a numpy array."""
    tmp = np.zeros(array.shape, dtype=np.float64)
    out = np.zeros(array.shape, dtype=np.float64)

    tmp[1:, :] = array[:-1, :]
    tmp, out = _absolute_subtract(array, tmp, out)
    tmp[:-1, :] = array[1:, :]
    tmp, out = _absolute_subtract(array, tmp, out)

    tmp[:, 1:] = array[:, :-1]
    tmp, out = _absolute_subtract(array, tmp, out)
    tmp[:, :-1] = array[:, 1:]
    tmp, out = _absolute_subtract(array, tmp, out)

    return out


def _absolute_subtract(array, tmp, out):
    tmp = np.abs(array - tmp)
    out = np.maximum(tmp, out)
    tmp = tmp * 0.
    return tmp, out


def gwcs_blot(median_model, blot_img, interp='poly5', sinscl=1.0):
    """
    Resample the output/resampled image to recreate an input image based on
    the input image's world coordinate system

    Parameters
    ----------
    median_model : ~jwst.datamodels.DataModel

    blot_img : datamodel
        Datamodel containing header and WCS to define the 'blotted' image

    interp : str, optional
        The type of interpolation used in the resampling. The
        possible values are "nearest" (nearest neighbor interpolation),
        "linear" (bilinear interpolation), "poly3" (cubic polynomial
        interpolation), "poly5" (quintic polynomial interpolation),
        "sinc" (sinc interpolation), "lan3" (3rd order Lanczos
        interpolation), and "lan5" (5th order Lanczos interpolation).

    sincscl : float, optional
        The scaling factor for sinc interpolation.
    """
    blot_wcs = blot_img.meta.wcs

    # Compute the mapping between the input and output pixel coordinates
    pixmap = calc_gwcs_pixmap(blot_wcs, median_model.meta.wcs, blot_img.data.shape)
    log.debug("Pixmap shape: {}".format(pixmap[:, :, 0].shape))
    log.debug("Sci shape: {}".format(blot_img.data.shape))

    pix_ratio = 1
    log.info('Blotting {} <-- {}'.format(blot_img.data.shape, median_model.data.shape))

    outsci = np.zeros(blot_img.shape, dtype=np.float32)

    # Currently tblot cannot handle nans in the pixmap, so we need to give some
    # other value.  -1 is not optimal and may have side effects.  But this is
    # what we've been doing up until now, so more investigation is needed
    # before a change is made.  Preferably, fix tblot in drizzle.
    pixmap[np.isnan(pixmap)] = -1
    tblot(median_model.data, pixmap, outsci, scale=pix_ratio, kscale=1.0,
                   interp=interp, exptime=1.0, misval=0.0, sinscl=sinscl)

    return outsci
