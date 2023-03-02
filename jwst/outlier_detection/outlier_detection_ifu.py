"""Class definition for performing outlier detection on IFU data."""

import numpy as np
from stsci.image import median
from astropy.stats import sigma_clipped_stats

from stdatamodels.jwst import datamodels

from .outlier_detection import OutlierDetection
from ..cube_build.cube_build_step import CubeBuildStep
from ..cube_build import blot_cube_build


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["OutlierDetectionIFU"]


class OutlierDetectionIFU(OutlierDetection):
    """Sub-class defined for performing outlier detection on IFU data.

    This is the controlling routine for the outlier detection process.
    It loads and sets the various input data and parameters needed by
    the various functions and then controls the operation of this process
    through all the steps used for the detection.

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input ModelContainer and merges
         them with any user-provided values
      2. Resamples all input images into IFUCubeModel observations.
      3. Creates a median image from all IFUCubeModels.
      4. Blot median image using CubeBlot to match
          each original input ImageModel.
      5. Perform statistical comparison between blotted image and original
         image to identify outliers.
      6. Updates input ImageModel DQ arrays with mask of detected outliers.

    """

    default_suffix = 's3d'

    def __init__(self, input_models, reffiles=None, **pars):
        """Initialize class for IFU data processing.

        Parameters
        ----------
        input_models : ~jwst.datamodels.ModelContainer, str
            list of data models as ModelContainer or ASN file,
            one data model for each input 2-D ImageModel

        drizzled_models : list of objects
            ModelContainer containing drizzled grouped input images

        reffiles : dict of ~jwst.datamodels.DataModel
            Dictionary of datamodels.  Keys are reffile_types.


        """
        OutlierDetection.__init__(self, input_models,
                                  reffiles=reffiles, **pars)

    def _find_ifu_coverage(self):
        self.channels = []
        self.subchannels = []
        self.band_name = []
        self.gratings = []
        self.ifu_band1 = []
        self.ifu_band2 = []
        self.instrument = self.input_models[0].meta.instrument.name.upper()
        n = len(self.input_models)
        for i in range(n):
            if self.instrument == 'MIRI':
                this_channel = self.input_models[i].meta.instrument.channel
                this_subchannel = self.input_models[i].meta.instrument.band.lower()
                nc = len(this_channel)
                for k in range(nc):
                    self.band_name.append('ch' + this_channel[k] + '_' + this_subchannel)
            elif self.instrument == 'NIRSPEC':
                self.gratings.append(self.input_models[i].meta.instrument.grating.lower())
            else:
                # add error
                raise ErrorWrongInstrument('Instrument must be MIRI or NIRSPEC')
        if self.instrument == 'MIRI':
            band_no_repeat = list(set(self.band_name))
            bands = ['short','medium','long','short-medium','short-long','medium-short','medium-long','long-short','long-medium']
            channels = ['1','2','3','4']
            for this_band in band_no_repeat:
                for j in channels:
                    check_channel = 'ch' + j
                    if check_channel in this_band:
                        self.channels.append(j)
                for k in bands:
                    compare_band = this_band[4:]  # remove the channel part of the name to see which band we have
                    if k == compare_band:
                        self.subchannels.append(k)
            self.ifu_band1 = self.channels
            self.ifu_band2 = self.subchannels
        elif self.instrument == 'NIRSPEC':
            self.gratings = list(set(self.gratings))
            self.ifu_band1 = self.gratings
            self.ifu_band2 = self.gratings  # not used in NIRSpec

    def _convert_inputs(self):
        self.input_models = self.inputs
        self.converted = False

    def do_detection(self):
        """Flag outlier pixels in DQ of input images."""
        self._convert_inputs()
        self._find_ifu_coverage()

        self.build_suffix(**self.outlierpars)

        save_intermediate_results = \
            self.outlierpars['save_intermediate_results']

        # start by creating copies of the input data to place the separate
        # data in after blotting the median-combined cubes for each channel
        self.blot_models = self.inputs.copy()
        for model in self.blot_models:
            # replace arrays with all zeros to accommodate blotted data
            model.data = np.zeros(model.data.shape, dtype=model.data.dtype)

        # Create the resampled/mosaic images for each group of exposures
        #
        exptype = self.input_models[0].meta.exposure.type
        log.info("Performing IFU outlier_detection for exptype {}".format(
                 exptype))
        num_bands = len(self.ifu_band1)
        for i in range(num_bands):
            select1 = self.ifu_band1[i]
            select2 = self.ifu_band2[i]

            if self.instrument == 'MIRI':
                cubestep = CubeBuildStep(channel=select1, band=select2, weighting='emsm',
                                         single=True)

            if self.instrument == 'NIRSPEC':
                cubestep = CubeBuildStep(grating=select1, weighting='emsm',
                                         single=True)

            single_IFUCube_result = cubestep.process(self.input_models)

            for model in single_IFUCube_result:
                model.meta.filename = self.make_output_path(
                    basepath=model.meta.filename,
                    suffix='outlier_s3d'
                )
                if save_intermediate_results:
                    log.info("Writing out (single) IFU cube {}".format(model.meta.filename))
                    model.save(model.meta.filename)

            # Initialize intermediate products used in the outlier detection
            median_model = datamodels.IFUCubeModel(
                init=single_IFUCube_result[0].data.shape)
            median_model.meta = single_IFUCube_result[0].meta

            if self.instrument == 'MIRI':
                median_model.meta.filename = self.make_output_path(
                    basepath=self.input_models.meta.asn_table.products[0].name,
                    suffix='ch{}_{}_median_s3d'.format(select1,select2))
            else:
                median_model.meta.filename = self.make_output_path(
                    basepath=self.input_models.meta.asn_table.products[0].name,
                    suffix='{}_median_s3d'.format(select1))
            # Perform median combination on set of drizzled mosaics
            median_model.data = self.create_median(single_IFUCube_result)

            if save_intermediate_results:
                log.info("Writing out MEDIAN image to: {}".format(
                    median_model.meta.filename))
                median_model.save(median_model.meta.filename)

            # Blot the median image back to recreate each input image specified
            # in the original input list/ASN/ModelContainer
            #
            # need to override with IFU-specific version of blot for
            # each channel/grating this will need to combine the multiple
            # channels (MIRI) of data into a single frame to match the
            # original input...
            self.blot_median(median_model)

        if save_intermediate_results:
            log.info("Writing out BLOT images...")

            for model in self.blot_models:
                model.meta.filename = self.make_output_path(
                    basepath=model.meta.filename, suffix='blot')

                log.info("Blotted files {}".format(model.meta.filename))
                model.save(model.meta.filename)
        # Perform outlier detection using statistical comparisons between
        # each original input image and the blotted version of the
        # median image of all channels
        self.detect_outliers(self.blot_models)

        # clean-up (just to be explicit about being finished
        # with these results)
        self.blot_models = None
        del median_model

    def create_median(self, resampled_models):
        """IFU-specific version of create_median."""
        resampled_sci = [i.data for i in resampled_models]
        resampled_wht = [i.weightmap for i in resampled_models]
        nlow = self.outlierpars.get('nlow', 0)
        nhigh = self.outlierpars.get('nhigh', 0)
        maskpt = self.outlierpars.get('maskpt', 0.7)
        badmasks = []
        for w in resampled_wht:
            # Due to a bug in numpy.nanmean, need to check
            # for a completely zero array
            if not np.any(w):
                mean_weight = 0.
            else:
                mean_weight, _, _ = sigma_clipped_stats(
                    w, sigma=3.0, mask_value=0.
                )
            weight_threshold = mean_weight * maskpt
            # Mask pixels were weight falls below MASKPT percent of
            #    the mean weight
            mask = np.less(w, weight_threshold)
            log.debug("Number of pixels with low weight: {}".format(
                np.sum(mask)))
            badmasks.append(mask)

        # Compute median of stack os images using BADMASKS to remove low weight
        # values
        median_image = median(resampled_sci, nlow=nlow, nhigh=nhigh,
                              badmasks=badmasks)
        # for i in range(len(resampled_sci)):
        #    data = resampled_sci[i]
        #    mask = badmasks[i]
        return median_image

    def blot_median(self, median_image):
        """IFU-specific version of blot_median."""
        cubeblot = blot_cube_build.CubeBlot(median_image, self.input_models)
        cubeblot.blot_info()
        blot_models, input_list_number = cubeblot.blot_images()
        for j in range(len(blot_models)):
            k = input_list_number[j]
            self.blot_models[k].data += blot_models[j].data
            self.blot_models[k].meta = blot_models[j].meta


class ErrorWrongInstrument(Exception):
    """ Raises an exception if the instrument is not MIRI or NIRSPEC
    """
    pass
