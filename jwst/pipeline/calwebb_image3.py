#!/usr/bin/env python
from jwst.stpipe import Pipeline
from .. import datamodels

# calwebb Image3 step imports
from jwst.resample import resample_step
from jwst.skymatch import skymatch_step
from jwst.outlier_detection import outlier_detection_step
from jwst.source_catalog import source_catalog_step
from jwst.tweakreg_catalog import tweakreg_catalog_step
from jwst.tweakreg import tweakreg_step

__version__ = "0.1"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class Image3Pipeline(Pipeline):
    """

    Image3Pipeline: Applies level 3 processing to imaging-mode data from
                    any JWST instrument.

    Included steps are:
    resample.

    """

    # Define alias to steps
    step_defs = {'resample': resample_step.ResampleStep,
                 'skymatch': skymatch_step.SkyMatchStep,
                 'outlier_detection': outlier_detection_step.OutlierDetectionStep,
                 'tweakreg': tweakreg_step.TweakRegStep,
                 'source_catalog': source_catalog_step.SourceCatalogStep,
                 'tweakreg_catalog': tweakreg_catalog_step.TweakregCatalogStep
                 }

    def process(self, input):

        log.info('Starting calwebb_image3 ...')

        input_models = datamodels.open(input)

        is_container = (type(input_models) == type(datamodels.ModelContainer()))
        if is_container and len(input_models.group_names) > 1:
            # perform full outlier_detection of ASN data
            log.info("Generating source catalogs for alignment...")
            input_models = self.tweakreg_catalog(input_models)
            log.info("Aligning input images...")
            input_models = self.tweakreg(input_models)
            log.info("Matching sky values across all input images...")
            input_models = self.skymatch(input_models)
            log.info("Performing outlier detection on input images...")
            input_models = self.outlier_detection(input_models)

            log.info("Resampling ASN to create combined product: {}".format(input_models.meta.resample.output))

        # Resample step always returns ModelContainer,
        # yet we only need the DataModel result
        output = self.resample(input_models)

        # create final source catalog from resampled output
        out_catalog = self.source_catalog(output)
        output.close()
        input_models.close()
        log.info('... ending calwebb_image3')

        return
