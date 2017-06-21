from __future__ import unicode_literals, absolute_import

import os

from ..stpipe import Pipeline
from .. import datamodels

from ..resample import resample_step
from ..skymatch import skymatch_step
from ..outlier_detection import outlier_detection_step
from ..source_catalog import source_catalog_step
from ..tweakreg_catalog import tweakreg_catalog_step
from ..tweakreg import tweakreg_step

__version__ = "0.7.0"


class Image3Pipeline(Pipeline):
    """
    Image3Pipeline: Applies level 3 processing to imaging-mode data from
                    any JWST instrument.

    Included steps are:
        tweakreg_catalog
        tweakreg
        skymatch
        outlier_detection
        resample
        source_catalog
    """

    spec = """
        suffix = string(default='i2d')
    """

    # Define alias to steps
    step_defs = {'tweakreg_catalog': tweakreg_catalog_step.TweakregCatalogStep,
                 'tweakreg': tweakreg_step.TweakRegStep,
                 'skymatch': skymatch_step.SkyMatchStep,
                 'outlier_detection': outlier_detection_step.OutlierDetectionStep,
                 'resample': resample_step.ResampleStep,
                 'source_catalog': source_catalog_step.SourceCatalogStep
                 }

    def process(self, input):
        """
        Run the Image3Pipeline

        Parameters
        ----------
        input: Level3 Association, or ModelContainer
            The exposures to process
        """

        self.log.info('Starting calwebb_image3 ...')

        input_models = datamodels.open(input)

        # Check if input is multiple exposures, as required by some steps
        is_container = isinstance(input_models, datamodels.ModelContainer)
        if is_container and len(input_models.group_names) > 1:

            self.log.info("Generating source catalogs for alignment...")
            input_models = self.tweakreg_catalog(input_models)

            self.log.info("Aligning input images...")
            input_models = self.tweakreg(input_models)

            # Clean up tweakreg catalogs which no are no longer needed
            for model in input_models:
                try:
                    catalog_name = model.meta.tweakreg_catalog.filename
                    os.remove(catalog_name)
                except:
                    pass

            self.log.info("Matching sky values across all input images...")
            input_models = self.skymatch(input_models)

            self.log.info("Performing outlier detection on input images...")
            input_models = self.outlier_detection(input_models)

            self.log.info("Writing Level 2c images with updated DQ arrays...")
            suffix_2c = 'cal-{}'.format(input_models.meta.asn_table.asn_id)
            for model in input_models:
                self.save_model(model, suffix=suffix_2c)

        self.log.info("Resampling images to final output...")
        output = self.resample(input_models)

        product = input_models.meta.asn_table.products[0].name + '.fits'
        output.meta.filename = product
        self.save_model(output, suffix=self.suffix)
        self.log.info('Saved resampled image to %s', output.meta.filename)

        self.log.info("Creating source catalog...")
        out_catalog = self.source_catalog(output)
        # NOTE: source_catalog step writes out the catalog in .ecsv format
        # In the future it would be nice if it was returned to the pipeline,
        # and then written here.  A datamodel for .ecsv might be required.

        return
