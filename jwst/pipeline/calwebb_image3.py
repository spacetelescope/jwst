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

__version__ = "0.7.1"


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

        # Check if input is single or multiple exposures
        is_container = isinstance(input_models, datamodels.ModelContainer)
        try:
            has_groups = len(input_models.group_names) > 1
        except:
            has_groups = False
        if is_container and has_groups:

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

            if input_models[0].meta.cal_step.outlier_detection == 'COMPLETE':
                self.log.info("Writing Level 2c images with updated DQ arrays...")
                # Set up Level 2c suffix to be used later
                asn_id = input_models.meta.asn_table.asn_id
                suffix_2c = '{}_{}'.format(asn_id, 'crf')
                for model in input_models:
                    self.save_model(model, suffix=suffix_2c)

        self.log.info("Resampling images to final result...")
        result = self.resample(input_models)

        try:
            result.meta.asn.pool_name = input_models.meta.asn_table.asn_pool
            result.meta.asn.table_name = input
        except:
            pass

        result.meta.filename = input_models.meta.asn_table.products[0].name
        self.save_model(result, suffix=self.suffix)

        self.log.info("Creating source catalog...")
        out_catalog = self.source_catalog(result)
        # NOTE: source_catalog step writes out the catalog in .ecsv format
        # In the future it would be nice if it was returned to the pipeline,
        # and then written here.  A datamodel for .ecsv might be required.

        return
