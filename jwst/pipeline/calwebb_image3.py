from ..stpipe import Pipeline
from .. import datamodels
from ..lib.exposure_types import is_moving_target

from ..assign_mtwcs import assign_mtwcs_step
from ..tweakreg import tweakreg_step
from ..skymatch import skymatch_step
from ..resample import resample_step
from ..outlier_detection import outlier_detection_step
from ..source_catalog import source_catalog_step

__all__ = ['Image3Pipeline']


class Image3Pipeline(Pipeline):
    """
    Image3Pipeline: Applies level 3 processing to imaging-mode data from
                    any JWST instrument.

    Included steps are:
        assign_mtwcs
        tweakreg
        skymatch
        outlier_detection
        resample
        source_catalog
    """

    class_alias = "calwebb_image3"

    spec = """
    """

    # Define alias to steps
    step_defs = {
        'assign_mtwcs': assign_mtwcs_step.AssignMTWcsStep,
        'tweakreg': tweakreg_step.TweakRegStep,
        'skymatch': skymatch_step.SkyMatchStep,
        'outlier_detection': outlier_detection_step.OutlierDetectionStep,
        'resample': resample_step.ResampleStep,
        'source_catalog': source_catalog_step.SourceCatalogStep
    }

    def process(self, input_data):
        """
        Run the Image3Pipeline

        Parameters
        ----------
        input_data: Level3 Association, or ~jwst.datamodels.ModelContainer
            The exposures to process
        """
        self.log.info('Starting calwebb_image3 ...')

        # Only load science members from input ASN;
        # background and target-acq members are not needed.
        asn_exptypes = ['science']

        # Configure settings for saving results files
        self.outlier_detection.suffix = 'crf'
        self.outlier_detection.save_results = self.save_results

        self.resample.suffix = 'i2d'
        self.resample.save_results = self.save_results

        self.source_catalog.save_results = self.save_results

        with datamodels.open(input_data, asn_exptypes=asn_exptypes) as input_models:
            # If input is an association, set the output to the product name.
            if self.output_file is None:
                try:
                    self.output_file = input_models.meta.asn_table.products[0].name
                except (AttributeError, IndexError):
                    pass

            # Check if input is single or multiple exposures
            try:
                has_groups = len(input_models.group_names) >= 1
            except (AttributeError, TypeError, KeyError):
                has_groups = False

            if isinstance(input_models, datamodels.ModelContainer) and has_groups:
                if is_moving_target(input_models):
                    input_models = self.assign_mtwcs(input_models)
                else:
                    input_models = self.tweakreg(input_models)

                input_models = self.skymatch(input_models)
                input_models = self.outlier_detection(input_models)

            elif self.skymatch.skymethod == 'match':
                self.log.warning("Turning 'skymatch' step off for a single "
                                 "input image when 'skymethod' is 'match'")

            else:
                input_models = self.skymatch(input_models)

            result = self.resample(input_models)
            if isinstance(result, datamodels.ImageModel) and result.meta.cal_step.resample == 'COMPLETE':
                self.source_catalog(result)
