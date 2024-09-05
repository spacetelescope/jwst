from collections.abc import Sequence
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelLibrary

from ..stpipe import Pipeline
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
    in_memory = boolean(default=True)  # If False, preserve memory using temporary files at the expense of runtime
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
        input_data: Level3 Association, or ~jwst.datamodels.ModelLibrary
            The exposures to process
        """
        self.log.info('Starting calwebb_image3 ...')

        # Configure settings for saving results files
        self.outlier_detection.suffix = 'crf'
        self.outlier_detection.mode = 'imaging'
        self.outlier_detection.save_results = self.save_results

        self.resample.suffix = 'i2d'
        self.resample.save_results = self.save_results

        self.source_catalog.save_results = self.save_results

        # Only load science members from input ASN;
        # background and target-acq members are not needed.
        input_models = self._load_input_as_library(input_data)

        if (self.output_file is None) and ("name" in input_models.asn["products"][0]):
            # If input is an association, set the output to the product name.
            self.output_file = input_models.asn["products"][0]["name"]

        # Check if input is single or multiple exposures
        has_groups = len(input_models.group_names) >= 1

        if has_groups:
            with input_models:
                model = input_models.borrow(0)
                is_moving = is_moving_target(model)
                input_models.shelve(model, 0, modify=False)
            if is_moving:
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
        del input_models
        if isinstance(result, datamodels.ImageModel) and result.meta.cal_step.resample == 'COMPLETE':
            self.source_catalog(result)


    def _load_input_as_library(self, input):
        """
        Load any valid input type into a ModelLibrary, including
        single datamodels, associations, ModelLibrary instances, and
        filenames pointing to those types.
        """

        if isinstance(input, ModelLibrary):
            return input

        if isinstance(input, (str, dict)):
            try:
                # Try opening input as an association
                return ModelLibrary(input, asn_exptypes=['science'], on_disk=not self.in_memory)
            except OSError:
                # Try opening input as a single cal file
                input = datamodels.open(input)
                input = [input,]
                return ModelLibrary(input, asn_exptypes=['science'], on_disk=not self.in_memory)
        elif isinstance(input, Sequence):
            return ModelLibrary(input, asn_exptypes=['science'], on_disk=not self.in_memory)
        elif isinstance(input, datamodels.JwstDataModel):
            return ModelLibrary([input], asn_exptypes=['science'], on_disk=not self.in_memory)
        else:
            raise TypeError(f"Input type {type(input)} not supported.")