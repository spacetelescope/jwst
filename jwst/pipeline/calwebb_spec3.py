#!/usr/bin/env python

from .. import datamodels
from ..associations.lib.rules_level3_base import format_product
from ..exp_to_source import multislit_to_container
from ..stpipe import Pipeline

# step imports
from ..cube_build import cube_build_step
from ..extract_1d import extract_1d_step
from ..mrs_imatch import mrs_imatch_step
from ..outlier_detection import outlier_detection_step
from ..resample import resample_spec_step

__version__ = '0.9.3'

# Group exposure types
MULTISOURCE_MODELS = ['MultiSlitModel']
IFU_EXPTYPES = ['MIR_MRS', 'NRS_IFU']


class Spec3Pipeline(Pipeline):
    """
    Spec3Pipeline: Processes JWST spectroscopic exposures from Level 2b to 3.

    Included steps are:
    MIRI MRS background matching (skymatch)
    outlier detection (outlier_detection)
    2-D spectroscopic resampling (resample_spec)
    3-D spectroscopic resampling (cube_build)
    1-D spectral extraction (extract_1d)
    """

    spec = """
    """

    # Define aliases to steps
    step_defs = {
        'mrs_imatch': mrs_imatch_step.MRSIMatchStep,
        'outlier_detection': outlier_detection_step.OutlierDetectionStep,
        'resample_spec': resample_spec_step.ResampleSpecStep,
        'cube_build': cube_build_step.CubeBuildStep,
        'extract_1d': extract_1d_step.Extract1dStep
    }

    # Main processing
    def process(self, input):
        """Entrypoint for this pipeline

        Parameters
        ----------
        input: str, Level3 Association, or DataModel
            The exposure or association of exposures to process
        """
        self.log.info('Starting calwebb_spec3 ...')

        # Retrieve the inputs:
        # could either be done via LoadAsAssociation and then manually
        # load input members into models and ModelContainer, or just
        # do a direct open of all members in ASN file, e.g.
        input_models = datamodels.open(input)

        # For the first round of development we will assume that the input
        # is ALWAYS an ASN. There's no use case for anyone ever running a
        # single exposure through.

        # Once data are loaded, store a few things for future use;
        # some of this is here only for the purpose of creating fake
        # products until the individual tasks work and do it themselves
        exptype = input_models[0].meta.exposure.type
        model_type = input_models[0].meta.model_type
        output_file = input_models.meta.asn_table.products[0].name
        self.output_file = output_file

        # `sources` is the list of astronomical sources that need be
        # processed. Each element is a ModelContainer, which contains
        # models for all exposures that belong to a single source.
        #
        # For JWST spectral modes, the input associations can contain
        # one of two types of collections. If the exposure type is
        # considered single-source, then the association contains only
        # exposures of that source.
        #
        # However, there are modes in which the exposures contain data
        # from multiple sources. In that case, the data must be
        # rearranged, collecting the exposures representing each
        # source into its own ModelContainer. This produces a list of
        # sources, each represented by a MultiExposureModel instead of
        # a single ModelContainer.
        sources = [input_models]
        if model_type in MULTISOURCE_MODELS:
            self.log.info('Convert from exposure-based to source-based data.')
            sources = [
                (name, model)
                for name, model in multislit_to_container(input_models).items()
            ]

        # Process each source
        for source in sources:

            # If each source is a SourceModelContainer
            # the output name needs to be updated with the source name.
            if isinstance(source, tuple):
                source_id, result = source
                self.output_file = format_product(
                    output_file, source_id=int(source_id)
                )
            else:
                result = source

            # The MultiExposureModel is a required output.
            if isinstance(result, datamodels.SourceModelContainer):
                self.save_model(result, 'cal')

            # Call the skymatch step for MIRI MRS data
            if exptype in ['MIR_MRS']:
                result = self.mrs_imatch(result)

            # Call outlier detection
            result = self.outlier_detection(result)

            # Resample time. Dependent on whether the data is IFU or
            # not.
            resample_complete = None
            if exptype in IFU_EXPTYPES:
                result = self.cube_build(result)
                try:
                    resample_complete = result[0].meta.cal_step.cube_build
                except AttributeError:
                    pass
            else:
                result = self.resample_spec(result)
                try:
                    resample_complete = result.meta.cal_step.resample
                except AttributeError:
                    pass

            # Do 1-D spectral extraction
            if resample_complete is not None and \
               resample_complete.upper() == 'COMPLETE':
                if exptype in IFU_EXPTYPES:
                    self.extract_1d.search_output_file = False
                result = self.extract_1d(result)
            else:
                self.log.warn(
                    'Resampling was not completed. Skipping extract_1d.'
                )

        # We're done
        self.log.info('Ending calwebb_spec3')
        return
