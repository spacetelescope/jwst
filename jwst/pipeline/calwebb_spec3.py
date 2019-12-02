#!/usr/bin/env python
from collections import defaultdict

from .. import datamodels
from ..associations.lib.rules_level3_base import format_product
from ..exp_to_source import multislit_to_container
from ..master_background.master_background_step import split_container
from ..stpipe import Pipeline

# step imports
from ..cube_build import cube_build_step
from ..extract_1d import extract_1d_step
from ..master_background import master_background_step
from ..mrs_imatch import mrs_imatch_step
from ..outlier_detection import outlier_detection_step
from ..resample import resample_spec_step
from ..combine_1d import combine_1d_step


__all__ = ['Spec3Pipeline']

# Group exposure types
MULTISOURCE_MODELS = ['MultiSlitModel']
IFU_EXPTYPES = ['MIR_MRS', 'NRS_IFU']
SLITLESS_TYPES = ['NIS_SOSS', 'NIS_WFSS', 'NRC_WFSS']


class Spec3Pipeline(Pipeline):
    """
    Spec3Pipeline: Processes JWST spectroscopic exposures from Level 2b to 3.

    Included steps are:
    master background subtraction (master_background)
    MIRI MRS background matching (mrs_imatch)
    outlier detection (outlier_detection)
    2-D spectroscopic resampling (resample_spec)
    3-D spectroscopic resampling (cube_build)
    1-D spectral extraction (extract_1d)
    1-D spectral combination (combine_1d)
    """

    spec = """
    """

    # Define aliases to steps
    step_defs = {
        'master_background': master_background_step.MasterBackgroundStep,
        'mrs_imatch': mrs_imatch_step.MRSIMatchStep,
        'outlier_detection': outlier_detection_step.OutlierDetectionStep,
        'resample_spec': resample_spec_step.ResampleSpecStep,
        'cube_build': cube_build_step.CubeBuildStep,
        'extract_1d': extract_1d_step.Extract1dStep,
        'combine_1d': combine_1d_step.Combine1dStep
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
        asn_exptypes = ['science', 'background']

        # Retrieve the inputs:
        # could either be done via LoadAsAssociation and then manually
        # load input members into models and ModelContainer, or just
        # do a direct open of all members in ASN file, e.g.
        input_models = datamodels.open(input, asn_exptypes=asn_exptypes)

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

        # Find all the member types in the product
        members_by_type = defaultdict(list)
        product = input_models.meta.asn_table.products[0].instance
        for member in product['members']:
            members_by_type[member['exptype'].lower()].append(member['expname'])

        # If background data are present, call the master background step
        if members_by_type['background']:
            source_models = self.master_background(input_models)
            source_models.meta.asn_table = input_models.meta.asn_table

            # If the step is skipped, do the container splitting that
            # would've been done in master_background
            if self.master_background.skip:
                source_models, bkg_models = split_container(input_models)
                del bkg_models  # we don't need the background members
        else:
            # The input didn't contain any background members,
            # so we use all the inputs in subsequent steps
            source_models = input_models

        # For JWST spectral modes, the input associations can contain
        # one of two types of collections: Single Object or Multi-Object.
        #
        # For Single Object, sometimes referred to as "target-based", data, the input
        # list of exposures contain data for just the one object.
        #
        # For Multi-Object, each exposure in the input list of exposures contain data
        # for all the objects.
        #
        # For uniformity of processing, the input is arranged by source.
        # For Single Object, the source id is simply the target id.
        # For Multi-Object, the data is rearranged from the exposure-centric models,
        # to a source-centric model where each model contains the part of the exposure
        # only relevant to that source.
        #
        # The rest of the processing is then done on a source-by-source basis.

        # Arrange the data in a source-based hierarchy.
        # `sources` will be a 2-tuple consisting of:
        #    (id, ModelContainer)
        #
        #     `id`: A `str` with the source name. If single-source,
        #     then this will have the target identifier
        #     `ModelContainer`: The list of data belonging to the source.
        #
        # If target-based, `id` will be `target`
        if model_type in MULTISOURCE_MODELS:
            # Multi-source information. Invert the data structure.
            self.log.info('Convert from exposure-based to source-based data.')
            sources = [
                (name, model)
                    for name, model in multislit_to_container(source_models).items()
                ]
        else:
            # Single-source. The source ID is simply the target name.
            sources = [('target', source_models)]

        # Process each source
        for source in sources:
            source_id, result = source

            # If multi-object data, reformat the output name.
            if source_id != 'target':
                self.output_file = format_product(
                    output_file, source_id=source_id.lower()
                )

            # The MultiExposureModel is a required output.
            if isinstance(result, datamodels.SourceModelContainer):
                self.save_model(result, 'cal')

            # Call the skymatch step for MIRI MRS data
            if exptype in ['MIR_MRS']:
                result = self.mrs_imatch(result)

            # Call outlier detection
            if exptype not in SLITLESS_TYPES:
                # Update the asn table name to the level 3 instance so that
                # the downstream products have the correct table name since
                # the _cal files are not saved they will not be updated
                for cal_array in result:
                    cal_array.meta.asn.table_name = result.meta.table_name
                result = self.outlier_detection(result)

                # Resample time. Dependent on whether the data is IFU or not.
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
            if exptype in SLITLESS_TYPES:

                # For slitless data, extract 1D spectra and then combine them

                if exptype in ['NIS_SOSS']:
                    # For NIRISS SOSS, don't save the extract_1d results,
                    # they're identical to the calwebb_spec2 x1d products
                    self.extract_1d.save_results = False

                result = self.extract_1d(result)
                result = self.combine_1d(result)

            elif resample_complete is not None and resample_complete.upper() == 'COMPLETE':

                # If 2D data were resampled and combined, just do a 1D extraction
                if exptype in IFU_EXPTYPES:
                    self.extract_1d.search_output_file = False
                result = self.extract_1d(result)

            else:
                self.log.warning(
                    'Resampling was not completed. Skipping extract_1d.'
                )

        # We're done
        self.log.info('Ending calwebb_spec3')
        return
