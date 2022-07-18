#!/usr/bin/env python
from collections import defaultdict
import os.path as op
import numpy as np

from .. import datamodels
from ..associations.lib.rules_level3_base import format_product
from ..exp_to_source import multislit_to_container
from ..master_background.master_background_step import split_container
from ..stpipe import Pipeline
from ..lib.exposure_types import is_moving_target

# step imports
from ..assign_mtwcs import assign_mtwcs_step
from ..cube_build import cube_build_step
from ..extract_1d import extract_1d_step
from ..master_background import master_background_step
from ..mrs_imatch import mrs_imatch_step
from ..outlier_detection import outlier_detection_step
from ..resample import resample_spec_step
from ..combine_1d import combine_1d_step
from ..photom import photom_step


__all__ = ['Spec3Pipeline']

# Group exposure types
MULTISOURCE_MODELS = ['MultiSlitModel']
IFU_EXPTYPES = ['MIR_MRS', 'NRS_IFU']
SLITLESS_TYPES = ['NIS_SOSS', 'NIS_WFSS', 'NRC_WFSS']


class Spec3Pipeline(Pipeline):
    """
    Spec3Pipeline: Processes JWST spectroscopic exposures from Level 2b to 3.

    Included steps are:
    assign moving target wcs (assign_mtwcs)
    master background subtraction (master_background)
    MIRI MRS background matching (mrs_imatch)
    outlier detection (outlier_detection)
    2-D spectroscopic resampling (resample_spec)
    3-D spectroscopic resampling (cube_build)
    1-D spectral extraction (extract_1d)
    Absolute Photometric Calibration (photom)
    1-D spectral combination (combine_1d)
    """

    class_alias = "calwebb_spec3"

    spec = """
    """

    # Define aliases to steps
    step_defs = {
        'assign_mtwcs': assign_mtwcs_step.AssignMTWcsStep,
        'master_background': master_background_step.MasterBackgroundStep,
        'mrs_imatch': mrs_imatch_step.MRSIMatchStep,
        'outlier_detection': outlier_detection_step.OutlierDetectionStep,
        'resample_spec': resample_spec_step.ResampleSpecStep,
        'cube_build': cube_build_step.CubeBuildStep,
        'extract_1d': extract_1d_step.Extract1dStep,
        'photom': photom_step.PhotomStep,
        'combine_1d': combine_1d_step.Combine1dStep
    }

    # Main processing
    def process(self, input):
        """Entrypoint for this pipeline

        Parameters
        ----------
        input: str, Level3 Association, or ~jwst.datamodels.DataModel
            The exposure or association of exposures to process
        """
        self.log.info('Starting calwebb_spec3 ...')
        asn_exptypes = ['science', 'background']

        # Setup sub-step defaults
        self.master_background.suffix = 'mbsub'
        self.mrs_imatch.suffix = 'mrs_imatch'
        self.outlier_detection.suffix = 'crf'
        self.outlier_detection.save_results = self.save_results
        self.resample_spec.suffix = 's2d'
        self.resample_spec.save_results = self.save_results
        self.cube_build.suffix = 's3d'
        self.cube_build.save_results = self.save_results
        self.extract_1d.suffix = 'x1d'
        self.extract_1d.save_results = self.save_results
        self.combine_1d.suffix = 'c1d'
        self.combine_1d.save_results = self.save_results

        # Retrieve the inputs:
        # could either be done via LoadAsAssociation and then manually
        # load input members into models and ModelContainer, or just
        # do a direct open of all members in ASN file, e.g.
        input_models = datamodels.open(input, asn_exptypes=asn_exptypes)

        # Immediately update the ASNTABLE keyword value in all inputs,
        # so that all outputs get the new value
        for model in input_models:
            model.meta.asn.table_name = op.basename(input_models.asn_table_name)

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

        if is_moving_target(input_models):
            self.log.info("Assigning WCS to a Moving Target exposure.")
            input_models = self.assign_mtwcs(input_models)

        # If background data are present, call the master background step
        if members_by_type['background']:
            source_models = self.master_background(input_models)
            source_models.meta.asn_table = input_models.meta.asn_table

            # If the step is skipped, do the container splitting that
            # would've been done in master_background
            if self.master_background.skip:
                source_models, bkg_models = split_container(input_models)
                # we don't need the background members
                bkg_models.close()
                del bkg_models
        else:
            # The input didn't contain any background members,
            # so we use all the inputs in subsequent steps
            source_models = input_models

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
        sources = [source_models]
        if model_type in MULTISOURCE_MODELS:
            self.log.info('Convert from exposure-based to source-based data.')
            sources = [
                (name, model)
                for name, model in multislit_to_container(source_models).items()
            ]

            # Check for negative and large source_id values
            if len(sources) > 99999:
                self.log.critical("Data contain more than 100,000 sources;"
                                  "filename does not support 6 digit source ids.")
                raise Exception

            available_src_ids = set(np.arange(99999) + 1)
            used_src_ids = set()
            for src in sources:
                src_id, model = src
                src_id = int(src_id)
                used_src_ids.add(src_id)
                if 0 < src_id <= 99999:
                    available_src_ids.remove(src_id)

            hotfixed_sources = []
            # now find and reset bad source_id values
            for src in sources:
                src_id, model = src
                src_id = int(src_id)
                # Replace ids that aren't positive 5-digit integers
                if src_id < 0 or src_id > 99999:
                    src_id_new = available_src_ids.pop()
                    self.log.info(f"Source ID {src_id} falls outside allowed range.")
                    self.log.info(f"Reassigning {src_id} to {str(src_id_new).zfill(5)}.")
                    # Replace source_id for each model in the SourceModelContainers
                    for contained_model in model:
                        contained_model.source_id = src_id_new
                    src_id = src_id_new
                hotfixed_sources.append((str(src_id), model))

            sources = hotfixed_sources

        # Process each source
        for source in sources:

            # If each source is a SourceModelContainer
            # the output name needs to be updated with the source name.
            if isinstance(source, tuple):
                source_id, result = source
                self.output_file = format_product(
                    output_file, source_id=source_id.lower()
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
            if exptype not in SLITLESS_TYPES:
                # Update the asn table name to the level 3 instance so that
                # the downstream products have the correct table name since
                # the _cal files are not saved they will not be updated
                for cal_array in result:
                    cal_array.meta.asn.table_name = op.basename(input_models.asn_table_name)
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
                    # instead run photom on the extract_1d results and save
                    # those instead.
                    self.extract_1d.save_results = False
                    result = self.extract_1d(result)

                    # SOSS F277W may return None - don't bother with that.
                    if result is not None:
                        self.photom.save_results = self.save_results
                        self.photom.suffix = 'x1d'
                        result = self.photom(result)
                else:
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

        input_models.close()

        self.log.info('Ending calwebb_spec3')
        return
