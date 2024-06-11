#!/usr/bin/env python
from collections import defaultdict
from functools import wraps
import os.path as op

from stdatamodels.jwst import datamodels

from jwst.datamodels import SourceModelContainer

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
from ..spectral_leak import spectral_leak_step
from ..pixel_replace import pixel_replace_step

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
        'pixel_replace': pixel_replace_step.PixelReplaceStep,
        'resample_spec': resample_spec_step.ResampleSpecStep,
        'cube_build': cube_build_step.CubeBuildStep,
        'extract_1d': extract_1d_step.Extract1dStep,
        'photom': photom_step.PhotomStep,
        'combine_1d': combine_1d_step.Combine1dStep,
        'spectral_leak': spectral_leak_step.SpectralLeakStep
    }

    # Main processing
    def process(self, input):
        """Entrypoint for this pipeline

        Parameters
        ----------
        input: str, Level3 Association, or ~jwst.datamodels.JwstDataModel
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
        self.spectral_leak.suffix = 'x1d'
        self.spectral_leak.save_results = self.save_results
        self.pixel_replace.suffix = 'pixel_replace'
        self.pixel_replace.output_use_model = True
        self.pixel_replace.save_results = self.save_results
        
        # Overriding the Step.save_model method for the following steps.
        # These steps save intermediate files, resulting in meta.filename
        # being modified. This can affect the filenames of subsequent
        # steps.
        self.outlier_detection.save_model = invariant_filename(self.outlier_detection.save_model)
        self.pixel_replace.save_model = invariant_filename(self.pixel_replace.save_model)

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

        # Process each source
        for source in sources:

            # If each source is a SourceModelContainer,
            # the output name needs to be updated based on the source ID,
            # and potentially also the slit name (for NIRSpec fixed-slit only).
            if isinstance(source, tuple):
                source_id, result = source
                # NIRSpec fixed-slit data
                if exptype == "NRS_FIXEDSLIT":
                    # Output file name is constructed using the source_id and the slit name
                    slit_name = self._create_nrsfs_slit_name(result)
                    srcid = f's{source_id:>09s}'
                    self.output_file = format_product(output_file, source_id=srcid, slit_name=slit_name)

                # NIRSpec MOS/MSA data
                elif exptype == "NRS_MSASPEC":
                    # Construct the specially formatted source_id to use in the output file
                    # name that separates source, background, and virtual slits
                    srcid = self._create_nrsmos_source_id(result)
                    self.output_file = format_product(output_file, source_id=srcid)
                    self.log.debug(f"output_file = {self.output_file}")

                else:
                    # All other types just use the source_id directly in the file name
                    srcid = f's{source_id:>09s}'
                    self.output_file = format_product(output_file, source_id=srcid)
            else:
                result = source

            # The MultiExposureModel is a required output.
            if isinstance(result, SourceModelContainer):
                self.save_model(result, 'cal')

            # Call the skymatch step for MIRI MRS data
            if exptype in ['MIR_MRS']:
                result = self.mrs_imatch(result)

            # Call outlier detection and pixel replacement
            if exptype not in SLITLESS_TYPES:
                # Update the asn table name to the level 3 instance so that
                # the downstream products have the correct table name since
                # the _cal files are not saved they will not be updated
                for cal_array in result:
                    cal_array.meta.asn.table_name = op.basename(input_models.asn_table_name)
                result = self.outlier_detection(result)

                # interpolate pixels that have a NaN value or are flagged
                # as DO_NOT_USE or NON_SCIENCE.
                result = self.pixel_replace(result)

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

                # interpolate pixels that have a NaN value or are flagged
                # as DO_NOT_USE or NON_SCIENCE
                result = self.pixel_replace(result)

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
                    if exptype in ['MIR_MRS']:
                        if not self.spectral_leak.skip:
                            self.extract_1d.save_results = False
                            self.spectral_leak.suffix = 'x1d'
                            self.spectral_leak.search_output_file = False
                            self.spectral_leak.save_results = self.save_results

                result = self.extract_1d(result)

                if exptype in ['MIR_MRS']:
                    result = self.spectral_leak(result)
            else:
                self.log.warning(
                    'Resampling was not completed. Skipping extract_1d.'
                )

        input_models.close()

        self.log.info('Ending calwebb_spec3')
        return

    def _create_nrsfs_slit_name(self, source_models):
        """Create the complete slit_name product field for NIRSpec fixed-slit products

        Each unique value of slit name within the list of input source models
        is appended to the final slit name string.
        """

        slit_names = []
        slit_names.append(source_models[0].name.lower())
        for i in range(len(source_models)):
            name = source_models[i].name.lower()
            if name not in slit_names:
                slit_names.append(name)
        slit_name = "-".join(slit_names)  # append slit names using a dash separator

        return slit_name

    def _create_nrsmos_source_id(self, source_models):
        """Create the complete source_id product field for NIRSpec MOS products.

        The source_id value gets a "s", "b", or "v" character prepended
        to uniquely identify source, background, and virtual slits.
        """

        # Get the original source name and ID from the input models
        source_name = source_models[0].source_name
        source_id = source_models[0].source_id

        # MOS background sources have "BKG" in the source name
        if "BKG" in source_name:
            # prepend "b" to the source_id number and format to 9 chars
            srcid = f'b{str(source_id):>09s}'
            self.log.debug(f"Source {source_name} is a MOS background slitlet: ID={srcid}")

        # MOS virtual sources have a negative source_id value
        elif source_id < 0:
            # prepend "v" to the source_id number and remove the leading negative sign
            srcid = f'v{str(source_id)[1:]:>09s}'
            self.log.debug(f"Source {source_name} is a MOS virtual slitlet: ID={srcid}")

        # Regular MOS sources
        else:
            # prepend "s" to the source_id number and format to 9 chars
            srcid = f's{str(source_id):>09s}'

        return srcid

# #########
# Utilities
# #########
def invariant_filename(save_model_func):
    """Restore meta.filename after save_model"""

    @wraps(save_model_func)
    def save_model(model, **kwargs):
        try:
            filename = model.meta.filename
        except AttributeError:
            filename = None

        result = save_model_func(model, **kwargs)

        if filename:
            model.meta.filename = filename

        return result

    return save_model
