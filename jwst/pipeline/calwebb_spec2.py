import os
from collections import defaultdict
import os.path as op
import traceback
import numpy as np

from stdatamodels.jwst import datamodels

from ..assign_wcs.util import NoDataOnDetectorError
from ..lib.exposure_types import is_nrs_ifu_flatlamp, is_nrs_ifu_linelamp, is_nrs_slit_linelamp
from ..stpipe import Pipeline

# step imports
from ..assign_wcs import assign_wcs_step
from ..background import background_step
from ..barshadow import barshadow_step
from ..cube_build import cube_build_step
from ..extract_1d import extract_1d_step
from ..extract_2d import extract_2d_step
from ..flatfield import flat_field_step
from ..fringe import fringe_step
from ..residual_fringe import residual_fringe_step
from ..imprint import imprint_step
from ..master_background import master_background_mos_step
from ..msaflagopen import msaflagopen_step
from ..pathloss import pathloss_step
from ..photom import photom_step
from ..pixel_replace import pixel_replace_step
from ..resample import resample_spec_step
from ..srctype import srctype_step
from ..straylight import straylight_step
from ..wavecorr import wavecorr_step
from ..wfss_contam import wfss_contam_step

__all__ = ['Spec2Pipeline']

# Classify various exposure types.
NRS_SLIT_TYPES = ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ', 'NRS_MSASPEC',
                  'NRS_LAMP', 'NRS_AUTOWAVE', 'NRS_AUTOFLAT']
WFSS_TYPES = ["NIS_WFSS", "NRC_GRISM", "NRC_WFSS"]
GRISM_TYPES = ['NRC_TSGRISM'] + WFSS_TYPES


class Spec2Pipeline(Pipeline):
    """
    Spec2Pipeline: Processes JWST spectroscopic exposures from Level 2a to 2b.
    Accepts a single exposure or an association as input.

    Included steps are:
    assign_wcs, background subtraction, NIRSpec MSA imprint subtraction,
    NIRSpec MSA bad shutter flagging, 2-D subwindow extraction, flat field,
    source type decision, straylight, fringe, residual_fringe, pathloss,
    barshadow,  photom, resample_spec, cube_build, and extract_1d.
    """

    class_alias = "calwebb_spec2"

    spec = """
        save_bsub = boolean(default=False)        # Save background-subtracted science
        fail_on_exception = boolean(default=True) # Fail if any product fails.
        save_wfss_esec = boolean(default=False)   # Save WFSS e-/sec image
    """

    # Define aliases to steps
    step_defs = {
        'bkg_subtract': background_step.BackgroundStep,
        'assign_wcs': assign_wcs_step.AssignWcsStep,
        'imprint_subtract': imprint_step.ImprintStep,
        'msa_flagging': msaflagopen_step.MSAFlagOpenStep,
        'extract_2d': extract_2d_step.Extract2dStep,
        'master_background_mos': master_background_mos_step.MasterBackgroundMosStep,
        'wavecorr': wavecorr_step.WavecorrStep,
        'flat_field': flat_field_step.FlatFieldStep,
        'srctype': srctype_step.SourceTypeStep,
        'straylight': straylight_step.StraylightStep,
        'fringe': fringe_step.FringeStep,
        'residual_fringe': residual_fringe_step.ResidualFringeStep,
        'pathloss': pathloss_step.PathLossStep,
        'barshadow': barshadow_step.BarShadowStep,
        'wfss_contam': wfss_contam_step.WfssContamStep,
        'photom': photom_step.PhotomStep,
        'pixel_replace': pixel_replace_step.PixelReplaceStep,
        'resample_spec': resample_spec_step.ResampleSpecStep,
        'cube_build': cube_build_step.CubeBuildStep,
        'extract_1d': extract_1d_step.Extract1dStep
    }

    # Main processing
    def process(self, data):
        """Entrypoint for this pipeline

        Parameters
        ----------
        input: str, Level2 Association, or ~jwst.datamodels.JwstDataModel
            The exposure or association of exposures to process
        """
        self.log.info('Starting calwebb_spec2 ...')

        # Setup step parameters required by the pipeline.
        self.resample_spec.save_results = self.save_results
        self.resample_spec.suffix = 's2d'
        self.cube_build.output_type = 'multi'
        self.cube_build.save_results = False
        self.cube_build.skip_dqflagging = True
        self.extract_1d.save_results = self.save_results

        # Retrieve the input(s)
        asn = self.load_as_level2_asn(data)

        # Each exposure is a product in the association.
        # Process each exposure.  Delay reporting failures until the end.
        results = []
        failures = []
        for product in asn['products']:
            self.log.info('Processing product {}'.format(product['name']))
            self.output_file = product['name']
            try:
                getattr(asn, 'filename')
            except AttributeError:
                asn.filename = "singleton"
            try:
                result = self.process_exposure_product(
                    product,
                    asn['asn_pool'],
                    asn.filename
                )
            except NoDataOnDetectorError as exception:
                # This error merits a special return
                # status if run from the command line.
                # Bump it up now.
                raise exception
            except Exception:
                traceback.print_exc()
                failures.append(traceback.format_exc())
            else:
                if result is not None:
                    results.append(result)

        if len(failures) > 0 and self.fail_on_exception:
            raise RuntimeError('\n'.join(failures))

        # We're done
        self.log.info('Ending calwebb_spec2')

        self.output_use_model = True
        self.suffix = False
        return results

    # Process each exposure
    def process_exposure_product(
            self,
            exp_product,
            pool_name=' ',
            asn_file=' '
    ):
        """Process an exposure found in the association product

        Parameters
        ---------
        exp_product: dict
            A Level2b association product.
        """

        # Find all the member types in the product
        members_by_type = defaultdict(list)
        for member in exp_product['members']:
            members_by_type[member['exptype'].lower()].append(member['expname'])

        # Get the science member. Technically there should only be
        # one. We'll just get the first one found.
        science_member = members_by_type['science']
        if len(science_member) != 1:
            self.log.warning(
                'Wrong number of science exposures found in {}'.format(
                    exp_product['name']
                )
            )
            self.log.warning('    Using only first one.')
        science_member = science_member[0]

        self.log.info('Working on input %s ...', science_member)
        with self.open_model(science_member) as science:
            exp_type = science.meta.exposure.type
            if isinstance(science, datamodels.CubeModel):
                multi_int = True
            else:
                multi_int = False

            # Suffixes are dependent on whether the science is multi-integration or not.
            if multi_int:
                suffix = 'calints'
                self.extract_1d.suffix = 'x1dints'
            else:
                suffix = 'cal'
                self.extract_1d.suffix = 'x1d'

            # Apply WCS info
            # check the datamodel to see if it's
            # a grism image, if so get the catalog
            # name from the asn and record it to the meta
            if exp_type in WFSS_TYPES:
                try:
                    science.meta.source_catalog = os.path.basename(members_by_type['sourcecat'][0])
                    self.log.info('Using sourcecat file {}'.format(science.meta.source_catalog))
                    science.meta.segmentation_map = os.path.basename(members_by_type['segmap'][0])
                    self.log.info('Using segmentation map {}'.format(science.meta.segmentation_map))
                    science.meta.direct_image = os.path.basename(members_by_type['direct_image'][0])
                    self.log.info('Using direct image {}'.format(science.meta.direct_image))
                except IndexError:
                    if science.meta.source_catalog is None:
                        raise IndexError("No source catalog specified in association or datamodel")

            # Decide on what steps can actually be accomplished based on the
            # provided input.

            self._step_verification(exp_type, science, members_by_type, multi_int)

            # Start processing the individual steps.
            # `assign_wcs` is the critical step. Without it, processing
            # cannot proceed.
            assign_wcs_exception = None
            try:
                calibrated = self.assign_wcs(science)
            except Exception as exception:
                assign_wcs_exception = exception
            if assign_wcs_exception is not None or \
               calibrated.meta.cal_step.assign_wcs != 'COMPLETE':
                message = (
                    'Assign_wcs processing was skipped.'
                    '\nAborting remaining processing for this exposure.'
                    '\nNo output product will be created.'
                )
                if self.assign_wcs.skip:
                    self.log.warning(message)
                    return
                else:
                    self.log.error(message)
                    if assign_wcs_exception is not None:
                        raise assign_wcs_exception
                    else:
                        raise RuntimeError('Cannot determine WCS.')

        # Steps whose order is the same for all types of input.

        # Leakcal subtraction (imprint)  occurs before background subtraction on a per-exposure basis.
        # If there is only one `imprint` member, this imprint exposure is subtracted from all the
        # science and background exposures.  Otherwise, there will be as many `imprint` members as
        # there are science plus background members.

        calibrated = self.imprint_subtract(calibrated, members_by_type['imprint'])

        # for each background image subtract an associated leak cal
        for i, bkg_file in enumerate(members_by_type['background']):
            bkg_imprint_sub = self.imprint_subtract(bkg_file, members_by_type['imprint'])
            members_by_type['background'][i] = bkg_imprint_sub

        calibrated = self.bkg_subtract(calibrated, members_by_type['background'])

        calibrated = self.msa_flagging(calibrated)

        # The order of the next few steps is tricky, depending on mode:
        # WFSS/Grism data need flat_field before extract_2d, but other modes
        # need extract_2d first. Furthermore, NIRSpec MOS and FS need
        # srctype and wavecorr before flat_field.
        if exp_type in GRISM_TYPES:
            calibrated = self._process_grism(calibrated)
        elif exp_type in NRS_SLIT_TYPES:
            calibrated = self._process_nirspec_slits(calibrated)
        elif exp_type == 'NIS_SOSS':
            calibrated = self._process_niriss_soss(calibrated)
        else:
            calibrated = self._process_common(calibrated)

        # Record ASN pool and table names in output
        calibrated.meta.asn.pool_name = pool_name
        calibrated.meta.asn.table_name = op.basename(asn_file)
        calibrated.meta.filename = self.make_output_path(suffix=suffix)

        # Produce a resampled product, either via resample_spec for
        # "regular" spectra or cube_build for IFU data. No resampled
        # product is produced for time-series modes.
        if exp_type in ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'MIR_LRS-FIXEDSLIT'] \
           and not isinstance(calibrated, datamodels.CubeModel):

            # Call the resample_spec step for 2D slit data
            resampled = calibrated.copy()
            resampled = self.resample_spec(resampled)

        elif is_nrs_slit_linelamp(calibrated):

            # Call resample_spec for NRS 2D line lamp slit data
            resampled = calibrated.copy()
            resampled = self.resample_spec(resampled)

        elif (exp_type in ['MIR_MRS', 'NRS_IFU']) or is_nrs_ifu_linelamp(calibrated):

            # Call the cube_build step for IFU data;
            # always create a single cube containing multiple
            # wavelength bands

            resampled = calibrated.copy()
            resampled = self.cube_build(resampled)
            if not self.cube_build.skip:
                self.save_model(resampled[0], 's3d')
        else:
            resampled = calibrated

        # Extract a 1D spectrum from the 2D/3D data
        if exp_type in ['MIR_MRS', 'NRS_IFU'] and self.cube_build.skip:
            # Skip extract_1d for IFU modes where no cube was built
            self.extract_1d.skip = True

        # SOSS data need to run photom on x1d products and optionally save the photom
        # output, while all other exptypes simply run extract_1d.
        if exp_type == 'NIS_SOSS':
            if multi_int:
                self.photom.suffix = 'x1dints'
            else:
                self.photom.suffix = 'x1d'
            self.extract_1d.save_results = False
            x1d = resampled.copy()
            x1d = self.extract_1d(x1d)

            # Possible that no fit was possible - if so, skip photom
            if x1d is not None:
                self.photom.save_results = self.save_results
                x1d = self.photom(x1d)
            else:
                self.log.warning("Extract_1d did not return a DataModel - skipping photom.")
        else:
            x1d = resampled.copy()
            x1d = self.extract_1d(x1d)

        resampled.close()
        if x1d is not None:
            x1d.close()

        # That's all folks
        self.log.info('Finished processing product {}'.format(exp_product['name']))

        return calibrated

    def _step_verification(self, exp_type, input, members_by_type, multi_int):
        """Verify whether requested steps can operate on the given data

        Though ideally this would all be controlled through the pipeline
        parameters, the desire to keep the number of config files down has
        pushed the logic into code.

        Once step and pipeline parameters are retrieved from CRDS, this
        logic can be removed.
        """

        # Check for image-to-image background subtraction can be done.
        if not self.bkg_subtract.skip:
            if exp_type in WFSS_TYPES or len(members_by_type['background']) > 0:

                if exp_type in WFSS_TYPES:
                    members_by_type['background'] = []           # will be overwritten by the step

                # Setup for saving
                self.bkg_subtract.suffix = 'bsub'
                if multi_int:
                    self.bkg_subtract.suffix = 'bsubints'

                # Backwards compatibility
                if self.save_bsub:
                    self.bkg_subtract.save_results = True
            else:
                self.log.debug('Science data does not allow direct background subtraction. Skipping "bkg_subtract".')
                self.bkg_subtract.skip = True

        # Check for imprint subtraction. If we have a background then we could have an imprint image associated with the
        # background.
        imprint = members_by_type['imprint']
        if not self.imprint_subtract.skip:
            if len(imprint) > 0 and (exp_type in ['NRS_MSASPEC', 'NRS_IFU']
                                     or is_nrs_ifu_flatlamp(input)):
                if len(imprint) > 1 and (exp_type in ['NRS_MSASPEC'] or is_nrs_ifu_flatlamp(input)):
                    self.log.warning('Wrong number of imprint members')
                    members_by_type['imprint'] = imprint[0]
            else:
                self.log.debug('Science data does not allow imprint processing. Skipping "imprint_subtraction".')
                self.imprint_subtract.skip = True

        # Check for NIRSpec MSA bad shutter flagging.
        if not self.msa_flagging.skip and exp_type not in ['NRS_MSASPEC', 'NRS_IFU', 'NRS_LAMP',
                                                           'NRS_AUTOFLAT', 'NRS_AUTOWAVE']:
            self.log.debug('Science data does not allow MSA flagging. Skipping "msa_flagging".')
            self.msa_flagging.skip = True

        # Check for straylight correction for MIRI MRS.
        if not self.straylight.skip and exp_type != 'MIR_MRS':
            self.log.debug('Science data does not allow stray light correction. Skipping "straylight".')
            self.straylight.skip = True

        # Check for residual_fringe correction for MIRI MRS.
        if not self.residual_fringe.skip and exp_type != 'MIR_MRS':
            self.log.debug('Science data does not allow residual fringe correction. Skipping "residual fringe".')
            self.residual_fringe.skip = True

        # Apply the fringe correction for MIRI MRS
        if not self.fringe.skip and exp_type != 'MIR_MRS':
            self.log.debug('Science data does not allow fringe correction. Skipping "fringe".')
            self.fringe.skip = True

        # Apply pathloss correction to MIRI LRS, NIRSpec, and NIRISS SOSS exposures
        if not self.pathloss.skip and exp_type not in ['MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT',
                                                       'NRS_MSASPEC', 'NRS_IFU', 'NIS_SOSS']:
            self.log.debug('Science data does not allow pathloss correction. Skipping "pathloss".')
            self.pathloss.skip = True

        # Apply barshadow correction to NIRSPEC MSA exposures
        if not self.barshadow.skip and exp_type != 'NRS_MSASPEC':
            self.log.debug('Science data does not allow barshadow correction. Skipping "barshadow".')
            self.barshadow.skip = True

        # Apply master background only to NIRSPEC MSA exposures
        if not self.master_background_mos.skip and exp_type != 'NRS_MSASPEC':
            self.log.debug('Science data does not allow master background correction. Skipping "master_background_mos".')
            self.master_background_mos.skip = True

        # Apply WFSS contamination correction only to WFSS exposures
        if not self.wfss_contam.skip and exp_type not in WFSS_TYPES:
            self.log.debug('Science data does not allow WFSS contamination correction. Skipping "wfss_contam".')
            self.wfss_contam.skip = True

    def _process_grism(self, data):
        """WFSS & Grism processing

        WFSS/Grism data need flat_field before extract_2d.
        """

        # Apply flat-field correction
        calibrated = self.flat_field(data)

        # Create and save a WFSS e-/sec image, if requested
        if self.save_wfss_esec:
            self.log.info('Creating WFSS e-/sec product')

            # Find and load the gain reference file that we need
            gain_filename = self.get_reference_file(calibrated, 'gain')
            self.log.info('Using GAIN reference file %s', gain_filename)
            with datamodels.GainModel(gain_filename) as gain_model:

                # Always use the full-frame version of the gain ref file,
                # even the science data are taken with a subarray
                gain_image = gain_model.data

                # Compute the simple mean of the gain image, excluding reference pixels.
                # The gain ref file doesn't have a DQ array that can be used to
                # mask bad values, so manually exclude NaN's and gain <= 0.
                gain_image[gain_image <= 0.] = np.nan
                mean_gain = np.nanmean(gain_image[4:-4, 4:-4])
                self.log.info('mean gain = %s', mean_gain)

                # Apply gain to the intermediate WFSS image
                wfss_esec = calibrated.copy()
                mean_gain_sqr = mean_gain ** 2
                wfss_esec.data *= mean_gain
                wfss_esec.var_poisson *= mean_gain_sqr
                wfss_esec.var_rnoise *= mean_gain_sqr
                wfss_esec.var_flat *= mean_gain_sqr
                wfss_esec.err = np.sqrt(wfss_esec.var_poisson + wfss_esec.var_rnoise + wfss_esec.var_flat)

                # Save the WFSS e-/sec image
                self.save_model(wfss_esec, suffix='esec', force=True)
                del wfss_esec

        # Continue with remaining calibration steps, using the original
        # DN/sec image
        calibrated = self.extract_2d(calibrated)
        calibrated = self.srctype(calibrated)
        calibrated = self.straylight(calibrated)
        calibrated = self.fringe(calibrated)
        calibrated = self.pathloss(calibrated)
        calibrated = self.barshadow(calibrated)
        calibrated = self.wfss_contam(calibrated)
        calibrated = self.photom(calibrated)
        calibrated = self.pixel_replace(calibrated)

        return calibrated

    def _process_nirspec_slits(self, data):
        """Process NIRSpec

        NIRSpec MOS and FS need srctype and wavecorr before flat_field.
        Also have to deal with master background operations.
        """
        calibrated = self.extract_2d(data)
        calibrated = self.srctype(calibrated)
        calibrated = self.master_background_mos(calibrated)
        calibrated = self.wavecorr(calibrated)
        calibrated = self.flat_field(calibrated)
        calibrated = self.pathloss(calibrated)
        calibrated = self.barshadow(calibrated)
        calibrated = self.photom(calibrated)
        calibrated = self.pixel_replace(calibrated)

        return calibrated

    def _process_niriss_soss(self, data):
        """Process SOSS

        New SOSS extraction requires input to extract_1d step in units
        of DN/s, with photom step to be run afterwards.
        """
        calibrated = self.srctype(data)
        calibrated = self.flat_field(calibrated)
        calibrated = self.straylight(calibrated)
        calibrated = self.fringe(calibrated)
        calibrated = self.pathloss(calibrated)
        calibrated = self.barshadow(calibrated)
        calibrated = self.pixel_replace(calibrated)

        return calibrated

    def _process_common(self, data):
        """Common spectral processing"""
        calibrated = self.srctype(data)
        calibrated = self.straylight(calibrated)
        calibrated = self.flat_field(calibrated)
        calibrated = self.fringe(calibrated)
        calibrated = self.pathloss(calibrated)
        calibrated = self.barshadow(calibrated)
        calibrated = self.photom(calibrated)
        calibrated = self.residual_fringe(calibrated)  # only run on MIRI_MRS data
        calibrated = self.pixel_replace(calibrated)

        return calibrated
