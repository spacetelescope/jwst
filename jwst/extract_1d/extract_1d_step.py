from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer, SourceModelContainer

from ..stpipe import Step
from . import extract
from .soss_extract import soss_extract
from .ifu import ifu_extract1d

__all__ = ["Extract1dStep"]


class Extract1dStep(Step):
    """Extract a 1-d spectrum from 2-d data

    Attributes
    ----------
    subtract_background : bool or None
        A flag which indicates whether the background should be subtracted.
        If None, the value in the extract_1d reference file will be used.
        If not None, this parameter overrides the value in the
        extract_1d reference file.

    apply_apcorr : bool
        Switch to select whether to apply an APERTURE correction during
        the Extract1dStep. Default is True.

    use_source_posn : bool or None
        If True, the source and background extraction positions specified in
        the extract1d reference file (or the default position, if there is no
        reference file) will be shifted to account for the computed position
        of the source in the data.  If None (the default), the values in the
        reference file will be used. Aperture offset is determined by computing
        the pixel location of the source based on its RA and Dec. It does not
        make sense to apply aperture offsets for extended sources, so this
        parameter can be overridden (set to False) internally by the step.

    smoothing_length : int or None
        If not None, the background regions (if any) will be smoothed
        with a boxcar function of this width along the dispersion
        direction.  This should be an odd integer.

    bkg_fit : str
        A string indicating the type of fitting to be applied to
        background values in each column (or row, if the dispersion is
        vertical). Allowed values are `poly`, `mean`, and `median`.
        Default is `None`.

    bkg_order : int or None
        If not None, a polynomial with order `bkg_order` will be fit to
        each column (or row, if the dispersion direction is vertical)
        of the background region or regions.  For a given column (row),
        one polynomial will be fit to all background regions.  The
        polynomial will be evaluated at each pixel of the source
        extraction region(s) along the column (row), and the fitted value
        will be subtracted from the data value at that pixel.
        If both `smoothing_length` and `bkg_order` are not None, the
        boxcar smoothing will be done first.

    log_increment : int
        if `log_increment` is greater than 0 (the default is 50) and the
        input data are multi-integration (which can be CubeModel or
        SlitModel), a message will be written to the log with log level
        INFO every `log_increment` integrations.  This is intended to
        provide progress information when invoking the step interactively.

    save_profile : bool
        If True, the spatial profile containing the extraction aperture
        is saved to disk.  Ignored for IFU and NIRISS SOSS extractions.

    save_scene_model : bool
        If True, a model of the 2D flux as defined by the extraction aperture
        is saved to disk.  Ignored for IFU and NIRISS SOSS extractions.

    center_xy : int or None
        A list of 2 pixel coordinate values at which to place the center
        of the IFU extraction aperture, overriding any centering done by the step.
        Two values, in x,y order, are used for extraction from IFU cubes.
        Default is None.

    ifu_autocen : bool
        Switch to turn on auto-centering for point source spectral extraction
        in IFU mode.  Default is False.

    bkg_sigma_clip : float
        Background sigma clipping value to use on background to remove outliers
        and maximize the quality of the 1d spectrum. Used for IFU mode only.

    ifu_rfcorr : bool
        Switch to select whether or not to apply a 1d residual fringe correction
        for MIRI MRS IFU spectra.  Default is False.

    ifu_set_srctype : str
        For MIRI MRS IFU data override srctype and set it to either POINT or EXTENDED.

    ifu_rscale : float
        For MRS IFU data a value for changing the extraction radius. The value provided is
        the number of PSF FWHMs to use for the extraction radius. Values accepted are between
        0.5 to 3.0. The default extraction size is set to 2 * FWHM. Values below 2 will result
        in a smaller radius, a value of 2 results in no change to the radius and a value above
        2 results in a larger extraction radius.

    ifu_covar_scale : float
        Scaling factor by which to multiply the ERR values in extracted spectra to account
        for covariance between adjacent spaxels in the IFU data cube.

    soss_atoca : bool, default=False
        Switch to toggle extraction of SOSS data with the ATOCA algorithm.
        WARNING: ATOCA results not fully validated, and require the photom step
        be turned off. Default is False, meaning SOSS data use box extraction.

    soss_threshold : float
        Threshold value above which a pixel will be included when modeling the SOSS
        trace in ATOCA. Default is 0.01.

    soss_n_os : int
        Oversampling factor of the underlying wavelength grid when modeling the SOSS
        trace in ATOCA. Default is 2.

    soss_wave_grid_in : str or SossWaveGrid or None
        Filename or SossWaveGrid containing the wavelength grid used by ATOCA
        to model each pixel valid pixel of the detector. If not given, the grid is determined
        based on an estimate of the flux (soss_estimate), the relative tolerance (soss_rtol)
        required on each pixel model and the maximum grid size (soss_max_grid_size).

    soss_wave_grid_out : str or None
        Filename to hold the wavelength grid calculated by ATOCA.

    soss_estimate : str or SpecModel or None
        Filename or SpecModel of the estimate of the target flux. The estimate must
        be a SpecModel with wavelength and flux values.

    soss_rtol : float
        The relative tolerance needed on a pixel model. It is used to determine the sampling
        of the soss_wave_grid when not directly given.

    soss_max_grid_size: int
        Maximum grid size allowed. It is used when soss_wave_grid is not provided
        to make sure the computation time or the memory used stays reasonable.

    soss_tikfac : float
        The regularization factor used for extraction in ATOCA. If left to default
        value of None, ATOCA will find an optimized value.

    soss_width : float
        Aperture width used to extract the SOSS spectrum from the decontaminated
        trace in ATOCA. Default is 40.

    soss_bad_pix : str
        Method used to handle bad pixels, accepts either "model" or "masking". Default
        method is "model".

    soss_modelname : str
        Filename for optional model output of ATOCA traces and pixel weights.

    """

    class_alias = "extract_1d"

    spec = """
    subtract_background = boolean(default=None)  # subtract background?
    apply_apcorr = boolean(default=True)  # apply aperture corrections?

    use_source_posn = boolean(default=None)  # use source coords to center extractions?
    smoothing_length = integer(default=None)  # background smoothing size
    bkg_fit = option("poly", "mean", "median", None, default=None)  # background fitting type
    bkg_order = integer(default=None, min=0)  # order of background polynomial fit
    log_increment = integer(default=50)  # increment for multi-integration log messages
    save_profile = boolean(default=False)  # save spatial profile to disk
    save_scene_model = boolean(default=False)  # save flux model to disk
    
    center_xy = float_list(min=2, max=2, default=None)  # IFU extraction x/y center
    ifu_autocen = boolean(default=False) # Auto source centering for IFU point source data.
    bkg_sigma_clip = float(default=3.0)  # background sigma clipping threshold for IFU
    ifu_rfcorr = boolean(default=False) # Apply 1d residual fringe correction
    ifu_set_srctype = option("POINT", "EXTENDED", None, default=None) # user-supplied source type
    ifu_rscale = float(default=None, min=0.5, max=3) # Radius in terms of PSF FWHM to scale extraction radii
    ifu_covar_scale = float(default=1.0) # Scaling factor to apply to errors to account for IFU cube covariance
    
    soss_atoca = boolean(default=True)  # use ATOCA algorithm
    soss_threshold = float(default=1e-2)  # TODO: threshold could be removed from inputs. Its use is too specific now.
    soss_n_os = integer(default=2)  # minimum oversampling factor of the underlying wavelength grid used when modeling trace.
    soss_wave_grid_in = input_file(default = None)  # Input wavelength grid used to model the detector
    soss_wave_grid_out = string(default = None)  # Output wavelength grid solution filename
    soss_estimate = input_file(default = None)  # Estimate used to generate the wavelength grid
    soss_rtol = float(default=1.0e-4)  # Relative tolerance needed on a pixel model
    soss_max_grid_size = integer(default=20000)  # Maximum grid size, if wave_grid not specified
    soss_tikfac = float(default=None)  # regularization factor for NIRISS SOSS extraction
    soss_width = float(default=40.)  # aperture width used to extract the 1D spectrum from the de-contaminated trace.
    soss_bad_pix = option("model", "masking", default="masking")  # method used to handle bad pixels
    soss_modelname = output_file(default = None)  # Filename for optional model output of traces and pixel weights
    """

    reference_file_types = ['extract1d', 'apcorr', 'pastasoss', 'specprofile', 'speckernel']

    def _get_extract_reference_files_by_mode(self, model, exp_type):
        """Get extraction reference files with special handling by exposure type."""
        if isinstance(model, ModelContainer):
            model = model[0]

        if exp_type not in extract.WFSS_EXPTYPES:
            extract_ref = self.get_reference_file(model, 'extract1d')
        else:
            extract_ref = 'N/A'
        if extract_ref != 'N/A':
            self.log.info(f'Using EXTRACT1D reference file {extract_ref}')

        if self.apply_apcorr:
            apcorr_ref = self.get_reference_file(model, 'apcorr')
        else:
            apcorr_ref = 'N/A'
        if apcorr_ref != 'N/A':
            self.log.info(f'Using APCORR file {apcorr_ref}')

        return extract_ref, apcorr_ref

    def _extract_soss(self, model):
        """Extract NIRISS SOSS spectra."""
        # Set the filter configuration
        if model.meta.instrument.filter == 'CLEAR':
            self.log.info('Exposure is through the GR700XD + CLEAR (science).')
            soss_filter = 'CLEAR'
        else:
            self.log.error('The SOSS extraction is implemented for the CLEAR filter only. '
                           f'Requested filter is {model.meta.instrument.filter}.')
            self.log.error('extract_1d will be skipped.')
            model.meta.cal_step.extract_1d = 'SKIPPED'
            return model

        # Set the subarray mode being processed
        if model.meta.subarray.name == 'SUBSTRIP256':
            self.log.info('Exposure is in the SUBSTRIP256 subarray.')
            self.log.info('Traces 1 and 2 will be modelled and decontaminated before extraction.')
            subarray = 'SUBSTRIP256'
        elif model.meta.subarray.name == 'SUBSTRIP96':
            self.log.info('Exposure is in the SUBSTRIP96 subarray.')
            self.log.info('Traces of orders 1 and 2 will be modelled but only order 1 '
                          'will be decontaminated before extraction.')
            subarray = 'SUBSTRIP96'
        else:
            self.log.error('The SOSS extraction is implemented for the SUBSTRIP256 '
                           'and SUBSTRIP96 subarrays only. Subarray is currently '
                           f'{model.meta.subarray.name}.')
            self.log.error('Extract1dStep will be skipped.')
            model.meta.cal_step.extract_1d = 'SKIPPED'
            return model

        # Load reference files.
        pastasoss_ref_name = self.get_reference_file(model, 'pastasoss')
        specprofile_ref_name = self.get_reference_file(model, 'specprofile')
        speckernel_ref_name = self.get_reference_file(model, 'speckernel')

        # Build SOSS kwargs dictionary.
        soss_kwargs = dict()
        soss_kwargs['threshold'] = self.soss_threshold
        soss_kwargs['n_os'] = self.soss_n_os
        soss_kwargs['tikfac'] = self.soss_tikfac
        soss_kwargs['width'] = self.soss_width
        soss_kwargs['bad_pix'] = self.soss_bad_pix
        soss_kwargs['subtract_background'] = self.subtract_background
        soss_kwargs['rtol'] = self.soss_rtol
        soss_kwargs['max_grid_size'] = self.soss_max_grid_size
        soss_kwargs['wave_grid_in'] = self.soss_wave_grid_in
        soss_kwargs['wave_grid_out'] = self.soss_wave_grid_out
        soss_kwargs['estimate'] = self.soss_estimate
        soss_kwargs['atoca'] = self.soss_atoca
        # Set flag to output the model and the tikhonov tests
        soss_kwargs['model'] = True if self.soss_modelname else False

        # Run the extraction.
        result, ref_outputs, atoca_outputs = soss_extract.run_extract1d(
            model,
            pastasoss_ref_name,
            specprofile_ref_name,
            speckernel_ref_name,
            subarray,
            soss_filter,
            soss_kwargs)

        # Set the step flag to complete
        if result is None:
            return None
        else:
            result.meta.cal_step.extract_1d = 'COMPLETE'
            result.meta.target.source_type = None

            model.close()

            if self.soss_modelname:
                soss_modelname = self.make_output_path(
                    basepath=self.soss_modelname,
                    suffix='SossExtractModel'
                )
                ref_outputs.save(soss_modelname)

        if self.soss_modelname:
            soss_modelname = self.make_output_path(
                basepath=self.soss_modelname,
                suffix='AtocaSpectra'
            )
            atoca_outputs.save(soss_modelname)

        return result

    def _extract_ifu(self, model, exp_type, extract_ref, apcorr_ref):
        """Extract IFU spectra from a single datamodel."""
        source_type = model.meta.target.source_type
        if self.ifu_set_srctype is not None and exp_type == 'MIR_MRS':
            source_type = self.ifu_set_srctype
            self.log.info(f"Overriding source type and setting it to {self.ifu_set_srctype}")

        result = ifu_extract1d(
            model, extract_ref, source_type, self.subtract_background,
            self.bkg_sigma_clip, apcorr_ref, self.center_xy,
            self.ifu_autocen, self.ifu_rfcorr, self.ifu_rscale,
            self.ifu_covar_scale
        )
        return result

    def _save_intermediate(self, intermediate_model, suffix):
        """Save an intermediate output file."""
        if isinstance(intermediate_model, ModelContainer):
            # Save the profile with the slit name + suffix 'profile'
            for model in intermediate_model:
                slit = str(model.name).lower()
                output_path = self.make_output_path(suffix=f'{slit}_{suffix}')
                self.log.info(f"Saving {suffix} {output_path}")
                model.save(output_path)
        else:
            # Only one profile - just use the suffix 'profile'
            output_path = self.make_output_path(suffix=suffix)
            self.log.info(f"Saving {suffix} {output_path}")
            intermediate_model.save(output_path)
        intermediate_model.close()

    def process(self, input):
        """Execute the step.

        Parameters
        ----------
        input: JWST data model

        Returns
        -------
        JWST data model
            This will be `input_model` if the step was skipped; otherwise,
            it will be a model containing 1-D extracted spectra.
        """

        # Open the input and figure out what type of model it is
        if isinstance(input, ModelContainer):
            input_model = input
        else:
            input_model = datamodels.open(input)

        if isinstance(input_model, (datamodels.CubeModel, datamodels.ImageModel,
                                    datamodels.SlitModel, datamodels.IFUCubeModel,
                                    ModelContainer, SourceModelContainer)):
            # Acceptable input type, just log it
            self.log.debug(f'Input is a {str(type(input_model))}.')
        elif isinstance(input_model, datamodels.MultiSlitModel):
            # If input is multislit, with 3D calints, skip the step
            self.log.debug('Input is a MultiSlitModel')
            if len((input_model[0]).shape) == 3:
                self.log.warning('3D input is unsupported; step will be skipped')
                input_model.meta.cal_step.extract_1d = 'SKIPPED'
                return input_model
        else:
            self.log.error(f'Input is a {str(type(input_model))}, ')
            self.log.error('which was not expected for extract_1d.')
            self.log.error('The extract_1d step will be skipped.')
            input_model.meta.cal_step.extract_1d = 'SKIPPED'
            return input_model

        if not isinstance(input_model, ModelContainer):
            exp_type = input_model.meta.exposure.type

            # Make the input iterable
            input_model = [input_model]
        else:
            exp_type = input_model[0].meta.exposure.type
        self.log.debug(f"Input for EXP_TYPE {exp_type} contains {len(input_model)} items")

        if len(input_model) > 1 and exp_type in extract.WFSS_EXPTYPES:
            # For WFSS level-3, the input is a single entry of a
            # SourceContainer, which contains a list of multiple
            # SlitModels for a single source. Send the whole container
            # into extract1d and put all results in a single product.
            input_model = [input_model]

        if exp_type == 'NIS_SOSS':
            # Data is NIRISS SOSS observation, use its own extraction routines
            self.log.info(
                'Input is a NIRISS SOSS observation, the specialized SOSS '
                'extraction (ATOCA) will be used.')

            # There is only one input model for this mode
            model = input_model[0]
            result = self._extract_soss(model)

        else:
            result = ModelContainer()
            for model in input_model:
                # Get the reference file names
                extract_ref, apcorr_ref = self._get_extract_reference_files_by_mode(
                    model, exp_type)

                profile = None
                scene_model = None
                if isinstance(model, datamodels.IFUCubeModel):
                    # Call the IFU specific extraction routine
                    extracted = self._extract_ifu(model, exp_type, extract_ref, apcorr_ref)
                else:
                    # Call the general extraction routine
                    extracted, profile, scene_model = extract.run_extract1d(
                        model,
                        extract_ref,
                        apcorr_ref,
                        self.smoothing_length,
                        self.bkg_fit,
                        self.bkg_order,
                        self.log_increment,
                        self.subtract_background,
                        self.use_source_posn,
                        self.save_profile,
                        self.save_scene_model
                    )

                # Set the step flag to complete in each model
                extracted.meta.cal_step.extract_1d = 'COMPLETE'
                result.append(extracted)
                del extracted

                # Save profile if needed
                if self.save_profile and profile is not None:
                    self._save_intermediate(profile, 'profile')

                # Save model if needed
                if self.save_scene_model and scene_model is not None:
                    self._save_intermediate(scene_model, 'scene_model')

            # If only one result, return the model instead of the container
            if len(result) == 1:
                result = result[0]

        return result
