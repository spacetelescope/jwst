from ..stpipe import Step
from .. import datamodels
from . import extract
from .soss_extract import soss_extract

__all__ = ["Extract1dStep"]


class Extract1dStep(Step):
    """Extract a 1-d spectrum from 2-d data

    Attributes
    ----------
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

    bkg_sigma_clip : float
        Background sigma clipping value to use on background to remove outliers
        and maximize the quality of the 1d spectrum

    log_increment : int
        if `log_increment` is greater than 0 (the default is 50) and the
        input data are multi-integration (which can be CubeModel or
        SlitModel), a message will be written to the log with log level
        INFO every `log_increment` integrations.  This is intended to
        provide progress information when invoking the step interactively.

    subtract_background : bool or None
        A flag which indicates whether the background should be subtracted.
        If None, the value in the extract_1d reference file will be used.
        If not None, this parameter overrides the value in the
        extract_1d reference file.

    use_source_posn : bool or None
        If True, the source and background extraction positions specified in
        the extract1d reference file (or the default position, if there is no
        reference file) will be shifted to account for the computed position
        of the source in the data.  If None (the default), the values in the
        reference file will be used. Aperture offset is determined by computing
        the pixel location of the source based on its RA and Dec. It does not
        make sense to apply aperture offsets for extended sources, so this
        parameter can be overridden (set to False) internally by the step.

    center_xy : int or None
        A list of 2 pixel coordinate values at which to place the center
        of the IFU extraction aperture, overriding any centering done by the step.
        Two values, in x,y order, are used for extraction from IFU cubes.
        Default is None.

    apply_apcorr : bool
        Switch to select whether or not to apply an APERTURE correction during
        the Extract1dStep. Default is True

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

    soss_transform : list[float]
        Rotation applied to the reference files to match the observation orientation.
        Default is None.

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

    soss_estimate : str or SpecModel or None
        Filename or SpecModel of the estimate of the target flux. The estimate must
        be a SpecModel with wavelength and flux values.

    soss_wave_grid_in : str or SossWaveGrid or None
        Filename or SossWaveGrid containing the wavelength grid used by ATOCA
        to model each pixel valid pixel of the detector. If not given, the grid is determined
        based on an estimate of the flux (soss_estimate), the relative tolerance (soss_rtol)
        required on each pixel model and the maximum grid size (soss_max_grid_size).

    soss_wave_grid_out : str or None
        Filename to hold the wavelength grid calculated by ATOCA.

    soss_rtol : float
        The relative tolerance needed on a pixel model. It is used to determine the sampling
        of the soss_wave_grid when not directly given.

    soss_max_grid_size: int
        Maximum grid size allowed. It is used when soss_wave_grid is not provided
        to make sure the computation time or the memory used stays reasonable.
    """

    class_alias = "extract_1d"

    spec = """
    smoothing_length = integer(default=None)  # background smoothing size
    bkg_fit = option("poly", "mean", "median", None, default=None)  # background fitting type
    bkg_order = integer(default=None, min=0)  # order of background polynomial fit
    bkg_sigma_clip = float(default=3.0)  # background sigma clipping threshold
    log_increment = integer(default=50)  # increment for multi-integration log messages
    subtract_background = boolean(default=None)  # subtract background?
    use_source_posn = boolean(default=None)  # use source coords to center extractions?
    center_xy = int_list(min=2, max=2, default=None)  # IFU extraction x/y center
    apply_apcorr = boolean(default=True)  # apply aperture corrections?
    soss_atoca = boolean(default=True)  # use ATOCA algorithm
    soss_threshold = float(default=1e-2)  # TODO: threshold could be removed from inputs. Its use is too specific now.
    soss_n_os = integer(default=2)  # minimum oversampling factor of the underlying wavelength grid used when modeling trace.
    soss_wave_grid_in = input_file(default = None)  # Input wavelength grid used to model the detector
    soss_wave_grid_out = string(default = None)  # Output wavelength grid solution filename
    soss_estimate = input_file(default = None)  # Estimate used to generate the wavelength grid
    soss_rtol = float(default=1.0e-4)  # Relative tolerance needed on a pixel model
    soss_max_grid_size = integer(default=20000)  # Maximum grid size, if wave_grid not specified
    soss_transform = list(default=None, min=3, max=3)  # rotation applied to the ref files to match observation.
    soss_tikfac = float(default=None)  # regularization factor for NIRISS SOSS extraction
    soss_width = float(default=40.)  # aperture width used to extract the 1D spectrum from the de-contaminated trace.
    soss_bad_pix = option("model", "masking", default="masking")  # method used to handle bad pixels
    soss_modelname = output_file(default = None)  # Filename for optional model output of traces and pixel weights
    """

    reference_file_types = ['extract1d', 'apcorr', 'wavemap', 'spectrace', 'specprofile', 'speckernel']

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
        input_model = datamodels.open(input)

        was_source_model = False  # default value
        if isinstance(input_model, datamodels.CubeModel):
            # It's a 3-D multi-integration model
            self.log.debug('Input is a CubeModel for a multiple integ. file')
        elif isinstance(input_model, datamodels.ImageModel):
            # It's a single 2-D image. This could be a resampled 2-D image
            self.log.debug('Input is an ImageModel')
        elif isinstance(input_model, datamodels.SourceModelContainer):
            self.log.debug('Input is a SourceModelContainer')
            was_source_model = True
        elif isinstance(input_model, datamodels.ModelContainer):
            self.log.debug('Input is a ModelContainer')
        elif isinstance(input_model, datamodels.MultiSlitModel):
            self.log.debug('Input is a MultiSlitModel')
        elif isinstance(input_model, datamodels.MultiExposureModel):
            self.log.warning('Input is a MultiExposureModel, '
                             'which is not currently supported')
        elif isinstance(input_model, datamodels.IFUCubeModel):
            self.log.debug('Input is an IFUCubeModel')
        elif isinstance(input_model, datamodels.SlitModel):
            # NRS_BRIGHTOBJ and MIRI LRS fixed-slit (resampled) modes
            self.log.debug('Input is a SlitModel')
        else:
            self.log.error(f'Input is a {str(type(input_model))}, ')
            self.log.error('which was not expected for extract_1d')
            self.log.error('extract_1d will be skipped.')
            input_model.meta.cal_step.extract_1d = 'SKIPPED'
            return input_model

        # ______________________________________________________________________
        # Do the extraction for ModelContainer - this might only be WFSS data
        if isinstance(input_model, datamodels.ModelContainer):

            # This is the branch  WFSS data take
            if len(input_model) > 1:
                self.log.debug(f"Input contains {len(input_model)} items")

                # --------------------------------------------------------------
                # Data is WFSS
                if input_model[0].meta.exposure.type in extract.WFSS_EXPTYPES:

                    # For WFSS level-3, the input is a single entry of a
                    # SourceContainer, which contains a list of multiple
                    # SlitModels for a single source. Send the whole list
                    # into extract1d and put all results in a single product.
                    apcorr_ref = (
                        self.get_reference_file(input_model[0], 'apcorr') if self.apply_apcorr is True else 'N/A'
                    )

                    if apcorr_ref == 'N/A':
                        self.log.info('APCORR reference file name is "N/A"')
                        self.log.info('APCORR will NOT be applied')
                    else:
                        self.log.info(f'Using APCORR file {apcorr_ref}')

                    extract_ref = 'N/A'
                    self.log.info('No EXTRACT1D reference file will be used')

                    result = extract.run_extract1d(
                        input_model,
                        extract_ref,
                        apcorr_ref,
                        self.smoothing_length,
                        self.bkg_fit,
                        self.bkg_order,
                        self.bkg_sigma_clip,
                        self.log_increment,
                        self.subtract_background,
                        self.use_source_posn,
                        self.center_xy,
                        was_source_model=was_source_model
                    )
                    # Set the step flag to complete
                    result.meta.cal_step.extract_1d = 'COMPLETE'

                # --------------------------------------------------------------
                # Data is a ModelContainer but is not WFSS
                else:
                    result = datamodels.ModelContainer()
                    for model in input_model:
                        # Get the reference file names
                        extract_ref = self.get_reference_file(model, 'extract1d')
                        self.log.info(f'Using EXTRACT1D reference file {extract_ref}')

                        apcorr_ref = self.get_reference_file(model, 'apcorr') if self.apply_apcorr is True else 'N/A'

                        if apcorr_ref == 'N/A':
                            self.log.info('APCORR reference file name is "N/A"')
                            self.log.info('APCORR will NOT be applied')
                        else:
                            self.log.info(f'Using APCORR file {apcorr_ref}')

                        temp = extract.run_extract1d(
                            model,
                            extract_ref,
                            apcorr_ref,
                            self.smoothing_length,
                            self.bkg_fit,
                            self.bkg_order,
                            self.bkg_sigma_clip,
                            self.log_increment,
                            self.subtract_background,
                            self.use_source_posn,
                            self.center_xy,
                            was_source_model=was_source_model,
                        )
                        # Set the step flag to complete in each MultiSpecModel
                        temp.meta.cal_step.extract_1d = 'COMPLETE'
                        result.append(temp)
                        del temp
            # ------------------------------------------------------------------------
            # Still in ModelContainer type, but only 1 model
            elif len(input_model) == 1:
                if input_model[0].meta.exposure.type in extract.WFSS_EXPTYPES:
                    extract_ref = 'N/A'
                    self.log.info('No EXTRACT1D reference file will be used')
                else:
                    # Get the extract1d reference file name for the one model in input
                    extract_ref = self.get_reference_file(input_model[0], 'extract1d')
                    self.log.info(f'Using EXTRACT1D reference file {extract_ref}')

                apcorr_ref = self.get_reference_file(input_model[0], 'apcorr') if self.apply_apcorr is True else 'N/A'

                if apcorr_ref == 'N/A':
                    self.log.info('APCORR reference file name is "N/A"')
                    self.log.info('APCORR will NOT be applied')
                else:
                    self.log.info(f'Using APCORR file {apcorr_ref}')

                result = extract.run_extract1d(
                    input_model[0],
                    extract_ref,
                    apcorr_ref,
                    self.smoothing_length,
                    self.bkg_fit,
                    self.bkg_order,
                    self.bkg_sigma_clip,
                    self.log_increment,
                    self.subtract_background,
                    self.use_source_posn,
                    self.center_xy,
                    was_source_model=was_source_model,
                )

                # Set the step flag to complete
                result.meta.cal_step.extract_1d = 'COMPLETE'
            else:
                self.log.error('Input model is empty;')
                self.log.error('extract_1d will be skipped.')
                return input_model

        # ______________________________________________________________________
        # Data that is not a ModelContainer (IFUCube and other single models)
        else:
            # Data is NRISS SOSS observation.
            if input_model.meta.exposure.type == 'NIS_SOSS':

                self.log.info(
                    'Input is a NIRISS SOSS observation, the specialized SOSS extraction (ATOCA) will be used.')

                # Set the filter configuration
                if input_model.meta.instrument.filter == 'CLEAR':
                    self.log.info('Exposure is through the GR700XD + CLEAR (science).')
                    soss_filter = 'CLEAR'
                elif input_model.meta.instrument.filter == 'F277W':
                    self.log.info('Exposure is through the GR700XD + F277W (calibration).')
                    soss_filter = 'F277W'
                else:
                    self.log.error('The SOSS extraction is implemented for the CLEAR or F277W filters only.')
                    self.log.error('extract_1d will be skipped.')
                    input_model.meta.cal_step.extract_1d = 'SKIPPED'
                    return input_model

                # Set the subarray mode being processed
                if input_model.meta.subarray.name == 'SUBSTRIP256':
                    self.log.info('Exposure is in the SUBSTRIP256 subarray.')
                    self.log.info('Traces 1 and 2 will be modelled and decontaminated before extraction.')
                    subarray = 'SUBSTRIP256'
                elif input_model.meta.subarray.name == 'FULL':
                    self.log.info('Exposure is in the FULL subarray.')
                    self.log.info('Traces 1 and 2 will be modelled and decontaminated before extraction.')
                    subarray = 'FULL'
                elif input_model.meta.subarray.name == 'SUBSTRIP96':
                    self.log.info('Exposure is in the SUBSTRIP96 subarray.')
                    self.log.info('Traces of orders 1 and 2 will be modelled but only order 1'
                                  ' will be decontaminated before extraction.')
                    subarray = 'SUBSTRIP96'
                else:
                    self.log.error('The SOSS extraction is implemented for the SUBSTRIP256,'
                                   ' SUBSTRIP96 and FULL subarray only.')
                    self.log.error('extract_1d will be skipped.')
                    input_model.meta.cal_step.extract_1d = 'SKIPPED'
                    return input_model

                # Load reference files.
                spectrace_ref_name = self.get_reference_file(input_model, 'spectrace')
                wavemap_ref_name = self.get_reference_file(input_model, 'wavemap')
                specprofile_ref_name = self.get_reference_file(input_model, 'specprofile')
                speckernel_ref_name = self.get_reference_file(input_model, 'speckernel')

                # Build SOSS kwargs dictionary.
                soss_kwargs = dict()
                soss_kwargs['threshold'] = self.soss_threshold
                soss_kwargs['n_os'] = self.soss_n_os
                soss_kwargs['tikfac'] = self.soss_tikfac
                soss_kwargs['width'] = self.soss_width
                soss_kwargs['bad_pix'] = self.soss_bad_pix
                soss_kwargs['transform'] = self.soss_transform
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
                    input_model,
                    spectrace_ref_name,
                    wavemap_ref_name,
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

                    input_model.close()

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
            else:
                # Get the reference file names
                if input_model.meta.exposure.type in extract.WFSS_EXPTYPES:
                    extract_ref = 'N/A'
                    self.log.info('No EXTRACT1D reference file will be used')
                else:
                    extract_ref = self.get_reference_file(input_model, 'extract1d')
                    self.log.info(f'Using EXTRACT1D reference file {extract_ref}')

                apcorr_ref = self.get_reference_file(input_model, 'apcorr') if self.apply_apcorr is True else 'N/A'

                if apcorr_ref == 'N/A':
                    self.log.info('APCORR reference file name is "N/A"')
                    self.log.info('APCORR will NOT be applied')
                else:
                    self.log.info(f'Using APCORR file {apcorr_ref}')

                result = extract.run_extract1d(
                    input_model,
                    extract_ref,
                    apcorr_ref,
                    self.smoothing_length,
                    self.bkg_fit,
                    self.bkg_order,
                    self.bkg_sigma_clip,
                    self.log_increment,
                    self.subtract_background,
                    self.use_source_posn,
                    self.center_xy,
                    was_source_model=False,
                )

                # Set the step flag to complete
                result.meta.cal_step.extract_1d = 'COMPLETE'

        input_model.close()

        return result
