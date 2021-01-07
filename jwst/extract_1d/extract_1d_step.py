from ..stpipe import Step
from .. import datamodels
from . import extract


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
        Default is `poly`.

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
        parameter can be overriden (set to False) internally by the step.

    apply_apcorr : bool
        Switch to select whether or not to apply an APERTURE correction during
        the Extract1dStep. Default is True
    """

    spec = """
    smoothing_length = integer(default=None)  # background smoothing size
    bkg_fit = option("poly", "mean", "median", default="poly")  # background fitting type
    bkg_order = integer(default=None, min=0)  # order of background polynomial fit
    log_increment = integer(default=50)  # increment for multi-integration log messages
    subtract_background = boolean(default=None)  # subtract background?
    use_source_posn = boolean(default=None)  # use source coords to center extractions?
    apply_apcorr = boolean(default=True)  # apply aperture corrections?
    """

    reference_file_types = ['extract1d', 'apcorr']

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

        was_source_model = False                 # default value
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
            # NRS_BRIGHTOBJ mode
            self.log.debug('Input is a SlitModel')
        else:
            self.log.error(f'Input is a {str(type(input_model))}, ')
            self.log.error('which was not expected for extract_1d')
            self.log.error('extract_1d will be skipped.')
            input_model.meta.cal_step.extract_1d = 'SKIPPED'
            return input_model


        #______________________________________________________________________
        # Do the extraction for ModelContainer - this might only be WFSS data
        if isinstance(input_model, datamodels.ModelContainer):

            # This is the branch  WFSS data take
            if len(input_model) > 1:
                self.log.debug(f"Input contains {len(input_model)} items")

                #--------------------------------------------------------------
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
                        self.log_increment,
                        self.subtract_background,
                        self.use_source_posn,
                        was_source_model=was_source_model
                    )
                    # Set the step flag to complete
                    result.meta.cal_step.extract_1d = 'COMPLETE'

                #--------------------------------------------------------------
                # Data is a ModelContainer but is not WFSS
                    result.meta.filetype = '1d spectrum'
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
                            self.log_increment,
                            self.subtract_background,
                            self.use_source_posn,
                            was_source_model=was_source_model,
                            )
                        # Set the step flag to complete in each MultiSpecModel
                        temp.meta.cal_step.extract_1d = 'COMPLETE'
                        temp.meta.filetype = '1d spectrum'
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
                    self.log_increment,
                    self.subtract_background,
                    self.use_source_posn,
                    was_source_model=was_source_model,
                )

                # Set the step flag to complete
                result.meta.cal_step.extract_1d = 'COMPLETE'
                result.meta.filetype = '1d spectrum'
            else:
                self.log.error('Input model is empty;')
                self.log.error('extract_1d will be skipped.')
                return input_model

        #______________________________________________________________________
        # Data that is not a ModelContainer (IFUCube and other single models)
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
                self.log_increment,
                self.subtract_background,
                self.use_source_posn,
                was_source_model=False,
            )

            # Set the step flag to complete
            result.meta.cal_step.extract_1d = 'COMPLETE'
            result.meta.filetype = '1d spectrum'

        input_model.close()

        return result
