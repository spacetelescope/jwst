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

    bkg_order : int or None
        If not None, a polynomial with order `bkg_order` will be fit to
        each column (or row if the dispersion direction is horizontal)
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

    apply_nod_offset : bool or None
        If True, the source and background positions specified in the
        reference file (or the default position, if there is no reference
        file) will be shifted to account for nod and/or dither offset.  If
        None (the default), the value in the reference file will be used,
        or it will be set to True if it is not specified in the reference
        file.  This offset is determined by finding the location in the data
        corresponding to the target position (keywords TARG_RA and TARG_DEC).
        For NIRSpec fixed-slit or MOS data, there must be different target
        coordinates for each slit for which a nod correction might be
        needed, and this is not implemented yet in extract_1d.
        It also doesn't make sense to apply a nod/dither offset for an
        extended target, so this flag can internally be overridden (set to
        False) for extended targets.
    """

    spec = """
    # Boxcar smoothing width for background regions.
    smoothing_length = integer(default=None)
    # Order of polynomial fit to one column (or row if the dispersion
    # direction is vertical) of background regions.
    bkg_order = integer(default=None, min=0)
    # Log a progress message when processing multi-integration data.
    log_increment = integer(default=50)
    # Flag indicating whether the background should be subtracted.
    subtract_background = boolean(default=None)
    # If True, the locations of the target and background regions will be
    # shifted to correct for the computed nod/dither offset.
    # Currently this offset is not applied for NIRSpec fixed-slit or
    # MOS (MSA) data), or for WFSS data.
    apply_nod_offset = boolean(default=None)
    """

    reference_file_types = ['extract1d']

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
            self.log.error('Input is a %s,', str(type(input_model)))
            self.log.error('which was not expected for extract_1d.')
            self.log.error('extract_1d will be skipped.')
            input_model.meta.cal_step.extract_1d = 'SKIPPED'
            return input_model

        # Do the extraction
        if isinstance(input_model, datamodels.ModelContainer):
            if len(input_model) > 1:
                self.log.debug("Input contains %d items", len(input_model))
                result = datamodels.ModelContainer()
                for model in input_model:
                    # This is a flag for do_extract1d.
                    if model.meta.exposure.type in extract.WFSS_EXPTYPES:
                        ref_file = 'N/A'
                        self.log.info('No EXTRACT1D reference file '
                                      'will be used')
                    else:
                        # Get the reference file name
                        ref_file = self.get_reference_file(model, 'extract1d')
                        self.log.info('Using EXTRACT1D reference file %s',
                                      ref_file)
                    temp = extract.run_extract1d(model, ref_file,
                                                 self.smoothing_length,
                                                 self.bkg_order,
                                                 self.log_increment,
                                                 self.subtract_background,
                                                 self.apply_nod_offset,
                                                 was_source_model=was_source_model)
                    # Set the step flag to complete in each MultiSpecModel
                    temp.meta.cal_step.extract_1d = 'COMPLETE'
                    result.append(temp)
                    del temp
            elif len(input_model) == 1:
                if input_model[0].meta.exposure.type in extract.WFSS_EXPTYPES:
                    ref_file = 'N/A'
                    self.log.info('No EXTRACT1D reference file will be used')
                else:
                    # Get the reference file name for the one model in input
                    ref_file = self.get_reference_file(input_model[0],
                                                            'extract1d')
                    self.log.info('Using EXTRACT1D reference file %s',
                                  ref_file)
                result = extract.run_extract1d(input_model[0], ref_file,
                                               self.smoothing_length,
                                               self.bkg_order,
                                               self.log_increment,
                                               self.subtract_background,
                                               self.apply_nod_offset,
                                               was_source_model=was_source_model)
                # Set the step flag to complete
                result.meta.cal_step.extract_1d = 'COMPLETE'
            else:
                self.log.error('Input model is empty;')
                self.log.error('extract_1d will be skipped.')
                return input_model
        else:
            # Get the reference file name
            if input_model.meta.exposure.type in extract.WFSS_EXPTYPES:
                ref_file = 'N/A'
                self.log.info('No EXTRACT1D reference file will be used')
            else:
                ref_file = self.get_reference_file(input_model, 'extract1d')
                self.log.info('Using EXTRACT1D reference file %s', ref_file)
            result = extract.run_extract1d(input_model, ref_file,
                                           self.smoothing_length,
                                           self.bkg_order,
                                           self.log_increment,
                                           self.subtract_background,
                                           self.apply_nod_offset,
                                           was_source_model=False)
            # Set the step flag to complete
            result.meta.cal_step.extract_1d = 'COMPLETE'

        input_model.close()

        return result
