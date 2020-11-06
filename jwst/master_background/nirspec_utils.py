import logging

from jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def apply_master_background(source_model, bkg_model, inverse=False):
    """Subtract 2D master background signal from each source
    slitlet in the input MultiSlitModel.

    Parameters
    ----------
    source_model : `~jwst.datamodels.MultiSlitModel`
        The input data model containing all source slit instances.

    bkg_model : `~jwst.datamodels.MultiSlitModel`
        The data model containing 2D background slit instances.

    inverse : boolean
        Invert the math operations used to apply the background.

    Returns
    -------
    output_model: `~jwst.datamodels.MultiSlitModel`
        The output background-subtracted data model.
    """
    from .master_background_step import subtract_2d_background

    if inverse:
        log.info('Adding master background from each MOS source slitlet')
        bkg = bkg_model.copy()
        for slit in bkg.slits:
            slit.data *= -1.0
    else:
        log.info('Subtracting master background from each MOS source slitlet')
        bkg = bkg_model

    # This does a one-to-one subtraction of the data in each background
    # slit from the data in the corresponding source slit (i.e. the
    # two MultiSlitModels must have matching numbers of slit instances).
    # This may be changed in the future to only do the subtraction from
    # a certain subset of source slits.
    output_model = subtract_2d_background(source_model, bkg)

    return output_model


def map_to_science_slits(input_model, master_bkg):
    """Interpolate 1D master background spectrum to the 2D space
    of each source slitlet in the input MultiSlitModel.

    Parameters
    ----------
    input_model : `~jwst.datamodels.MultiSlitModel`
        The input data model containing all slit instances.

    master_bkg : `~jwst.datamodels.CombinedSpecModel`
        The 1D master background spectrum.

    Returns
    -------
    output_model: `~jwst.datamodels.MultiSlitModel`
        The output data model containing background signal.
    """
    from .expand_to_2d import expand_to_2d

    log.info('Interpolating 1D master background to all MOS 2D slitlets')

    # Loop over all input slits, creating 2D master background to
    # match each 2D slitlet cutout
    output_model = expand_to_2d(input_model, master_bkg)

    return output_model


def create_background_from_multislit(input_model):
    """Create a 1D master background spectrum from a set of
    calibrated background MOS slitlets in the input
    MultiSlitModel.

    Parameters
    ----------
    input_model : `~jwst.datamodels.MultiSlitModel`
        The input data model containing all slit instances.

    Returns
    -------
    master_bkg: `~jwst.datamodels.CombinedSpecModel`
        The 1D master background spectrum created from the inputs.
    """
    from ..resample import resample_spec_step
    from ..extract_1d import extract_1d_step
    from ..combine_1d.combine1d import combine_1d_spectra

    log.info('Creating MOS master background from background slitlets')

    # Copy dedicated background slitlets to a temporary model
    bkg_model = datamodels.MultiSlitModel()
    bkg_model.update(input_model)
    slits = []
    for slit in input_model.slits:
        if "background" in slit.source_name:
            log.info(f'Using background slitlet {slit.source_name}')
            slits.append(slit)

    if len(slits) == 0:
        log.warning('No background slitlets found; skipping master bkg correction')
        return None

    bkg_model.slits.extend(slits)

    # Apply resample_spec and extract_1d to all background slitlets
    log.info('Applying resampling and 1D extraction to background slits')
    resamp = resample_spec_step.ResampleSpecStep.call(bkg_model)
    x1d = extract_1d_step.Extract1dStep.call(resamp)

    # Call combine_1d to combine the 1D background spectra
    log.info('Combining 1D background spectra into master background')
    master_bkg = combine_1d_spectra(x1d, exptime_key='exposure_time')

    del bkg_model
    del resamp
    del x1d

    return master_bkg


def correct_nrs_ifu_bkg(input_model):
    """Apply point source vs. uniform source pathloss adjustments
    to a NIRSpec IFU 2D master background array.

    Parameters
    ----------
    input_model : `~jwst.datamodels.IFUImageModel`
        The input background data.

    Returns
    -------
    input_model : `~jwst.datamodels.IFUIMAGEModel`
        An updated (in place) version of the input with the data
        replaced by the corrected 2D background.
    """

    log.info('Applying point source pathloss updates to IFU background')

    # Try to load the appropriate pathloss correction arrays
    try:
        pl_point = input_model.getarray_noinit('pathloss_point')
    except AttributeError:
        log.warning('Pathloss_point array not found in input')
        log.warning('Skipping pathloss background updates')
        return input_model

    try:
        pl_uniform = input_model.getarray_noinit('pathloss_uniform')
    except AttributeError:
        log.warning('Pathloss_uniform array not found in input')
        log.warning('Skipping pathloss background updates')
        return input_model

    # Apply the corrections
    input_model.data *= (pl_uniform / pl_point)

    return input_model


def correct_nrs_fs_bkg(input_model, primary_slit):
    """Apply point source vs. uniform source corrections
    to a NIRSpec Fixed-Slit 2D master background array.

    Parameters
    ----------
    input_model : `~jwst.datamodels.SlitModel`
        The input background data.

    primary_slit : bool
        Is this the primary slit in the exposure?

    Returns
    -------
    input_model : `~jwst.datamodels.SlitModel`
        An updated (in place) version of the input with the data
        replaced by the corrected 2D background.
    """
    log.info('Applying point source updates to FS background')

    # Try to load the appropriate pathloss correction arrays
    if 'pathloss_point' in input_model.instance:
        pl_point = getattr(input_model, 'pathloss_point')
    else:
        log.warning('pathloss_point array not found in input')
        log.warning('Skipping background updates')
        return input_model

    if 'pathloss_uniform' in input_model.instance:
        pl_uniform = getattr(input_model, 'pathloss_uniform')
    else:
        log.warning('pathloss_uniform array not found in input')
        log.warning('Skipping background updates')
        return input_model

    if primary_slit:
        # If processing the primary slit, we also need flatfield and
        # photom correction arrays
        if 'flatfield_point' in input_model.instance:
            ff_point = getattr(input_model, 'flatfield_point')
        else:
            log.warning('flatfield_point array not found in input')
            log.warning('Skipping background updates')
            return input_model

        if 'flatfield_uniform' in input_model.instance:
            ff_uniform = getattr(input_model, 'flatfield_uniform')
        else:
            log.warning('flatfield_uniform array not found in input')
            log.warning('Skipping background updates')
            return input_model

        if 'photom_point' in input_model.instance:
            ph_point = getattr(input_model, 'photom_point')
        else:
            log.warning('photom_point array not found in input')
            log.warning('Skipping background updates')
            return input_model

        if 'photom_uniform' in input_model.instance:
            ph_uniform = getattr(input_model, 'photom_uniform')
        else:
            log.warning('photom_uniform array not found in input')
            log.warning('Skipping background updates')
            return input_model

        # Apply the corrections for the primary slit
        input_model.data *= (pl_uniform / pl_point) * \
                            (ff_uniform / ff_point) * \
                            (ph_point / ph_uniform)

    else:
        # Apply the corrections for secondary slits
        input_model.data *= (pl_uniform / pl_point)

    return input_model
