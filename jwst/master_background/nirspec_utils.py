import logging

from jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


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

    log.info('Interpolating 1D master background to 2D slitlets')

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

    log.info('Creating master background from background slitlets')

    # Copy dedicated background slitlets to a temporary model
    bkg_model = datamodels.MultiSlitModel()
    bkg_model.update(input_model)
    slits = []
    for slit in input_model.slits:
        if "background" in slit.source_name:
            log.info(f'Using slitlet {slit.source_name}')
            slits.append(slit)

    if len(slits) == 0:
        log.warning('No background slitlets found; skipping master bkg correction')
        return None

    bkg_model.slits.extend(slits)

    # Apply resample_spec and extract_1d to all background slitlets
    resamp = resample_spec_step.ResampleSpecStep.call(bkg_model)
    x1d = extract_1d_step.Extract1dStep.call(resamp)

    # Call combine_1d to combine the 1D background spectra
    log.info('Combining background spectra into master background')
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
    input_model.data *= (pl_point / pl_uniform)

    return input_model
