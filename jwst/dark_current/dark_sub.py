#
#  Module for dark subtracting science data sets
#

import numpy as np
import logging
from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, dark_model, dark_output=None):
    """
    Short Summary
    -------------
    Execute all tasks for Dark Current Subtraction

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    dark_model: dark model object
        dark data

    dark_output: string
        file name in which to optionally save averaged dark data

    Returns
    -------
    output_model: data model object
        dark-subtracted science data

    """

    # Save some data params for easy use later
    instrument = input_model.meta.instrument.name
    sci_nints = input_model.data.shape[0]
    sci_ngroups = input_model.data.shape[1]
    sci_nframes = input_model.meta.exposure.nframes
    sci_groupgap = input_model.meta.exposure.groupgap

    if instrument == 'MIRI':
        drk_nints = dark_model.data.shape[0]
        drk_ngroups = dark_model.data.shape[1]
    else:
        drk_nints = 1
        drk_ngroups = dark_model.data.shape[0]
    drk_nframes = dark_model.meta.exposure.nframes
    drk_groupgap = dark_model.meta.exposure.groupgap

    log.info(
        'Science data nints=%d, ngroups=%d, nframes=%d, groupgap=%d',
        sci_nints, sci_ngroups, sci_nframes, sci_groupgap
    )
    log.info(
        'Dark data nints=%d, ngroups=%d, nframes=%d, groupgap=%d',
        drk_nints, drk_ngroups, drk_nframes, drk_groupgap
    )

    # Check that the number of groups in the science data does not exceed
    # the number of groups in the dark current array.
    sci_total_frames = sci_ngroups * (sci_nframes + sci_groupgap)
    drk_total_frames = drk_ngroups * (drk_nframes + drk_groupgap)
    if sci_total_frames > drk_total_frames:
        log.warning(
            "Not enough data in dark reference file to match to "
            "science data."
        )
        log.warning("Input will be returned without subtracting dark current.")
        input_model.meta.cal_step.dark_sub = 'SKIPPED'
        return input_model.copy()

    # Check that the value of nframes and groupgap in the dark
    # are not greater than those of the science data
    if drk_nframes > sci_nframes or drk_groupgap > sci_groupgap:
        log.warning(
            "The value of nframes or groupgap in the dark data is "
            "greater than that of the science data."
            "Input will be returned without subtracting dark current."
        )
        input_model.meta.cal_step.dark_sub = 'SKIPPED'
        return input_model.copy()

    # Replace NaN's in the dark with zeros
    dark_model.data[np.isnan(dark_model.data)] = 0.0

    # Check whether the dark and science data have matching
    # nframes and groupgap settings.
    if sci_nframes == drk_nframes and sci_groupgap == drk_groupgap:

        # They match, so we can subtract the dark ref file data directly
        output_model = subtract_dark(input_model, dark_model)

        # If the user requested to have the dark file saved,
        # save the reference model as this file. This will
        # ensure consistency from the user's standpoint
        if dark_output is not None:
            log.info('Writing dark current data to %s', dark_output)
            dark_model.save(dark_output)

    else:

        # Create a frame-averaged version of the dark data to match
        # the nframes and groupgap settings of the science data.
        # If the data are from MIRI, the darks are integration-dependent and
        # we average them with a seperate routine.

        if instrument == 'MIRI':
            averaged_dark = average_MIRIdark_frames(
                dark_model, sci_nints, sci_ngroups, sci_nframes, sci_groupgap
            )
        else:
            averaged_dark = average_dark_frames(
                dark_model, sci_ngroups, sci_nframes, sci_groupgap
            )

        # Save the frame-averaged dark data that was just created,
        # if requested by the user
        if dark_output is not None:
            log.info('Writing dark current data to %s', dark_output)
            averaged_dark.save(dark_output)

        # Subtract the frame-averaged dark data from the science data
        output_model = subtract_dark(input_model, averaged_dark)

        averaged_dark.close()

    output_model.meta.cal_step.dark_sub = 'COMPLETE'

    return output_model


def average_dark_frames(input_dark, ngroups, nframes, groupgap):
    """
    Averages the individual frames of data in a dark reference
    file to match the group structure of a science data set.
    This routine is not used for MIRI (see average_MIRIdark_frames)

    Parameters
    ----------
    input_dark: dark data model
        the input dark data

    ngroups: int
        number of groups in the science data set

    nframes: int
        number of frames per group in the science data set

    groupgap: int
        number of frames skipped between groups in the science data set

    Returns
    -------
    avg_dark: dark data model
        New dark object with averaged frames

    """

    # Create a model for the averaged dark data
    dny = input_dark.data.shape[1]
    dnx = input_dark.data.shape[2]
    avg_dark = datamodels.DarkModel((ngroups, dny, dnx))
    avg_dark.update(input_dark)

    # Do a direct copy of the 2-d DQ array into the new dark
    avg_dark.dq = input_dark.dq

    # Loop over the groups of the input science data, copying or
    # averaging the dark frames to match the group structure
    start = 0

    for group in range(ngroups):
        end = start + nframes

        # If there's only 1 frame per group, just copy the dark frames
        if nframes == 1:
            log.debug('copy dark frame %d', start)
            avg_dark.data[group] = input_dark.data[start]
            avg_dark.err[group] = input_dark.err[start]

        # Otherwise average nframes into a new group: take the mean of
        # the SCI arrays and the quadratic sum of the ERR arrays.
        else:
            log.debug('average dark frames %d to %d', start + 1, end)
            avg_dark.data[group] = input_dark.data[start:end].mean(axis=0)
            avg_dark.err[group] = np.sqrt(np.add.reduce(
                input_dark.err[start:end]**2, axis=0)) / (end - start)

        # Skip over unused frames
        start = end + groupgap

    # Reset some metadata values for the averaged dark
    avg_dark.meta.exposure.nframes = nframes
    avg_dark.meta.exposure.ngroups = ngroups
    avg_dark.meta.exposure.groupgap = groupgap

    return avg_dark


def average_MIRIdark_frames(input_dark, nints, ngroups, nframes, groupgap):
    """
    Averages the individual frames of data in a dark reference
    file to match the group structure of a science data set.
    MIRI needs a seperate routine because the darks are integration dependent.
    We need an average dark for each dark integration.

    Parameters
    ----------
    input_dark: dark data model
        the input dark data

    nints: int
        number of integrations in the science data

    ngroups: int
        number of groups in the science data set

    nframes: int
        number of frames per group in the science data set

    groupgap: int
        number of frames skipped between groups in the science data set

    Returns
    -------
    avg_dark: dark data model
        New dark object with averaged frames

    """

    # Create a model for the averaged dark data
    dint = input_dark.data.shape[0]
    dny = input_dark.data.shape[2]
    dnx = input_dark.data.shape[3]
    avg_dark = datamodels.DarkMIRIModel((dint, ngroups, dny, dnx))
    avg_dark.update(input_dark)

    # Do a direct copy of the 2-d DQ array into the new dark
    avg_dark.dq = input_dark.dq

    # check if the number of integrations in dark reference file
    # is less than science data, if so then we only need to find the
    # average for num_ints integrations (if science data only has
    # 1 integration then there is no need to average the second integration
    # of the dark refefence file)
    num_ints = dint
    if(dint > nints):
        num_ints = nints

    # Loop over valid integrations in dark reference file
    # Loop over the groups of the input science data, copying or
    # averaging the dark frames to match the group structure

    for it in range(num_ints):
        start = 0
        for group in range(ngroups):
            end = start + nframes

            # If there's only 1 frame per group, just copy the dark frames
            if nframes == 1:
                log.debug('copy dark frame %d', start)
                avg_dark.data[it, group] = input_dark.data[it, start]
                avg_dark.err[it, group] = input_dark.err[it, start]

            # Otherwise average nframes into a new group: take the mean of
            # the SCI arrays and the quadratic sum of the ERR arrays.
            else:
                log.debug('average dark frames %d to %d', start + 1, end)
                avg_dark.data[it, group] = input_dark.data[it, start:end].mean(axis=0)
                avg_dark.err[it, group] = np.sqrt(np.add.reduce(
                    input_dark.err[it, start:end]**2, axis=0)) / (end - start)

        # Skip over unused frames
            start = end + groupgap

    # Reset some metadata values for the averaged dark
    avg_dark.meta.exposure.nframes = nframes
    avg_dark.meta.exposure.ngroups = ngroups
    avg_dark.meta.exposure.groupgap = groupgap

    return avg_dark


def subtract_dark(input, dark):
    """
    Subtracts dark current data from science arrays, combines
    error arrays in quadrature, and updates data quality array based on
    DQ flags in the dark arrays.

    Parameters
    ----------
    input: data model object
        the input science data

    dark: dark model object
        the dark current data

    Returns
    -------
    output: data model object
        dark-subtracted science data

    """

    instrument = input.meta.instrument.name
    if instrument == 'MIRI':
        dark_nints = dark.data.shape[0]
    else:
        dark_nints = 1

    log.debug("subtract_dark: nints=%d, ngroups=%d, size=%d,%d",
              input.data.shape[0], input.data.shape[1],
              input.data.shape[2], input.data.shape[3])

    # Create output as a copy of the input science data model
    output = input.copy()

    if instrument == 'MIRI':
        # MIRI dark reference file has a DQ plane for each integration,
        # so we collapse the dark DQ planes into a single 2-D array
        darkdq = dark.dq[0, 0, :, :].copy()
        for i in range(1, dark_nints):
            darkdq = np.bitwise_or(darkdq, dark.dq[i, 0, :, :])
    else:
        # All other instruments have a single 2D dark DQ array
        darkdq = dark.dq

    # Combine the dark and science DQ data
    output.pixeldq = np.bitwise_or(input.pixeldq, darkdq)

    # loop over all integrations and groups in input science data
    for i in range(input.data.shape[0]):

        if instrument == 'MIRI':
            if i < dark_nints:
                dark_int = dark.data[i]
            else:
                dark_int = dark.data[dark_nints - 1]

        for j in range(input.data.shape[1]):
            # subtract the SCI arrays
            if instrument == 'MIRI':
                output.data[i, j] -= dark_int[j]
            else:
                output.data[i, j] -= dark.data[j]

            # combine the ERR arrays in quadrature
            # NOTE: currently stubbed out until ERR handling is decided
            #output.err[i,j] = np.sqrt(
            #           output.err[i,j]**2 + dark.err[j]**2)

    return output
