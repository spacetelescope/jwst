"""Module for dark subtracting science data sets."""

import copy
import logging
import warnings

import numpy as np

from jwst.dark_current import dark_class

log = logging.getLogger(__name__)


__all__ = [
    "average_dark_frames_3d",
    "average_dark_frames_4d",
    "do_correction",
    "do_correction_data",
    "extrapolate_dark",
    "subtract_dark",
]


def do_correction(input_model, dark_model, dark_output=None):
    """
    Compute dark current subtraction.

    Convert models to classes to remove internal dependence on data models.
    After the creation of internal classes, the dark current subtraction is
    done.

    Parameters
    ----------
    input_model : data model
        Input model assumed to be like a JWST RampModel.

    dark_model : data model
        Input model assumed to be like a JWST DarkModel.

    dark_output : str
        A path to denote where to save the dark model, if desired.

    Returns
    -------
    output_data : `~jwst.dark_current.dark_class.ScienceData`
        Dark-subtracted science data

    averaged_dark : `~jwst.dark_current.dark_class.DarkData`
        New dark object with averaged frames
    """
    science_data = dark_class.ScienceData(input_model)
    dark_data = dark_class.DarkData(dark_model=dark_model)

    return do_correction_data(science_data, dark_data, dark_output)


def do_correction_data(science_data, dark_data, dark_output=None):
    """
    Execute all tasks for Dark Current Subtraction.

    Parameters
    ----------
    science_data : `~jwst.dark_current.dark_class.ScienceData`
        Science data to be corrected.

    dark_data : `~jwst.dark_current.dark_class.DarkData`
        Dark data.

    dark_output : str
        File name in which to optionally save averaged dark data.

    Returns
    -------
    output_data : `~jwst.dark_current.dark_class.ScienceData`
        Dark-subtracted science data.

    averaged_dark : `~jwst.dark_current.dark_class.DarkData`
        New dark object with averaged frames.
    """
    # Save some data params for easy use later
    sci_nints = science_data.data.shape[0]
    sci_ngroups = science_data.data.shape[1]
    sci_nframes = science_data.exp_nframes
    sci_groupgap = science_data.exp_groupgap

    if len(dark_data.data.shape) == 4:
        drk_nints = dark_data.data.shape[0]
        drk_ngroups = dark_data.data.shape[1]
    else:
        drk_nints = 1
        drk_ngroups = dark_data.data.shape[0]

    drk_nframes = dark_data.exp_nframes
    drk_groupgap = dark_data.exp_groupgap

    log.info(
        "Science data nints=%d, ngroups=%d, nframes=%d, groupgap=%d",
        sci_nints,
        sci_ngroups,
        sci_nframes,
        sci_groupgap,
    )
    log.info(
        "Dark data nints=%d, ngroups=%d, nframes=%d, groupgap=%d",
        drk_nints,
        drk_ngroups,
        drk_nframes,
        drk_groupgap,
    )

    # Check that the number of groups in the science data does not exceed
    # the number of groups in the dark current array.
    sci_total_frames = sci_ngroups * sci_nframes + (sci_ngroups - 1) * sci_groupgap
    drk_total_frames = drk_ngroups * drk_nframes + (drk_ngroups - 1) * drk_groupgap
    if sci_total_frames > drk_total_frames:
        log.warning("Not enough data in dark reference file to match to science data.")
        log.warning("Dark reference will be extrapolated to match length of science file.")

        frames_to_extrapolate = sci_total_frames - drk_total_frames
        # Find number of new groups required from above calculation of total frames.
        groups_to_extrapolate = np.ceil(
            (frames_to_extrapolate + drk_groupgap) / (drk_nframes + drk_groupgap)
        ).astype(int)
        extrapolate_dark(dark_data, groups_to_extrapolate)

    # Check that the value of nframes and groupgap in the dark
    # are not greater than those of the science data
    if drk_nframes > sci_nframes or drk_groupgap > sci_groupgap:
        log.warning(
            "The value of nframes or groupgap in the dark data is "
            "greater than that of the science data."
            "Input will be returned without subtracting dark current."
        )
        science_data.cal_step = "SKIPPED"
        out_data = copy.deepcopy(science_data)
        return out_data, None

    # Replace NaN's in the dark with zeros
    dark_data.data[np.isnan(dark_data.data)] = 0.0

    # Check whether the dark and science data have matching
    # nframes and groupgap settings.
    averaged_dark = None
    if sci_nframes == drk_nframes and sci_groupgap == drk_groupgap:
        # They match, so we can subtract the dark ref file data directly
        output_data = subtract_dark(science_data, dark_data)

        # If the user requested to have the dark file saved,
        # save the reference model as this file. This will
        # ensure consistency from the user's standpoint
        if dark_output is not None:
            averaged_dark = dark_data
            averaged_dark.save = True
            averaged_dark.output_name = dark_output

    else:
        # Create a frame-averaged version of the dark data to match
        # the nframes and groupgap settings of the science data.
        # If the data are from JWST/MIRI, the darks are integration-dependent
        # and we average them with a separate routine.

        if len(dark_data.data.shape) == 4:  # only MIRI uses 4-D darks
            averaged_dark = average_dark_frames_4d(
                dark_data, sci_nints, sci_ngroups, sci_nframes, sci_groupgap
            )
        else:
            averaged_dark = average_dark_frames_3d(
                dark_data, sci_ngroups, sci_nframes, sci_groupgap
            )

        # Save the frame-averaged dark data that was just created,
        # if requested by the user
        if dark_output is not None:
            averaged_dark.save = True
            averaged_dark.output_name = dark_output

        # Subtract the frame-averaged dark data from the science data
        output_data = subtract_dark(science_data, averaged_dark)

    output_data.cal_step = "COMPLETE"

    return output_data, averaged_dark


def average_dark_frames_3d(dark_data, ngroups, nframes, groupgap):
    """
    Average 3D dark frames.

    Averages the individual frames of data in a dark reference
    file to match the group structure of a science data set.
    This routine is not used for JWST/MIRI (see average_dark_frames_4d).

    Parameters
    ----------
    dark_data : `~jwst.dark_current.dark_class.DarkData`
        The input dark data.

    ngroups : int
        Number of groups in the science data set.

    nframes : int
        Number of frames per group in the science data set.

    groupgap : int
        Number of frames skipped between groups in the science data set.

    Returns
    -------
    avg_dark : `~jwst.dark_current.dark_class.DarkData`
        New dark object with averaged frames.
    """
    # Create a model for the averaged dark data
    dny = dark_data.data.shape[1]
    dnx = dark_data.data.shape[2]
    dims = (ngroups, dny, dnx)
    avg_dark = dark_class.DarkData(dims=dims)

    # Do a direct copy of the 2-d DQ array into the new dark
    avg_dark.groupdq = dark_data.groupdq

    # Loop over the groups of the input science data, copying or
    # averaging the dark frames to match the group structure
    start = 0

    for group in range(ngroups):
        end = start + nframes

        # If there's only 1 frame per group, just copy the dark frames
        if nframes == 1:
            log.debug("copy dark frame %d", start)
            avg_dark.data[group] = dark_data.data[start]

        # Otherwise generate new group by taking the mean of the SCI array nframes.
        else:
            log.debug("average dark frames %d to %d", start + 1, end)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "Mean of empty slice", RuntimeWarning)
                warnings.filterwarnings("ignore", "invalid value", RuntimeWarning)
                avg_dark.data[group] = dark_data.data[start:end].mean(axis=0)

        # Skip over unused frames
        start = end + groupgap

    # Reset some metadata values for the averaged dark
    avg_dark.exp_nframes = nframes
    avg_dark.exp_ngroups = ngroups
    avg_dark.exp_groupgap = groupgap

    return avg_dark


def average_dark_frames_4d(dark_data, nints, ngroups, nframes, groupgap):
    """
    Average 4D dark frames.

    Averages the individual frames of data in a dark reference
    file to match the group structure of a science data set.
    JWST/MIRI needs a separate routine because the darks are
    integration-dependent and hence 4D in shape, instead of 3D.
    An average dark is created for each integration.

    Parameters
    ----------
    dark_data : `~jwst.dark_current.dark_class.DarkData`
        The input dark data.

    nints : int
        Number of integrations in the science data.

    ngroups : int
        Number of groups in the science data set.

    nframes : int
        Number of frames per group in the science data set.

    groupgap : int
        Number of frames skipped between groups in the science data set.

    Returns
    -------
    avg_dark : `~jwst.dark_current.dark_class.DarkData`
        New dark object with averaged frames.
    """
    # Create a model for the averaged dark data
    dint = dark_data.data.shape[0]
    dny = dark_data.data.shape[2]
    dnx = dark_data.data.shape[3]
    dims = (dint, ngroups, dny, dnx)
    avg_dark = dark_class.DarkData(dims)

    # Do a direct copy of the 2-d DQ array into the new dark
    avg_dark.groupdq = dark_data.groupdq

    # check if the number of integrations in dark reference file
    # is less than science data.  If so then we only need to find the
    # average for dint integrations, otherwise we find the average for
    # nints.  (if science data only has 1 integration then there is no
    # need to average the second integration of the dark refefence file)
    num_ints = dint
    if dint > nints:
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
                log.debug("copy dark frame %d", start)
                avg_dark.data[it, group] = dark_data.data[it, start]

            # Otherwise average nframes into a new group: take the mean of
            # the SCI arrays and the quadratic sum of the ERR arrays.
            else:
                log.debug("average dark frames %d to %d", start + 1, end)
                avg_dark.data[it, group] = dark_data.data[it, start:end].mean(axis=0)

            # Skip over unused frames
            start = end + groupgap

    # Reset some metadata values for the averaged dark
    avg_dark.exp_nframes = nframes
    avg_dark.exp_ngroups = ngroups
    avg_dark.exp_groupgap = groupgap

    return avg_dark


def subtract_dark(science_data, dark_data):
    """
    Subtract dark current data from science arrays.

    Also updates data quality array based on DQ flags in the dark arrays.

    Parameters
    ----------
    science_data : `~jwst.dark_current.dark_class.ScienceData`
        The input science data.

    dark_data : `~jwst.dark_current.dark_class.DarkData`
        The dark current data.

    Returns
    -------
    output : `~jwst.dark_current.dark_class.ScienceData`
        Dark-subtracted science data.
    """
    # The integration start number is only needed for JWST/MIRI data.
    # It defaults to 1 if the keyword is not in the science data.
    int_start = 1 if science_data.exp_intstart is None else science_data.exp_intstart

    # Determine the number of integrations contained in the dark reference file
    dark_nints = dark_data.data.shape[0] if len(dark_data.data.shape) == 4 else 1

    log.debug(
        "subtract_dark: nints=%d, ngroups=%d, size=%d,%d",
        science_data.data.shape[0],
        science_data.data.shape[1],
        science_data.data.shape[2],
        science_data.data.shape[3],
    )

    # Create output as a copy of the input science data model
    output = copy.deepcopy(science_data)

    if len(dark_data.data.shape) == 4:
        # MIRI dark reference file has a DQ plane for each integration,
        # so we collapse the dark DQ planes into a single 2-D array
        darkdq = dark_data.groupdq[0, 0, :, :].copy()
        for i in range(1, dark_nints):
            darkdq = np.bitwise_or(darkdq, dark_data.groupdq[i, 0, :, :])
    else:
        # All other instruments have a single 2D dark DQ array
        darkdq = dark_data.groupdq

    # Propagate the dark DQ data into the science DQ data
    output.pixeldq = np.bitwise_or(science_data.pixeldq, darkdq)

    # Loop over all integrations in input science data
    for i in range(science_data.data.shape[0]):
        if len(dark_data.data.shape) == 4:  # MIRI data
            # Apply the first dark_nints-1 integrations from the dark ref file
            # to the first few science integrations. There's an additional
            # check of the starting integration number in case the science
            # data are segmented.
            #
            # else
            #
            # For science integrations beyond the number of
            # dark integrations, use the last dark integration
            dark_sci = (
                dark_data.data[i] if i < dark_nints and int_start == 1 else dark_data.data[-1]
            )
        else:
            # Use single-integration dark data
            dark_sci = dark_data.data

        # Loop over all groups in this integration
        for j in range(science_data.data.shape[1]):
            # subtract the dark from the science data
            output.data[i, j] -= dark_sci[j]

    return output


def extrapolate_dark(dark_data, ngroups):
    """
    Extrapolate a dark data array to cover the specified additional number of groups.

    Modifies the dark_data input in-place; dark_data.data may be 4-dimensional, as MIRI
    uses multi-integration darks, or 3-dimensional single integration darks for NIR instruments.

    Parameters
    ----------
    dark_data : `~stdatamodels.jwst.datamodels.DarkModel` or \
                `~stdatamodels.jwst.datamodels.DarkMIRIModel`
        The dark datamodel to be extrapolated.
    ngroups : int
        The number of groups to add to the dark array's ngroups.
    """

    def _extrapolate_int(arr, ngroups):
        """
        Extrapolate using rate derived from difference of last two groups.

        Parameters
        ----------
        arr : ndarray
            Dark array to extrapolate.

        ngroups : int
            The number of groups to add to the dark array's ngroups.

        Returns
        -------
        arr : ndarray
            Extrapolated array.
        """
        rate_arr = arr[-1, :, :] - arr[-2, :, :]
        new_groups = (
            np.broadcast_to(np.arange(1, ngroups + 1).reshape(-1, 1, 1), (ngroups, *rate_arr.shape))
            * rate_arr
            + arr[-1]
        )
        return np.concatenate((arr, new_groups))

    if len(dark_data.data.shape) == 4:
        nints, init_groups, nx, ny = dark_data.data.shape
        tmp_dark = dark_data.data
        dark_data.data = np.full((nints, init_groups + ngroups, nx, ny), np.nan)
        for i in range(nints):
            dark_data.data[i] = _extrapolate_int(tmp_dark[i], ngroups)
    else:
        dark_data.data = _extrapolate_int(dark_data.data, ngroups)

    dark_data.exp_ngroups += ngroups
