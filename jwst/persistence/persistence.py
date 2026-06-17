#
#  Module for correcting for persistence

import logging

import numpy as np
from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)

SCALEFACTOR = 2.0
"""This factor is to account for the difference in gain for charges freed
from traps, compared with photon-generated charges.
"""

__all__ = ["DataSet"]


class DataSet:
    """
    Input dataset to which persistence will be applied.

    Attributes
    ----------
    output_obj : JWST data model
        A copy of the input model.  This will be modified in-place.

    save_persistence : bool
        If True, the persistence that was subtracted will be written to an
        output file.

    persistence_time : int
        The number of seconds for a persistence window to persist.

    persistence_array : string or None
        If not None, then it is the path to a file containing a persistence array.

    persistence_dnu : boolean
        When flagging PERSISTENCE, if true, then flag as DO_NOT_USE as well.
    """

    def __init__(
        self,
        output_obj,
        save_persistence,
        persistence_time=None,
        persistence_array=None,
        persistence_dnu=None,
    ):
        """
        Assign values to attributes.

        Parameters
        ----------
        output_obj : JWST data model
            Copy of input data model object

        save_persistence : bool
            If True, the persistence that was subtracted will be written
            to an output file.
        """
        log.debug("save_persistence = %s", str(save_persistence))

        self.output_obj = output_obj
        self.save_persistence = save_persistence
        self.output_pers = None

        self.persistence_time = persistence_time
        self.persistence_array = persistence_array
        self.persistence_dnu = persistence_dnu

    def do_all(self):
        """
        Execute all tasks for persistence correction.

        Returns
        -------
        output_obj : data model
            The persistence-corrected science data, a RampModel object.

        output_pers :  data model or None
            A RampModel object, giving the value of persistence that
            was subtracted from each pixel of each group of each
            integration.

        skipped : bool
            This will be True if the step has been skipped.
        """
        # Initial value, indicates that processing was done successfully.
        skipped = False

        shape = self.output_obj.data.shape
        if len(shape) != 4:
            log.warning(f"Don't understand shape {shape} of input data, skipping...")
            skipped = True
            return self.output_obj, skipped

        (nints, ngroups, nrows, ncols) = shape

        epoch_time = mjd_to_epoch(self.output_obj.meta.exposure.start_time)
        integration_time = self.output_obj.meta.exposure.integration_time
        group_time = self.output_obj.meta.exposure.group_time
        sat_array = np.zeros(shape=(nrows, ncols), dtype=np.uint32)

        for integ in range(nints):
            if self.output_obj.int_times is not None and len(self.output_obj.int_times) > integ:
                # The index into int_times is integration, then time type. [integ][1] is the
                # int_start_MJD_UTC type for the integration 'integ'.
                current_time = mjd_to_epoch(self.output_obj.int_times[integ]["int_start_MJD_UTC"])
            else:
                # Exposure start time is used if int_times is not available.
                current_time = epoch_time + integ * integration_time

            for group in range(ngroups):
                current_time = current_time + group_time
                if self.persistence_time is not None:
                    self.process_persistence_flagging(current_time, sat_array, integ, group)
            sat_array[:, :] = 0

        return self.output_obj, skipped

    def process_persistence_flagging(self, current_time, sat_array, integ, group):
        """
        Flag groups that are within a persistence window.

        The structure of the persistence_array is as follows:
            1. Zero entries indicate no persistence flagging.
            2. Non-zero entries indicate the end time of the persistence flagging window.

        Since a non-zero entry indicates a persistence flagging window has been found for
        that pixel. The entry is the epoch time of the end of that window for a pixel. To
        check to see if any groups are within a persistence flagging window, it is first
        decided if there is a non-zero entry in the persistence_array for that pixel. If
        non-zero, make sure the current_time is between the start and end times of the
        persistence window in the persistence_array.

        Additionally, the current time may be the beginning of a persistence window if it
        is determined to be the first group in a ramp to be saturated, in which case the
        end time of the window is current_time + persistence_time.

        Parameters
        ----------
        current_time : float
            The epoch time of the current integration and group.

        sat_array : ndarray
            The saturation count for each pixel.

        integ : int
            The current integration being processed.

        group : int
            The current group being processed.
        """
        # Any persistence_array time earlier than current time gets set to zero. The
        # persistence window has ended for that pixel because the current time is
        # after the end of the persistence window.
        self.persistence_array[self.persistence_array < current_time] = 0.0

        # Calculate any first saturation points. Any group found to be the first
        # saturated group in a ramp is the beginning of a persistence window.
        gdq_plane = self.output_obj.groupdq[integ, group, :, :]
        sat_loc = np.bitwise_and(gdq_plane, dqflags.group["SATURATED"])
        sat_array[sat_loc > 0] += 1
        self.persistence_array[sat_array == 1] = current_time + self.persistence_time

        # This prevents 'backwards flagging'.
        # Subtracting the persistence_time gives the beginning of the window.
        # If the current time occurs before this time, then the current group is
        #     outside the window and subtracting it will be a positive number.
        # Raise an exception if a backwards flagging situation arrives, as this is
        #     an invalid state.
        start_plane = self.persistence_array - (self.persistence_time + current_time)
        start_plane[self.persistence_array == 0.0] = 0.0
        if np.any(start_plane > 0.0):
            raise ValueError("Invalid persistence array, due to backwards flagging.")

        # Set persistence flag for any group in persistence window
        if self.persistence_dnu:
            flag = dqflags.group["DO_NOT_USE"] | dqflags.group["PERSISTENCE"]
        else:
            flag = dqflags.group["PERSISTENCE"]

        gdq_plane[self.persistence_array > 0.0] |= flag
        self.output_obj.groupdq[integ, group, :, :] = gdq_plane


def mjd_to_epoch(mjd):
    """
    Convert Modified Julian Date to Unix epoch time.

    Parameters
    ----------
    mjd : float
        Modified Julian Date

    Returns
    -------
    float
        Unix epoch time (seconds since January 1, 1970 00:00:00 UTC)
    """
    # MJD 40587.0 corresponds to Unix epoch (January 1, 1970 00:00:00 UTC)
    # 86400 seconds per day
    epoch_time = (mjd - 40587.0) * 86400.0
    return epoch_time


def epoch_to_mjd(epoch):
    """
    Convert Unix epoch time to Modified Julian Date.

    Parameters
    ----------
    epoch : float
        Unix epoch time (seconds since January 1, 1970 00:00:00 UTC)

    Returns
    -------
    float
        Modified Julian Date
    """
    # MJD 40587.0 corresponds to Unix epoch (January 1, 1970 00:00:00 UTC)
    # 86400 seconds per day
    mjd = epoch / 86400.0 + 40587.0
    return mjd
