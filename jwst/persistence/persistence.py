#
#  Module for correcting for persistence

import datetime
import logging
import math

import numpy as np
from stdatamodels.jwst import datamodels
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

    persistence_time : 
    persistence_array : 
    persistence_dnu : 
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
            return self.output_obj, None, None, skipped

        (nints, ngroups, nrows, ncols) = shape

        # XXX Probably change this somewhat
        # Initialize the optional output.
        '''
        if self.save_persistence:
            self.output_pers = self.output_obj.copy()
            self.output_pers.data[:, :, :, :] = 0.0
            self.output_pers.pixeldq = None
            self.output_pers.groupdq = None
            self.output_pers.err = None
        else:
            self.output_pers = None
        '''

        # XXX Use different start time. This may not be the start time of the
        #     exposure, but the beginning of the tasking, which could include
        #     setup time and other time outside the actual exposure.
        etime = datetime.datetime.fromisoformat(self.output_obj.meta.observation.date_beg)
        epoch_time = etime.timestamp()
        integration_time = self.output_obj.meta.exposure.integration_time
        group_time = self.output_obj.meta.exposure.group_time
        sat_array = np.zeros(shape=(nrows, ncols), dtype=np.uint32)

        for integ in range(nints):
            if self.output_obj.int_times is not None and len(self.output_obj.int_times) > integ:
                # The index into int_times is integration, then time type. [integ][1] is the
                # int_start_MJD_UTC type for the integration 'integ'.
                current_time = mjd_to_epoch(self.output_obj.int_times[integ][1])
            else:
                current_time = epoch_time + integ * integration_time
            print(f"Integration {integ} start time (epoch): {current_time}, mjd: {epoch_to_mjd(current_time)}")
            for group in range(ngroups):
                current_time = current_time + group_time
                if self.persistence_time is not None:
                    self.process_persistence_flagging(current_time, sat_array, integ, group)

        if self.save_persistence:
            self.output_pers.update(self.output_obj, only="PRIMARY")

        return self.output_obj, None, self.output_pers, skipped

    def process_persistence_flagging(self, current_time, sat_array, integ, group):
        """
        Flag groups that are within a persistence window.

        The structure of the persistence_array is as follows:
            1. Zero entries indicate no persistence flagging.
            2. Non-zero entries indicate the end time of the persistence flagging window.

        Since a non-zero entry indicates a persistence flagging window has been found for
        that pixel. The entry is the epoch time of the end of that window for a pixel. If
        the current_time is less than an epoch time, then that group gets flagged with
        persistence. Any non-zero entry that occurs before the current_time indicates
        the persistence window for that pixel has ended, so the persistece_array entry
        for that pixel will be set to zero.

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
        # persistence window has ended for that pixel.
        self.persistence_array[self.persistence_array < current_time] = 0.0

        # Calculate any first saturation points
        gdq_plane = self.output_obj.groupdq[integ, group, :, :]
        sat_loc = np.bitwise_and(gdq_plane, dqflags.group["SATURATED"])
        sat_array[sat_loc > 0] += 1

        # Add first saturation points and the end window time to persistence_array
        self.persistence_array[sat_array==1] = current_time + self.persistence_time

        # Set persistence flag for any group in persistence window
        if self.persistence_dnu:
            flag = dqflags.group["DO_NOT_USE"] | dqflags.pixel["PERSISTENCE"]
        else:
            flag = dqflags.pixel["PERSISTENCE"]

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