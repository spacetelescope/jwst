from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
#
#  Module for correcting for persistence

import math
import numpy as np
import logging
from .. import datamodels
from .. datamodels import dqflags

from astropy.io import fits             # xxx test debug

log = logging.getLogger()
log.setLevel(logging.DEBUG)


def no_NaN(input_model, fill_value,
           zap_nan=False, zap_zero=False):
    """Replace NaNs and/or zeros with a fill value."""

    mask = None

    if zap_nan:
        mask = (np.isnan(input_model.data))

    if zap_zero:
        if mask is None:
            mask = (input_model.data == 0.)
        else:
            mask = np.logical_or(mask, (input_model.data == 0.))

    if mask is None or mask.sum(dtype=np.intp) == 0:
        return input_model
    else:
        temp = input_model.copy()
        temp.data[mask] = fill_value
        return temp


class DataSet():
    """
    Input dataset to which persistence will be applied

    Parameters
   ----------
    """
    def __init__(self, output_obj, traps_filled_model,
                 trap_density_model, traps_model,
                 pers_sat_model):
        """
        Short Summary
        -------------
        Set attributes

        Parameters
        ----------
        output_obj:
            copy of input data model object

        traps_filled_model: cube model, or None
            Image of trap state.  There will be one or more image planes,
            each of which corresponds to a trap "family," i.e. a set of
            traps with similar capture and decay properties.
            If this is None, the state will be initialized to an array
            of zeros, indicating that there are no filled traps.

        trap_density_model: image model
            Image of the total number of traps per pixel.

        traps_model: traps model
            Table giving parameters for traps.

        pers_sat_model: saturation model
            Persistence saturation limit reference file.
        """

        self.output_obj = output_obj
        self.traps_filled = traps_filled_model
        self.trap_density = trap_density_model
        self.traps_model = traps_model
        self.persistencesat = pers_sat_model


    def do_all(self):
        """
        Short Summary
        -------------
        Execute all tasks for persistence correction

        Parameters
        ----------

        Returns
        -------
        self.output_obj: ramp model
            persistence-corrected input data
        """

        # Initial value, indicates that processing was done successfully.
        skipped = False

        detector = self.output_obj.meta.instrument.detector

        shape = self.output_obj.data.shape
        if len(shape) != 4:
            log.warning("Don't understand shape %s of input data, "
                        "skipping ...", str(shape))
            skipped = True
            return (self.output_obj, None, skipped)

        # Read the table of capture and decay parameters.
        par = self.get_parameters()
        nfamilies = len(par[0])

        # Note that this might be a subarray.
        persistence = self.output_obj.data[0, 0, :, :] * 0.
        (nints, ngroups, ny, nx) = shape
        t_group = self.output_obj.meta.exposure.group_time
        # Time of one integration, in seconds.  This is the duration,
        # not necessarily the effective integration time.
        t_int = t_group * float(ngroups)

        if not self.traps_filled:
            have_traps_filled = False
            if detector == "MIRI":
                (det_ny, det_nx) = (1024, 1032)
            else:
                (det_ny, det_nx) = (2048, 2048)
            self.traps_filled = datamodels.CubeModel(
                        data=np.zeros((nfamilies, det_ny, det_nx),
                                      dtype=self.output_obj.data.dtype))
        else:
            have_traps_filled = True

        self.traps_filled = no_NaN(self.traps_filled, 0., zap_nan=True)
        self.trap_density = no_NaN(self.trap_density, 0., zap_nan=True)
        self.persistencesat = no_NaN(self.persistencesat, 1.e7, zap_nan=True)
        log.debug("xxx do_all:  k = 1, traps_filled:  %g %g %g",
                  self.traps_filled.data[1].min(),
                  self.traps_filled.data[1].mean(),
                  self.traps_filled.data[1].max())
        log.debug("xxx do_all:  persistencesat.mean = %g",
                  self.persistencesat.data.mean())

        log.debug("xxx do_all:  trap_density:  %g %g %g",
                  self.trap_density.data.min(),
                  self.trap_density.data.mean(),
                  self.trap_density.data.max())

        if have_traps_filled:   # was there an actual traps_filled file?
            log.debug("xxx do_all:  yes, we have a traps_filled file")
            # Decrease traps_filled by the number of traps that decayed
            # in the time (to_start) from the end of the traps_filled file
            # to the start of the current exposure.
            to_start = (self.output_obj.meta.exposure.start_time -
                        self.traps_filled.meta.exposure.end_time) * 86400.
            log.debug("xxx do_all:  to_start = %g seconds", to_start)
            for k in range(nfamilies):
                log.debug("xxx do_all:  k = %d, traps_filled:  %g %g %g",
                          k, self.traps_filled.data[k].min(),
                             self.traps_filled.data[k].mean(),
                             self.traps_filled.data[k].max())
                decay_param_k = self.get_decay_param(par, k)
                decay = self.compute_decay(self.traps_filled.data[k],
                                           decay_param_k, to_start)
                log.debug("xxx do_all:  k = %d, decay:  %g %g %g",
                          k, decay.min(), decay.mean(), decay.max())
                self.traps_filled.data[k, :, :] -= decay
                log.debug("xxx do_all:  traps_filled:  %g %g %g",
                          self.traps_filled.data[k, :, :].min(),
                          self.traps_filled.data[k, :, :].mean(),
                          self.traps_filled.data[k, :, :].max())

        # If the science image is a subarray, extract matching sections of
        # reference images.
        # We don't need to extract a subarray from self.traps_filled (which
        # is full-frame and 3-D, one image plane for each trap family), but
        # we can use this tuple of slice objects later, to update the
        # self.traps_filled data.
        save_slice = self.get_slice(self.traps_filled, self.output_obj)
        if self.ref_matches_sci(self.traps_filled, save_slice):
            is_subarray = False
            save_slice = None
        else:
            is_subarray = True
        log.debug("xxx is_subarray = %s", str(is_subarray))

        slc = self.get_slice(self.trap_density, self.output_obj)
        if not self.ref_matches_sci(self.trap_density, slc):
            self.trap_density = self.get_subarray(self.trap_density, slc)

        slc = self.get_slice(self.persistencesat, self.output_obj)
        if not self.ref_matches_sci(self.persistencesat, slc):
            self.persistencesat = self.get_subarray(self.persistencesat, slc)

        # These are for saving info in output files for testing, and
        # save_decayed is also used for temporary storage (it's populated
        # in the loop over groups and used in the loop over integrations).
        # They're 3-D arrays, one plane for each trap family.  If the
        # input is a subarray, so are these.
        save_filled = np.zeros((nfamilies, shape[-2], shape[-1]),
                               dtype=self.output_obj.data.dtype)
        save_cr_filled = save_filled.copy()
        save_decayed = save_filled.copy()

        # self.traps_filled will be updated with each integration, to
        # account for charge capture and decay of traps.
        for integ in range(nints):
            delta_t = 0.
            # slope has to be computed early in the loop over integrations,
            # before the data are modified by subtracting persistence.
            # The slope is needed for computing charge captures.
            (grp_slope, slope) = self.compute_slope(integ)
            for group in range(ngroups):
                # Time from beginning of integration to end of current group.
                delta_t += t_group
                persistence[:, :] = 0.          # initialize
                for k in range(nfamilies):
                    decay_param_k = self.get_decay_param(par, k)
                    # This will be the size of the full detector.
                    if group == 0:              # xxx test debug
                        log.debug("xxx do_all/compute_decay:  "
                                  "k = %d, traps_filled:  %g %g %g",
                                  k, self.traps_filled.data[k].min(),
                                     self.traps_filled.data[k].mean(),
                                     self.traps_filled.data[k].max())
                    decayed = self.compute_decay(self.traps_filled.data[k],
                                                 decay_param_k, delta_t)
                    if group == 0:              # xxx test debug
                        log.debug("xxx do_all/compute_decay:  "
                                  "decayed:  %g %g %g",
                                  decayed.min(), decayed.mean(), decayed.max())
                    if is_subarray:
                        persistence += decayed[save_slice[0], save_slice[1]]
                    else:
                        persistence += decayed
                    if group == ngroups - 1:
                        if is_subarray:
                            save_decayed[k, :, :] = decayed[save_slice[0],
                                                            save_slice[1]]
                        else:
                            save_decayed[k, :, :] = decayed
                # Persistence was computed in DN.
                self.output_obj.data[integ, group, :, :] -= persistence

            # Update traps_filled with the number of traps that captured
            # a charge or decayed during the current integration.
            for k in range(nfamilies):
                capture_param_k = self.get_capture_param(par, k)
                # This may be a subarray.
                filled = self.compute_capture(capture_param_k,
                                              self.trap_density.data,
                                              slope, t_int)
                cr_filled = self.delta_fcn_capture(capture_param_k,
                                                   self.trap_density.data,
                                                   grp_slope,
                                                   integ, ngroups, t_group)
                save_filled[k, :, :] = filled.copy()
                save_cr_filled[k, :, :] = cr_filled.copy()
                filled += cr_filled
                del cr_filled
                if is_subarray:
                    self.traps_filled.data[k, save_slice[0], save_slice[1]] += \
                                (filled - save_decayed[k, :, :])
                else:
                    self.traps_filled.data[k, :, :] += \
                                (filled - save_decayed[k, :, :])
                log.debug("yyy k = %d:  %g, %g, %g, %g, %g, %g --> %g",
                          k, capture_param_k[0], capture_param_k[1],
                          capture_param_k[2],
                          self.trap_density.data[900, 900],
                          slope[900, 900], t_int,
                          filled[900, 900])

            # Write to files for testing and debugging.         # xxx
            fits.writeto("filled_{}.fits".format(integ),
                         data=save_filled, overwrite=True)
            fits.writeto("cr_filled_{}.fits".format(integ),
                         data=save_cr_filled, overwrite=True)
            fits.writeto("decayed_{}.fits".format(integ),
                         data=save_decayed, overwrite=True)

        # Update the start and end times (and other stuff) in the
        # traps_filled image to the times for the current exposure.
        self.traps_filled.update(self.output_obj, only="PRIMARY")

        # The file name of traps_filled is now the name of the science file
        # that was passed as input to this step.  This is good, since we're
        # going to write traps_filled out, and the output name should be
        # related to the output science file name.  However, let's append
        # a dummy suffix, to preserve the real last suffix.
        n = len(self.traps_filled.meta.filename)
        self.traps_filled.meta.filename = \
                self.traps_filled.meta.filename[0:n - 5] + "_zzzz.fits"

        return (self.output_obj, self.traps_filled, skipped)

        # Set meta.filename to the name of the input file, so that the
        # save_model() method of PersistenceStep will use that name (but
        # with a different suffix) as the output traps_filled file name,
        # to (hopefully clearly) link them.
        self.traps_filled.meta.filename = self.output_obj.meta.filename


    def get_slice(self, ref, sci):
        """Find the 2-D slice for a reference file.

        Parameters
        ----------
        ref: data model
            A reference image.

        sci: data model
            The science data.

        Returns
        -------
        tuple
            A two-element tuple of slice objects, the Y and X slices that
            can be used to extract a subarray from the reference file.
        """

        ref_shape = ref.shape
        ref_nx = ref_shape[-1]
        ref_ny = ref_shape[-2]
        sci_shape = sci.shape
        sci_nx = sci_shape[-1]
        sci_ny = sci_shape[-2]

        # These are limits of slices.
        sci_x1 = sci.meta.subarray.xstart - 1
        sci_x2 = sci_x1 + sci_nx
        sci_y1 = sci.meta.subarray.ystart - 1
        sci_y2 = sci_y1 + sci_ny
        log.debug("sci xstart=%d, xstop=%d, ystart=%d, ystop=%d",
                  sci_x1, sci_x2, sci_y1, sci_y2)

        ref_x1 = ref.meta.subarray.xstart - 1
        ref_y1 = ref.meta.subarray.ystart - 1
        if ref_x1 is None:
            ref_x1 = sci_x1
        if ref_y1 is None:
            ref_y1 = sci_y1
        ref_x2 = ref_x1 + ref_nx
        ref_y2 = ref_y1 + ref_ny
        log.debug("ref xstart=%d, xstop=%d, ystart=%d, ystop=%d",
                  ref_x1, ref_x2, ref_y1, ref_y2)

        # Compute the slicing indexes
        xstart = sci_x1 - ref_x1
        ystart = sci_y1 - ref_y1
        xstop = xstart + sci_nx
        ystop = ystart + sci_ny
        log.debug("ref slice %d:%d, %d:%d", ystart, ystop, xstart, xstop)

        # Check for errors in the slice indexes
        if (xstart < 0 or ystart < 0 or
            xstop > ref.data.shape[-1] or ystop > ref.data.shape[-2]):
            log.error("Science and reference file arrays not compatible")
            raise ValueError("Can't extract matching subarray from "
                             "reference data")

        return (slice(ystart, ystop), slice(xstart, xstop))


    def ref_matches_sci(self, ref, slc):
        """Test whether ref and sci are full-frame or the same subarray.

        Parameters
        ----------
        ref: data model
            A reference image.

        slc: tuple of two slice objects
            The Y and X slices that can be used to extract a subarray from
            the reference file.  This is returned by function get_slice.

        Returns
        -------
        """

        slc_test_y = slice(0, ref.shape[-2])
        slc_test_x = slice(0, ref.shape[-1])

        if slc[0] == slc_test_y and slc[1] == slc_test_x:
            return True
        else:
            return False


    def get_subarray(self, ref, slc):
        """Extract a subarray from a reference file.

        Parameters
        ----------
        ref: data model
            A reference image.

        slc: tuple of two slice objects
            The Y and X slices that can be used to extract a subarray from
            the reference file.  This is returned by function get_slice.

        Returns
        -------
        refsub: data model
            `refsub` is the subarray extracted from `ref`.  `refsub` will
            sometimes be 2-D and other times 3-D.
        """

        # The xxx
        refsub = ref.copy()
        refsub.data = ref.data[..., slc[0], slc[1]].copy()
        if ref.__hasattr__(err):
            refsub.err = ref.err[..., slc[0], slc[1]].copy()
        if ref.__hasattr__(dq):
            refsub.dq = ref.dq[..., slc[0], slc[1]].copy()

        return refsub


    def get_parameters(self):
        """Read capture and decay parameters from a reference table.

        Parameters
        ----------

        Returns
        -------
        tuple (par0, par1, par2, par3) of 1-D ndarrays
            par0, par1, par2 are trap capture columns
            par3 is a trap decay column
        """

        data = self.traps_model.traps_table
        par0 = data["capture0"].copy()
        par1 = data["capture1"].copy()
        par2 = data["capture2"].copy()
        par3 = data["decay_param"].copy()

        return (par0, par1, par2, par3)


    def compute_slope(self, integ):
        """Compute an estimate of the slope of the ramp for each pixel.

        We need a value for the slope that will not include cosmic-ray
        jumps, and groups that are flagged as saturated must also not
        contribute to the slope.  The approach will be to find the
        difference between adjacent groups, set the difference to a very
        large value if the group was saturated, then sort the differences
        along the axis for the ramp.  All the saturated groups will then
        be at the (high) end.  Cosmic-ray-affected groups should be just
        below the saturated groups (unless the effect of a CR was very
        small, in which case it won't matter where it is).  We can therefore
        ignore all the CR jumps and saturated values just by knowing how
        many of them there are.  The average of the remaining differences
        is the slope.

        Parameters
        ----------
        integ: int
            The number (index) of the current integration.

        Returns
        -------
        (grp_slope, slope): tuple of 2-D ndarrays
            Both arrays give the ramp slope at each pixel, but the units
            differ.  `grp_slope` is the ramp slope in units of counts per
            group, while `slope` is the ramp slope in units of (counts per
            persistence saturation limit) per second.  The reason for
            keeping both arrays is that the persistence saturation limit
            (which could be, for example, of order 90000) can differ from
            one pixel to another.
        """

        (_, ngroups, ny, nx) = self.output_obj.shape
        if ngroups < 2:
            # These values are supposed to be arrays, but the way they're
            # used (in compute_capture and delta_fcn_capture), they will
            # either broadcast to the array shape (in compute_capture) or
            # they will not actually be referenced when ngroups < 2.
            return (0., 0.)

        gdqflags = dqflags.group

        data = self.output_obj.data[integ, :, :, :]
        gdq = self.output_obj.groupdq[integ, :, :, :]

        # This assumes that the jump step has already been run, so that
        # CR jumps will have been flagged in the groupdq extension.
        # n_cr is a 2-D array of the number of cosmic-ray hits per pixel.
        cr = (np.bitwise_and(gdq, gdqflags["JUMP_DET"]) > 0).astype(np.int32)
        n_cr = cr.sum(axis=0)
        del cr

        # In the absence of saturation and jumps, these differences would
        # be approximately constant for a given pixel.
        diff = data[1:, :, :] - data[0:-1, :, :]
        max_diff = diff.max()
        # Larger than any actual data value
        huge_diff = max(2. * max_diff, 1.e5)
        del max_diff

        # This assumes that the saturation step has already been run.
        # n_sat is a 2-D array of the number of saturated groups per pixel.
        s_mask = (np.bitwise_and(gdq[1:, :, :], gdqflags["SATURATED"]) > 0)
        diff[s_mask] = huge_diff
        sat = s_mask.astype(np.int32)
        n_sat = sat.sum(axis=0)
        del s_mask, sat

        sdiff = np.sort(diff, axis=0)
        # Upper limit of good data (slice notation, i.e. exclusive).
        upper = diff.shape[0] - (n_cr + n_sat)
        really_bad = np.where(upper < 0)
        upper[really_bad] = 0           # for the comparison with indx below
        del really_bad
        temp = np.arange(ngroups - 1, dtype=np.int32)
        temp2 = np.repeat(temp, nx * ny)
        indx = temp2.reshape((ngroups - 1, ny, nx))
        mask = np.where(indx >= upper)
        sdiff[mask] = 0.                # zero values won't affect the sum
        del temp, temp2, indx, mask
        bad = np.where(upper <= 0)
        upper[bad] = 1                  # so we can divide by upper
        # xxx begin test debug
        print("xxx sdiff:", flush=True)
        for j in range(944, 949):
            for i in range(940, 945):
                print("{} {}  {}".format(j, i, sdiff[0:5, j, i]), flush=True)
        # xxx end test debug
        # grp_slope has units of DN / group
        grp_slope = np.sum(sdiff, axis=0) / upper.astype(np.float32)
        grp_slope[bad] = 0.
        del bad
        print("xxx slope[944:949, 940:945]:", flush=True)       # xxx debug
        print(grp_slope[944:949, 940:945], flush=True)      # xxx test debug

        # units = (DN / persistence_saturation_limit) / second
        slope = grp_slope / (self.persistencesat.data *
                             self.output_obj.meta.exposure.group_time)
        log.debug("xxx compute_slope:  grp_slope = %g %g %g",
                  grp_slope.min(), grp_slope.mean(), grp_slope.max())
        log.debug("xxx compute_slope:  slope = %g %g %g",
                  slope.min(), slope.mean(), slope.max())
        fits.writeto("slope_{}.fits".format(integ),
                     data=slope, overwrite=True)

        return (grp_slope, slope)


    def get_capture_param(self, par, k):
        """Extract capture parameters for the current trap family.

        Parameters
        ----------
        par: tuple of ndarray
            These were read from the traps reference table.  Each element
            of the tuple is a column from the table.  Each row of the
            table is for a different trap family.

        k: int
            Index of the current trap family

        Returns
        -------
        tuple of float
            This includes just the capture parameters, and only the values
            for the current trap family.
        """

        (par0, par1, par2) = par[0:3]

        return (par0[k], par1[k], par2[k])


    def get_decay_param(self, par, k):
        """Extract decay parameter(s) for the current trap family.

        Parameters
        ----------
        par: tuple of ndarray
            These were read from the traps reference table.  Each element
            of the tuple is a column from the table.  Each row of the
            table is for a different trap family.

        k: int
            Index of the current trap family

        Returns
        -------
        float
            This is just the decay parameter, and only the value for the
            current trap family.
        """

        par3 = par[3]

        return par3[k]


    def compute_capture(self, capture_param_k, trap_density, slope, delta_t):
        """Compute the number of traps that will be filled in time delta_t.

        This is based on Michael Regan's predictrampcapture3.pro.

        Parameters
        ----------
        capture_param_k: tuple
            Three values read from a reference table.  These will be from
            three separate columns but just one row; the row corresponds
            to the current trap family.  (The _k in the variable name
            refers to the index of the trap family.)

        trap_density: 2-D ndarray
            Image of the total number of traps per pixel.

        slope: 2-D ndarray
            Array of the slope of the ramp at each pixel.  The slope was
            computed from the pixel values that were not saturated and were
            not affected by jumps, based on flags in the groupdq extension.
            The unit is (fraction of the persistence saturation limit)
            per second.

        delta_t: float
            The time interval (unit = second) over which the charge
            capture is to be computed.

        Returns
        -------
        2-D ndarray
            The computed traps_filled at the end of the integration.
        """

        (par0, par1, par2) = capture_param_k

        t1 = delta_t**2 / 2.
        t2 = (-delta_t / par1) * math.exp(par1 * delta_t)
        t3 = (1.0 / par1**2) * math.exp(par1 * delta_t)
        t4 = -1.0 / par1**2
        traps_filled = (trap_density * slope**2 *
                        (par0 * (t1 + t2 + t3 + t4) + delta_t**2 * par2 / 2.))

        return traps_filled


    def delta_fcn_capture(self, capture_param_k, trap_density,
                          grp_slope, integ, ngroups, t_group):
        """Compute capture due to cosmic-ray jumps.

        If a cosmic-ray hit was in group number g (meaning that
        data[integ, g, y, x] was higher than expected), then
        delta_t = (ngroups - g - 0.5) * t_group
        is the time from the CR-affected group to the end of the
        integration, assuming that the CR hit in the middle (timewise)
        of the group.
        cr_filled = trap_density * par0 * jump * (1. - exp(par1 * delta_t))

        Parameters
        ----------
        capture_param_k: tuple
            Three values read from a reference table.  These will be from
            three separate columns but just one row; the row corresponds
            to the current trap family.  (The _k in the variable name
            refers to the index of the trap family.)

        trap_density: 2-D ndarray
            Image of the total number of traps per pixel.

        grp_slope: 2-D ndarray
            Array of the slope of the ramp at each pixel.  The slope was
            computed from the pixel values that were not saturated and were
            not affected by jumps, based on flags in the groupdq extension.
            The unit is counts (DN) per group.

        integ: int
            Integration number.

        ngroups: int
            Total number of groups in the integration.

        t_group: float
            The time (seconds) from the start of one group to the start
            of the next group.

        Returns
        -------
        2-D ndarray
            The computed cr_filled at the end of the integration.
        """

        (par0, par1, _) = capture_param_k
        # cr_filled will be incremented group-by-group, depending on
        # where cosmic rays were found in each group.
        cr_filled = np.zeros_like(trap_density)
        data = self.output_obj.data[integ, :, :, :]
        gdq = self.output_obj.groupdq[integ, :, :, :]
        gdqflags = dqflags.group

        # If there's a CR hit in the first group, we can't determine its
        # amplitude, so skip the first group.  Loop over all subsequent
        # groups, checking (via the groupdq extension) for cosmic-ray hits.
        for group in range(1, ngroups):
            cr_flag = np.where(np.bitwise_and(gdq[group, :, :],
                                              gdqflags['JUMP_DET']) > 0)
            if len(cr_flag[0]) > 0:
                delta_t = (float(ngroups - group) - 0.5) * t_group
                # z and z_prev are index arrays.
                z = np.empty(len(cr_flag[0]), dtype=np.intp)
                z[:] = group
                z_prev = z - 1
                # jump is a 1-D array, just for the CRs in the current group.
                jump = ((data[z, cr_flag[0], cr_flag[1]] -
                         data[z_prev, cr_flag[0], cr_flag[1]]) -
                        grp_slope[cr_flag])
                jump = np.where(jump < 0., 0., jump)
                cr_filled[cr_flag] += (trap_density[cr_flag] * par0 * jump *
                                       (1. - math.exp(par1 * delta_t)))

        return cr_filled


    def compute_decay(self, traps_filled, decay_param, delta_t):
        """Compute the number of trap decays.

        This is based on Michael Regan's trapdecay.pro.

        Parameters
        ----------
        traps_filled: 2-D ndarray
            This is an image of the number of filled traps in each pixel
            for the current trap family.

        decay_param: float
            The decay parameter.  This is negative, but otherwise it's
            the e-folding time for trap decay for the current trap family.

        delta_t: float
            The time interval (unit = second) over which the trap decay
            is to be computed.

        Returns
        -------
        decayed: 2-D ndarray
            Image of the computed number of trap decays for each pixel,
            for the current trap family.
        """

        if decay_param > 0.:
            log.warning("compute_decay:  decay_param is %g; "
                        "a negative value was expected.", decay_param)
            decay_param = -decay_param
            log.warning("compute_decay:  Using decay_param = %g instead.",
                        decay_param)
        decayed = np.where(traps_filled > 0.,
                           traps_filled * (1. - math.exp(decay_param * delta_t)),
                           0.)
        return decayed
