#
#  Module for correcting for persistence

import math
import numpy as np
import logging
from .. import datamodels
from .. datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

SCALEFACTOR = 2.
"""This factor is to account for the difference in gain for charges freed
from traps, compared with photon-generated charges.
"""


def no_NaN(input_model, fill_value,
           zap_nan=False, zap_zero=False):
    """Replace NaNs and/or zeros with a fill value.

    Parameters
    ----------
    input_model : JWST data model
        The input will typically be a reference file model.

    fill_value : float
        Use this value to replace NaNs and/or zeros.

    zap_nan : bool
        If True, replace NaNs in a copy of `input_model.data`.
        The default is False.

    zap_zero : bool
        If True, replace zeros in a copy of `input_model.data`.
        The default is False.

    Returns
    -------
    JWST data model
        A copy of `input_model` with NaNs and/or zeros in the `data`
        attribute replaced with `fill_value`.
    """

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
    """Input dataset to which persistence will be applied

    Attributes
    ----------
    output_obj : JWST data model
        A copy of the input model.  This will be modified in-place.

    traps_filled : JWST data model, `TrapsFilledModel`
        The trap state at some time prior to the current exposure

    flag_pers_cutoff : float or None
        If not None, pixels will be flagged if the value of persistence
        that was subtracted is larger than `flag_pers_cutoff`.

    save_persistence : bool
        If True, the persistence that was subtracted will be written to an
        output file.

    output_pers : JWST data model or None
        If `save_persistence` is True, the amount of persistence that was
        subtracted will be copied to the `data` attribute of a data model
        and written to a file.

    trap_density : JWST data model, `TrapDensityModel`
        Reference file, giving the total number of traps per pixel.

    trappars_model : JWST data model, `TrapParsModel`
        Reference file (table), giving parameters for traps.

    persistencesat : JWST data model, `PersistenceSatModel`
        Persistence saturation limit (full well) reference file.

    tframe : float
        The frame time, in seconds.

    tgroup : float
        The group time, in seconds.

    ngroups : int
        The number of groups per integration.

    nframes : int
        The number of frames per group.

    groupgap : int
        The number of dropped frames.  Currently not used.

    nresets : int
        The number of resets (frames) at the beginning of each integration.
    """

    def __init__(self, output_obj, input_traps_filled,
                 flag_pers_cutoff, save_persistence,
                 trap_density_model, trappars_model,
                 persat_model):
        """Assign values to attributes.

        Parameters
        ----------
        output_obj : JWST data model
            copy of input data model object

        input_traps_filled : cube model or None
            Image of trap state.  There will be one or more image planes,
            each of which corresponds to a trap "family," i.e. a set of
            traps with similar capture and decay properties.
            If this is None, the state will be initialized to an array
            of zeros, indicating that there are no filled traps.

        flag_pers_cutoff : float or None
            If not None, pixels will be flagged (with what? xxx) if the
            value of persistence that was subtracted is larger than
            `flag_pers_cutoff`.

        save_persistence : bool
            If True, the persistence that was subtracted will be written
            to an output file.

        trap_density_model : image model
            Image (reference file) of the total number of traps per pixel.

        trappars_model : traps model
            Table (reference file) giving parameters for traps.

        persat_model : persistence saturation model
            Persistence saturation limit (full well) reference file.
        """

        log.debug("input_traps_filled = %s", str(input_traps_filled))
        log.debug("flag_pers_cutoff = %g", flag_pers_cutoff)
        log.debug("save_persistence = %s", str(save_persistence))

        self.output_obj = output_obj
        self.traps_filled = input_traps_filled
        self.flag_pers_cutoff = flag_pers_cutoff
        self.save_persistence = save_persistence
        self.output_pers = None

        self.trap_density = trap_density_model
        self.trappars_model = trappars_model
        self.persistencesat = persat_model

        # These will be populated from metadata.
        self.tframe = 0.
        self.tgroup = 0.
        self.ngroups = 0
        self.nframes = 0
        self.groupgap = 0
        self.nresets = 0

    def do_all(self):
        """Execute all tasks for persistence correction

        Returns
        -------
        output_obj : data model
            The persistence-corrected science data, a RampModel object.

        traps_filled : data model or None
            A TrapsFilledModel object, giving the number of traps that
            are filled in each pixel at the end of the exposure; there
            will be one plane for each trap family.

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

        # Read the table of capture and decay parameters.
        par = self.get_parameters()
        nfamilies = len(par[0])
        if nfamilies <= 0:
            log.error("The trappars reference table is empty!")

        # Note that this might be a subarray.
        persistence = np.zeros((shape[-2], shape[-1]), dtype=np.float64)
        (nints, ngroups, ny, nx) = shape
        t_group = self.output_obj.meta.exposure.group_time

        # The trap density image is full-frame, so use it to get the shape of
        # a full-frame array.  traps_filled is also full-frame, but we can't
        # rely on that because there might not be an input traps_filled file.
        (det_ny, det_nx) = self.trap_density.data.shape

        is_subarray = (shape[-2] != det_ny or shape[-1] != det_nx)
        if is_subarray:
            log.debug("The input is a subarray.")
        else:
            log.debug("The input is not a subarray.")

        if not self.traps_filled:
            have_traps_filled = False
            self.traps_filled = datamodels.TrapsFilledModel(
                data=np.zeros((nfamilies, det_ny, det_nx),
                              dtype=self.output_obj.data.dtype)
            )
            self.traps_filled.meta.subarray.xstart = 1
            self.traps_filled.meta.subarray.ystart = 1
            self.traps_filled.meta.subarray.xsize = det_nx
            self.traps_filled.meta.subarray.ysize = det_ny
        else:
            have_traps_filled = True
            self.traps_filled = no_NaN(self.traps_filled, 0., zap_nan=True)

        self.trap_density = no_NaN(self.trap_density, 0., zap_nan=True)
        self.persistencesat = no_NaN(self.persistencesat, 1.e7, zap_nan=True)

        if have_traps_filled:   # was there an actual traps_filled file?
            # Decrease traps_filled by the number of traps that decayed
            # in the time (to_start) from the end of the traps_filled file
            # to the start of the current exposure.
            # Note that to_start includes the reset (if any) at the
            # beginning of the exposure, because meta.exposure.start_time
            # is the time when the actual exposure started, not the time
            # at the beginning of the reset prior to the exposure.
            # Decays during the time covered by the reset (again, if any)
            # between integrations will be taken care of in the loop
            # over integrations.
            to_start = (self.output_obj.meta.exposure.start_time
                        - self.traps_filled.meta.exposure.end_time) * 86400.
            log.debug("Decay time for previous traps-filled file = %g s",
                      to_start)
            for k in range(nfamilies):
                decay_param_k = self.get_decay_param(par, k)
                decay = self.compute_decay(self.traps_filled.data[k],
                                           decay_param_k, to_start)
                self.traps_filled.data[k, :, :] -= decay
                del decay

        """
        These will be full-frame:
            self.traps_filled           (nfamilies, det_ny, det_nx)
            self.trap_density (before extracting subarray)
            self.persistencesat (before extracting subarray)
            decayed                     (nfamilies, det_ny, det_nx)
            decayed_in_group            (det_ny, det_nx)

        These will be subarrays if the input object is a subarray:
            self.output_obj             (nints, ngroups, ny, nx)
            self.trap_density (after extracting subarray)
            self.persistencesat (after extracting subarray)
            persistence
            self.output_pers
            filled
            cr_filled
        """

        # If the science image is a subarray, extract matching sections of
        # reference images unless they match the science image.

        # We don't need to extract a subarray from self.traps_filled (which
        # is full-frame and 3-D, one image plane for each trap family), but
        # we can use this tuple of slice objects later, e.g. to update the
        # self.traps_filled data.
        save_slice = self.get_slice(self.traps_filled, self.output_obj)

        slc = self.get_slice(self.trap_density, self.output_obj)
        if not self.ref_matches_sci(self.trap_density, slc):
            self.trap_density = self.get_subarray(self.trap_density, slc)

        # Trap density is used for computing the number of traps that
        # capture a charge.  If some pixels in the trap_density reference
        # file are flagged as bad, set the corresponding data values to
        # zero, so the computed number of traps will also be zero.
        if hasattr(self.trap_density, "dq"):
            mask = (np.bitwise_and(self.trap_density.dq,
                                   dqflags.pixel["DO_NOT_USE"]) > 0)
            self.trap_density.data[mask] = 0.
            del mask

        slc = self.get_slice(self.persistencesat, self.output_obj)
        if not self.ref_matches_sci(self.persistencesat, slc):
            self.persistencesat = self.get_subarray(self.persistencesat, slc)

        # Initialize the optional output.
        if self.save_persistence:
            self.output_pers = self.output_obj.copy()
            self.output_pers.data[:, :, :, :] = 0.
            self.output_pers.pixeldq = None
            self.output_pers.groupdq = None
            self.output_pers.err = None
        else:
            self.output_pers = None

        # Buffer for accumulating the number of decayed traps from the
        # start of an integration to the current group (in the loop below).
        decayed = np.zeros((nfamilies, det_ny, det_nx), dtype=np.float64)

        # self.traps_filled will be updated with each integration, to
        # account for charge capture and decay of traps.
        filled = -1                             # just to ensure that it exists
        for integ in range(nints):
            self.get_group_info(integ)          # self.tgroup, etc.
            decayed[:, :, :] = 0.               # initialize
            # slope has to be computed early in the loop over integrations,
            # before the data are modified by subtracting persistence.
            # The slope is needed for computing charge captures.
            (grp_slope, slope) = self.compute_slope(integ)
            for group in range(ngroups):
                persistence[:, :] = 0.          # initialize
                for k in range(nfamilies):
                    decay_param_k = self.get_decay_param(par, k)
                    # Compute and subtract the decays during the reset.
                    # Decays during the reset at the beginning of the
                    # first integration have already been accounted for.
                    if integ > 0 and group == 0 and self.nresets > 0:
                        reset_time = self.tframe * self.nresets
                        decay_during_reset = \
                            self.compute_decay(self.traps_filled.data[k],
                                               decay_param_k, reset_time)
                        self.traps_filled.data[k, :, :] -= decay_during_reset
                    # Decays during current group, for current trap family.
                    decayed_in_group = \
                        self.compute_decay(self.traps_filled.data[k],
                                           decay_param_k, t_group)
                    # Cumulative decay to the end of the current group.
                    decayed[k, :, :] += decayed_in_group
                    self.traps_filled.data[k, :, :] -= decayed_in_group
                    if is_subarray:
                        persistence += decayed[k, save_slice[0], save_slice[1]]
                    else:
                        persistence += decayed[k, :, :]
                    del decayed_in_group

                # Persistence was computed in DN.
                self.output_obj.data[integ, group, :, :] -= persistence
                if self.save_persistence:
                    self.output_pers.data[integ, group, :, :] = \
                        persistence.copy()
                if persistence.max() >= self.flag_pers_cutoff:
                    mask = (persistence >= self.flag_pers_cutoff)
                    self.output_obj.pixeldq[mask] |= dqflags.pixel['PERSISTENCE']

            # Update traps_filled with the number of traps that captured
            # a charge during the current integration.
            for k in range(nfamilies):
                capture_param_k = self.get_capture_param(par, k)
                # This may be a subarray.
                filled = self.predict_capture(capture_param_k,
                                              self.trap_density.data,
                                              integ, grp_slope, slope)
                if is_subarray:
                    self.traps_filled.data[k, save_slice[0],
                                           save_slice[1]] += filled
                else:
                    self.traps_filled.data[k, :, :] += filled

        del filled

        # Update the start and end times (and other stuff) in the
        # traps_filled image to the times for the current exposure.
        self.traps_filled.update(self.output_obj, only="PRIMARY")

        # meta.filename of traps_filled is now the name of the science file
        # that was passed as input to this step.  This is good, since we're
        # going to write traps_filled out, and the output name should be
        # related to the output science file name.

        if self.save_persistence:
            self.output_pers.update(self.output_obj, only="PRIMARY")

        return self.output_obj, self.traps_filled, self.output_pers, skipped

    def get_slice(self, ref, sci):
        """Find the 2-D slice for a reference file.

        Parameters
        ----------
        ref : data model
            A reference image.

        sci : data model
            The science data.

        Returns
        -------
        tuple of two slice objects
            The elements are the Y and X slices that are intended to be
            used to extract a subarray from the reference file.
        """

        sci_shape = sci.shape
        sci_nx = sci_shape[-1]
        sci_ny = sci_shape[-2]

        # These are limits of slices.
        sci_x1 = sci.meta.subarray.xstart - 1
        sci_y1 = sci.meta.subarray.ystart - 1

        ref_x1 = ref.meta.subarray.xstart - 1
        ref_y1 = ref.meta.subarray.ystart - 1
        if ref_x1 is None:
            ref_x1 = sci_x1
        if ref_y1 is None:
            ref_y1 = sci_y1

        # Compute the slicing indexes
        xstart = sci_x1 - ref_x1
        ystart = sci_y1 - ref_y1
        xstop = xstart + sci_nx
        ystop = ystart + sci_ny

        # Check for errors in the slice indexes
        if (xstart < 0 or ystart < 0 or
                xstop > ref.data.shape[-1] or ystop > ref.data.shape[-2]):
            log.error("Science and reference file arrays not compatible")
            raise ValueError("Can't extract matching subarray from "
                             "reference data")

        return slice(ystart, ystop), slice(xstart, xstop)

    def ref_matches_sci(self, ref, slc):
        """Test whether ref and sci cover the same area of the detector.

        Parameters
        ----------
        ref : data model
            Reference data.

        slc : tuple of two slice objects
            The Y and X slices that can be used to extract a subarray from
            the reference file.  This was returned by function `get_slice`.

        Returns
        -------
        bool
            True if both are full-frame or if they are the same subarray;
            False otherwise.
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
        ref : JWST data model
            Reference data.

        slc : tuple of two slice objects
            The Y and X slices that can be used to extract a subarray from
            the reference file.  This was returned by function `get_slice`.

        Returns
        -------
        refsub : data model
            `refsub` is a copy of `ref`, but the data attribute in `refsub`
            includes only the relevant slice.
        """

        refsub = ref.copy()
        # If the reference file data might have a dimension greater than
        # two, use this syntax:
        # refsub.data = ref.data[..., slc[0], slc[1]].copy()
        refsub.data = ref.data[slc[0], slc[1]].copy()
        if hasattr(ref, "dq"):
            refsub.dq = ref.dq[slc[0], slc[1]].copy()

        return refsub

    def get_parameters(self):
        """Read capture and decay parameters from a reference table.

        Returns
        -------
        par0 : ndarray
            Column "capture0" from the trappars table.

        par1 : ndarray
            Column "capture1" from the trappars table.

        par2 : ndarray
            Column "capture2" from the trappars table.

        par3 : ndarray
            Column "decay_param" from the trappars table.
        """

        data = self.trappars_model.trappars_table
        par0 = data["capture0"].copy()
        par1 = data["capture1"].copy()
        par2 = data["capture2"].copy()
        par3 = data["decay_param"].copy()

        return par0, par1, par2, par3

    def compute_slope(self, integ):
        """Compute an estimate of the slope of the ramp for each pixel.

        Extended Summary
        ----------------
        We need a value for the slope that will not include cosmic-ray
        jumps, and groups that are flagged as saturated must also not
        contribute to the slope.  The approach will be to find the
        difference between adjacent groups, set the difference to a very
        large value if the group was saturated, then sort the differences
        along the axis for the ramp.  All the saturated groups will then
        be at the high end.  Cosmic-ray-affected groups should be just
        below the saturated groups (unless the effect of a CR was very
        small, in which case it won't matter where it is).  We can therefore
        ignore all the CR jumps and saturated values just by knowing how
        many of them there are.  The average of the remaining differences
        is the slope.

        Two arrays are returned, both giving the slope at each pixel, but
        with different units for the slope.  The reason for returning both
        arrays is that there is not a single factor relating the two; the
        factor can differ from one pixel to another.

        Parameters
        ----------
        integ : int
            The number (index) of the current integration.

        Returns
        -------
        grp_slope : ndarray, 2-D
            The ramp slope in units of counts (DN) per group.

        slope : ndarray, 2-D
            The ramp slope in units of fraction of the persistence
            saturation limit per second.
        """

        (_, ngroups, ny, nx) = self.output_obj.shape
        if ngroups == 1:
            # This won't be accurate, because there's only one group.
            grp_slope = self.output_obj.data[integ, 0, :, :]
            if hasattr(self.persistencesat, "dq"):
                mask = (np.bitwise_and(self.persistencesat.dq,
                                       dqflags.pixel["DO_NOT_USE"]) > 0)
                if mask.sum() == 0:
                    mask = None
            else:
                mask = None
            if mask is None:
                slope = grp_slope / (self.persistencesat.data * self.tgroup)
            else:
                # Set to 1 so we don't divide by 0.
                persat = np.where(mask, 1., self.persistencesat.data)
                slope = np.where(mask, 0.,
                                 grp_slope / (persat * self.tgroup))
            return grp_slope, slope

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

        # grp_slope has units of DN / group
        grp_slope = np.sum(sdiff, axis=0) / upper.astype(np.float32)
        grp_slope[bad] = 0.
        del bad

        # slope will have units (DN / persistence_saturation_limit) / second,
        # where persistence_saturation_limit is in units of DN.
        if hasattr(self.persistencesat, "dq"):
            mask = (np.bitwise_and(self.persistencesat.dq,
                                   dqflags.pixel["DO_NOT_USE"]) > 0)
            if mask.sum() == 0:
                mask = None
        else:
            mask = None
        if mask is None:
            slope = grp_slope / (self.persistencesat.data * self.tgroup)
        else:
            # Set to 1 so we don't divide by 0.
            persat = np.where(mask, 1., self.persistencesat.data)
            slope = np.where(mask, 0.,
                             grp_slope / (persat * self.tgroup))

        return grp_slope, slope

    def get_capture_param(self, par, k):
        """Extract capture parameters for the current trap family.

        Parameters
        ----------
        par : tuple of ndarray
            These were read from the trap parameters reference table.
            Each element of the tuple is a column from the table.  Each
            row of the table is for a different trap family.

        k : int
            Index of the current trap family

        Returns
        -------
        tuple of three floats
            These are the capture parameters for the current trap family.
        """

        (par0, par1, par2) = par[0:3]

        return par0[k], par1[k], par2[k]

    def get_decay_param(self, par, k):
        """Extract decay parameter(s) for the current trap family.

        Parameters
        ----------
        par : tuple of ndarray
            These were read from the trap parameters reference table.
            Each element of the tuple is a column from the table.  Each
            row of the table is for a different trap family.

        k : int
            Index of the current trap family

        Returns
        -------
        float
            This is the decay parameter for the current trap family.
        """

        par3 = par[3]

        return par3[k]

    def get_group_info(self, integ):
        """Get some metadata.

        Extended Summary
        ----------------
        This method populates these attributes:
            self.tframe
            self.tgroup
            self.ngroups
            self.nframes
            self.groupgap
            self.nresets

        Parameters
        ----------
        integ : int
            Integration number.  This is needed because the number of
            resets before the first integration can be different from the
            number of resets between integrations.
        """

        shape = self.output_obj.data.shape

        self.tframe = self.output_obj.meta.exposure.frame_time
        self.tgroup = self.output_obj.meta.exposure.group_time
        ngroups = self.output_obj.meta.exposure.ngroups
        if ngroups != shape[-3]:
            log.warning("model.meta and data disagree about ngroups:")
            log.warning("  %d vs %d; will use %d.",
                        ngroups, shape[-3], shape[-3])
        self.ngroups = shape[-3]
        self.nframes = self.output_obj.meta.exposure.nframes
        self.groupgap = self.output_obj.meta.exposure.groupgap
        try:
            if integ == 0:
                self.nresets = self.output_obj.meta.exposure.nresets_at_start
            else:
                self.nresets = self.output_obj.meta.exposure.nresets_between_ints
        except AttributeError:
            if self.output_obj.meta.instrument.detector == "MIRI":
                self.nresets = 0
            else:
                self.nresets = 1

    def predict_capture(self, capture_param_k, trap_density, integ,
                        grp_slope, slope):
        """Compute the number of traps that will be filled in time dt.

        This is based on Michael Regan's trapcapturemodel.pro.

        Parameters
        ----------
        capture_param_k : tuple of three floats
            Three values read from a reference table.  These will be from
            three separate columns but just one row; the row corresponds
            to the current trap family.  (The _k in the variable name
            indicates that the values are for one trap family.)

        trap_density : ndarray, 2-D
            Image of the total number of traps per pixel.

        integ : int
            Integration number.

        grp_slope : ndarray, 2-D
            The slope of the ramp at each pixel, in units of counts (DN)
            per group.  See also `slope`.
            The slope was computed from the pixel values that were not
            saturated and were not affected by jumps, based on flags in
            the groupdq extension.

        slope : ndarray, 2-D
            The slope of the ramp at each pixel, in units of
            fraction of the persistence saturation limit per second.
            This is the same as `grp_slope` except for units.

        Returns
        -------
        ndarray, 2-D
            The computed traps_filled at the end of the integration.
        """

        data = self.output_obj.data[integ, :, :, :]

        t_frame = self.tframe
        t_group = self.tgroup
        ngroups = self.ngroups
        nresets = self.nresets

        # nresets (usually equal to 1) was included because, in Mike
        # Regan's words:  "you get an extra frame of soak due to the
        # full frame reset at the beginning of the integration."
        totaltime = ngroups * t_group + nresets * t_frame

        # Find pixels exceeding the persistence saturation limit (full well).
        pflag = (data > self.persistencesat.data)
        if hasattr(self.persistencesat, "dq"):
            mask = (np.bitwise_and(self.persistencesat.dq,
                                   dqflags.pixel["DO_NOT_USE"]) > 0)
            mshape = mask.shape
            mask = mask.reshape((1,) + mshape)
            m3 = np.repeat(mask, ngroups, 0)
            pflag[m3] = False
            del mask, m3

        # All of these are 2-D arrays.
        sat_count = pflag.sum(axis=0, dtype=np.intp)
        del pflag
        sattime = sat_count.astype(np.float64) * t_group
        dt = totaltime - sattime

        # Traps that were filled due to the linear portion of the ramp.
        ramp_traps_filled = self.predict_ramp_capture(
            capture_param_k,
            trap_density, slope, dt)

        filled = ramp_traps_filled.copy()
        mask = (sat_count > 0)
        any_saturated = np.any(mask)
        if any_saturated:
            # Traps that were filled due to the saturated portion of the ramp.
            filled[mask] = self.predict_saturation_capture(
                capture_param_k,
                trap_density[mask],
                ramp_traps_filled[mask],
                sattime[mask], sat_count[mask], ngroups)
        del sat_count, sattime, ramp_traps_filled, mask

        # Traps that were filled due to cosmic-ray jumps.
        cr_filled = self.delta_fcn_capture(
            capture_param_k,
            trap_density, integ,
            grp_slope, ngroups, t_group)
        filled += cr_filled

        return filled

    def predict_ramp_capture(self, capture_param_k, trap_density, slope, dt):
        """Compute the number of traps that will be filled in time dt.

        This is based on Michael Regan's predictrampcapture3.pro.

        Parameters
        ----------
        capture_param_k : tuple
            Three values read from a reference table.  These will be from
            three separate columns but just one row; the row corresponds
            to the current trap family.  (The _k in the variable name
            indicates that the values are for one trap family.)

        trap_density : ndarray, 2-D
            Image of the total number of traps per pixel.

        slope : ndarray, 2-D
            Array of the slope of the ramp at each pixel.  The slope was
            computed from the pixel values that were not saturated and were
            not affected by jumps, based on flags in the groupdq extension.
            The unit is fraction of the persistence saturation limit
            per second.

        dt : float
            The time interval (unit = second) over which the charge capture
            is to be computed.  This does not include saturated groups.

        Returns
        -------
        ndarray, 2-D
            The computed traps_filled at the end of the integration.
        """

        (par0, par1, par2) = capture_param_k
        if par1 == 0:
            log.error("Capture parameter is zero; parameters are %g, %g, %g",
                      par0, par1, par2)
            tau = 1.e10                 # arbitrary "big" number
        else:
            tau = 1. / abs(par1)

        traps_filled = (trap_density * slope**2
                        * (dt**2 * (par0 + par2) / 2.
                           + par0 * (dt * tau + tau**2) * np.exp(-dt / tau)
                           - par0 * tau**2))

        traps_filled *= SCALEFACTOR
        return traps_filled

    def predict_saturation_capture(self, capture_param_k, trap_density,
                                   incoming_filled_traps,
                                   sattime, sat_count, ngroups):
        """Compute number of traps filled due to saturated pixels.

        This is based on Michael Regan's predictsaturationcapture.pro.

        Extended Summary
        ----------------
        This should not be called for ramps that do not have any groups
        that exceed the persistence saturation limit.  One reason is that
        `incoming_filled_traps` can be so small that `exp_filled_traps`
        would be negative.

        `trap_density`, `incoming_filled_traps`, `sattime`, and `sat_count`
        were all 2-D arrays in the calling function `predict_capture`, but
        these arrays have been masked to select only ramps with at least
        one saturated group, so in this function these arrays are 1-D.

        Parameters
        ----------
        capture_param_k : tuple
            Three values read from a reference table.  These will be from
            three separate columns but just one row; the row corresponds
            to the current trap family.  (The _k in the variable name
            indicates that the values are for one trap family.)

        trap_density : ndarray
            Image of the total number of traps per pixel.

        incoming_filled_traps : ndarray
            Traps filled due to linear portion of the ramp.
            This may be modified in-place.

        sattime : ndarray
            Time (seconds) during which each pixel was saturated.

        sat_count : array, int
            For each pixel, the number of groups with value exceeding the
            persistence saturation limit.

        ngroups : int
            The number of groups in the ramp

        Returns
        -------
        ndarray, 2-D
            The computed traps_filled at the end of the integration.
        """

        (par0, par1, par2) = capture_param_k
        par1 = abs(par1)        # the minus sign will be specified explicitly

        # For each pixel that had no ramp before saturation, fill all the
        # instantaneous traps; otherwise, they were filled during the ramp.
        flag = (sat_count == ngroups)
        incoming_filled_traps[flag] = trap_density[flag] * par2

        # Find out how many exponential traps have already been filled.
        exp_filled_traps = incoming_filled_traps - trap_density * par2

        # Subtract the already filled traps from the total traps possible
        # to fill.
        empty_traps = trap_density * par0 - exp_filled_traps

        # Filling of the empty traps depends only on the exponential
        # component; the instantaneous component would have been filled
        # during the ramp or above.
        new_filled_traps = empty_traps * (1. - np.exp(-par1 * sattime))
        total_filled_traps = incoming_filled_traps + new_filled_traps

        return total_filled_traps

    def delta_fcn_capture(self, capture_param_k, trap_density, integ,
                          grp_slope, ngroups, t_group):
        """Compute number of traps filled due to cosmic-ray jumps.

        Extended Summary
        ----------------
        If a cosmic-ray hit was in group number g (meaning that
        `data[integ, g, y, x]` was higher than expected), then
        delta_t = (ngroups - g - 0.5) * t_group
        is the time from the CR-affected group to the end of the
        integration, assuming that the CR hit in the middle (timewise)
        of the group.
        cr_filled = trap_density * jump
                    * (par0 * (1 - exp(-delta_t / tau)) + par2)

        Parameters
        ----------
        capture_param_k : tuple
            Three values read from a reference table.  These will be from
            three separate columns but just one row; the row corresponds
            to the current trap family.  (The _k in the variable name
            indicates that the values are for one trap family.)

        trap_density : ndarray, 2-D
            Image of the total number of traps per pixel.

        integ : int
            Integration number.

        grp_slope : ndarray, 2-D
            Array of the slope of the ramp at each pixel.  The slope was
            computed from the pixel values that were not saturated and were
            not affected by jumps, based on flags in the groupdq extension.
            The unit is counts (DN) per group.

        ngroups : int
            Total number of groups in the integration.

        t_group : float
            The time (seconds) from the start of one group to the start
            of the next group.

        Returns
        -------
        ndarray, 2-D
            The computed cr_filled at the end of the integration.
        """

        (par0, par1, par2) = capture_param_k
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
                jump = ((data[z, cr_flag[0], cr_flag[1]]
                         - data[z_prev, cr_flag[0], cr_flag[1]])
                        - grp_slope[cr_flag])
                jump = np.where(jump < 0., 0., jump)
                cr_filled[cr_flag] += trap_density[cr_flag] * jump \
                    * (par0 * (1. - math.exp(par1 * delta_t))
                       + par2)

        cr_filled *= SCALEFACTOR
        return cr_filled

    def compute_decay(self, traps_filled, decay_param, delta_t):
        """Compute the number of trap decays.

        This is based on Michael Regan's trapdecay.pro.

        Parameters
        ----------
        traps_filled : ndarray, 2-D
            This is an image of the number of filled traps in each pixel
            for the current trap family.

        decay_param : float
            The decay parameter.  This is negative, but otherwise it's
            the reciprocal of the e-folding time for trap decay for the
            current trap family.

        delta_t : float
            The time interval (unit = second) over which the trap decay
            is to be computed.

        Returns
        -------
        decayed : ndarray, 2-D
            Image of the computed number of trap decays for each pixel,
            for the current trap family.
        """

        if decay_param == 0.:
            decayed = traps_filled * 0.
        else:
            tau = 1. / abs(decay_param)
            decayed = traps_filled * (1. - math.exp(-delta_t / tau))

        return decayed
