from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
#
#  Module for correcting for persistence

import numpy as np
import logging
from .. import datamodels

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
    def __init__(self, output_obj, traps_filled_model, trap_density_model,
                 traps_model, gain_model, sat_model):
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
            Table giving parameters for traps

        gain_model: image model
            Gain reference file

        sat_model: image model
            Saturation reference file
        """

        self.output_obj = output_obj
        self.traps_filled = traps_filled_model
        self.trap_density = trap_density_model
        self.traps_model = traps_model
        self.gain = gain_model
        self.saturation = sat_model

        # Initial values.
        self.save_section = None


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

        # Note that this might be a subarray.
        persistence = self.output_obj.data[0, 0, :, :] * 0.
        (nints, ngroups, ny, nx) = shape
        # xxx test debug t_group = self.output_obj.meta.exposure.group_time
        t_group = 1.
        # Time of one integration, in seconds.  This is the duration,
        # not necessarily the effective integration time.
        t_int = t_group * float(ngroups)

        if self.traps_filled is None:
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
        # Set meta.filename to the name of the input file, so that the
        # save_model() method of PersistenceStep will use that name (but
        # with a different suffix) as the output traps_filled file name,
        # to (hopefully clearly) link them.
        self.traps_filled.meta.filename = self.output_obj.meta.filename

        self.traps_filled = no_NaN(self.traps_filled, 0., zap_nan=True)
        self.gain = no_NaN(self.gain, 1., zap_nan=True, zap_zero=True)
        self.saturation = no_NaN(self.saturation, 1.e7, zap_nan=True)
        log.info("xxx do_all:  self.traps_filled.mean = %g",
                 self.traps_filled.data[1].mean())
        log.info("xxx do_all:  self.gain.mean = %g", self.gain.data[1].mean())
        log.info("xxx do_all:  self.saturation.mean = %g",
                 self.saturation.data[1].mean())

        log.info("xxx do_all:  self.trap_density.data.mean = %g",
                 self.trap_density.data.mean())

        # Read the table of capture and decay parameters.
        par = self.get_parameters()
        nfamilies = len(par[0])

        if have_traps_filled:   # was there an actual traps_filled file?
            log.info("xxx do_all:  yes, we have a traps_filled file")
            # Update traps_filled to the start of the exposure.  to_start
            # is the time difference in seconds from the traps_filled file
            # to the current exposure.
            to_start = (self.output_obj.meta.exposure.start_time -
                        self.traps_filled.meta.exposure.end_time) * 86400.
            # xxx test debug Convert to units of group time.
            to_start /= self.output_obj.meta.exposure.group_time  # xxx test debug
            log.info("xxx do_all:  to_start = %g 'groups'", to_start)
            for k in range(nfamilies):
                decay = self.compute_decay(self.traps_filled.data[k],
                                           par[3][k], to_start)
                log.info("xxx do_all:  k = %d, decay.mean = %g",
                         k, decay.mean())
                self.traps_filled.data[k, :, :] -= decay
                log.info("xxx do_all:  traps_filled.data.mean = %g",
                         self.traps_filled.data[k, :, :].mean())

        # If the input is a subarray, extract matching sections of
        # reference images.
        if self.ref_matches_sci(self.traps_filled, self.output_obj):
            is_subarray = False
            save_slice = None
            tr_filled = self.traps_filled.copy()
        else:
            # save_slice is the location of the subarray.
            (tr_filled, save_slice) = self.get_subarray(self.traps_filled,
                                                        self.output_obj)
            is_subarray = True
        log.info("xxx is_subarray = %s", str(is_subarray))
        # self.traps_filled and tr_filled are 3-D arrays, one image plane
        # for each trap family.
        # self.traps_filled is the size of the full detector.
        # tr_filled is the size of the science image, so it may be either
        # the full detector or a subarray.

        if not self.ref_matches_sci(self.trap_density, self.output_obj):
            (self.trap_density, _) = self.get_subarray(self.trap_density,
                                                       self.output_obj)

        if not self.ref_matches_sci(self.gain, self.output_obj):
            (self.gain, _) = self.get_subarray(self.gain, self.output_obj)

        if not self.ref_matches_sci(self.saturation, self.output_obj):
            (self.saturation, _) = self.get_subarray(self.saturation,
                                                     self.output_obj)

        # xxx These are for saving info in output files.  They're
        # 3-D arrays, one plane for each trap family.  If the input is a
        # subarray, so are these.
        temp_filled = tr_filled.data * 0.       # xxx test debug
        temp_decayed = temp_filled.copy()       # xxx test debug

        # self.traps_filled will be decreased with each group, to account
        # for the decay of traps.
        # tr_filled will be increased each group, to account for the
        # capture of charge by traps.
        # After processing all integrations, tr_filled will be added to
        # self.traps_filled, at the appropriate location if the science
        # image is a subarray.
        for integ in range(nints):
            slope = self.compute_slope(integ)
            temp_filled[:, :, :] = 0.           # xxx
            temp_decayed[:, :, :] = 0.          # xxx
            for group in range(ngroups):
                persistence[:, :] = 0.          # initialize
                for k in range(nfamilies):
                    capture_param_k = self.get_capture_param(par, k)
                    decay_param_k = self.get_decay_param(par, k)
                    # This may be a subarray.
                    filled = self.compute_capture(capture_param_k,
                                                  self.trap_density.data,
                                                  slope, t_group)
                    log.info("yyy k = %d:  %g, %g, %g, %g, %g, %g --> %g",
                             k, capture_param_k[0], capture_param_k[1],
                             capture_param_k[2],
                             self.trap_density.data[900, 900],
                             slope[900, 900], t_group,
                             filled[900, 900])
                    # This will be the size of the full detector.
                    decayed = self.compute_decay(self.traps_filled.data[k],
                                                 decay_param_k, t_group)
                    self.traps_filled.data[k, :, :] -= decayed
                    if is_subarray:
                        persistence += decayed[save_slice[0], save_slice[1]]
                    else:
                        persistence += decayed
                    # tr_filled accumulates throughout the entire exposure.
                    tr_filled.data[k, :, :] += filled
                    # These accumulate throughout each integration.
                    temp_filled[k, :, :] += filled              # xxx
                    if is_subarray:                             # xxx
                        temp_decayed[k, :, :] += decayed[save_slice[0],
                                                         save_slice[1]]
                    else:
                        temp_decayed[k, :, :] += decayed        # xxx
                # Persistence was computed in electrons, so convert to DN.
                self.output_obj.data[integ, group, :, :] -= \
                        persistence / self.gain.data
            # save for testing and debugging                    # xxx
            fits.writeto("filled_{}.fits".format(integ),        # xxx
                         data=temp_filled, overwrite=True)      # xxx
            fits.writeto("decayed_{}.fits".format(integ),       # xxx
                         data=temp_decayed, overwrite=True)     # xxx

        # Update self.traps_filled to include traps that have captured
        # charge during the current exposure.
        if is_subarray:
            for k in range(nfamilies):
                self.traps_filled.data[k, save_slice[0], save_slice[1]] += \
                        tr_filled.data[k, :, :].copy()
        else:
            for k in range(nfamilies):
                self.traps_filled.data[k, :, :] += \
                        tr_filled.data[k, :, :].copy()

        return (self.output_obj, self.traps_filled, skipped)


    def ref_matches_sci(self, ref, sci):
        """Test whether ref and sci are the same subarray."""

        # For the time being, just compare shapes.      xxx not finished
        if ref.shape[-2:] == sci.shape[-2:]:
            return True
        else:
            return False

        """
        xxx maybe use later
        # See if the reference and science model subarray parameters match
        if (ref.meta.subarray.xstart == sci.meta.subarray.xstart and
            ref.meta.subarray.xsize == sci.meta.subarray.xsize and
            ref.meta.subarray.ystart == sci.meta.subarray.ystart and
            ref.meta.subarray.ysize == sci.meta.subarray.ysize):
            return True
        else:
            return False
        """


    def get_subarray(self, ref, sci):
        """Extract a subarray from a reference file.

        Parameters
        ----------
        ref: data model
            A reference image.

        sci: data model
            The science data.

        Returns
        -------
        tuple (sub, s)
            `sub` is the subarray extracted from `ref`.  `sub` will
            sometimes be 2-D and other times 3-D.
            `s` is a two-element tuple, the Y and X slices that were used
            to extract the subarray.
        """

        ref_shape = ref.shape
        ref_nx = ref_shape[-1]
        ref_ny = ref_shape[-2]
        sci_shape = sci.shape
        sci_nx = sci_shape[-1]
        sci_ny = sci_shape[-2]

        # These are limits of slices.
        sci_x1 = sci.meta.subarray.xstart
        sci_x2 = sci_x1 + sci_nx
        sci_y1 = sci.meta.subarray.ystart
        sci_y2 = sci_y1 + sci_ny
        log.debug("sci xstart=%d, xstop=%d, ystart=%d, ystop=%d",
                  (sci_x1, sci_x2 - 1, sci_y1, sci_y2 - 1))

        ref_x1 = ref.meta.subarray.xstart
        ref_y1 = ref.meta.subarray.ystart
        if ref_x1 is None:
            ref_x1 = sci_x1
        if ref_y1 is None:
            ref_y1 = sci_y1
        ref_x2 = ref_x1 + ref_nx
        ref_y2 = ref_y1 + ref_ny
        log.debug("ref xstart=%d, xstop=%d, ystart=%d, ystop=%d",
                  (ref_x1, ref_x2 - 1, ref_y1, ref_y2 - 1))

        # Compute the slicing indexes
        xstart = sci_x1 - ref_x1
        ystart = sci_y1 - ref_y1
        xstop = xstart + sci_nx
        ystop = ystart + sci_ny
        log.debug("ref slice %d:%d, %d:%d", ystart, ystop, xstart, xstop)

        # Check for errors in the slice indexes
        if (xstart < 0 or ystart < 0 or
            xstop >= ref.data.shape[-1] or ystop >= ref.data.shape[-2]):
            log.error("Science and reference file arrays not compatible")
            raise ValueError("Can't extract matching subarray from "
                             "reference data")

        s = (slice(ystart, ystop), slice(xstart, xstop))

        # Slice the reference model arrays
        sub = ref.copy()
        sub.data = ref.data[..., s[0], s[1]].copy()
        if ref.__hasattr__(err):
            sub.err = ref.err[..., s[0], s[1]].copy()
        if ref.__hasattr__(dq):
            sub.dq = ref.dq[..., s[0], s[1]].copy()

        # Return the sliced reference model and the y and x slices
        return (sub, s)


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
        par0 = data.field("capture0").copy()
        par1 = data.field("capture1").copy()
        par2 = data.field("capture2").copy()
        par3 = data.field("decay_param").copy()

        return (par0, par1, par2, par3)


    def compute_slope(self, integ):
        """Compute the slope of the ramp at each pixel.

        Parameters
        ----------
        integ: int
            The number (index) of the current integration.

        Returns
        -------
        2-D ndarray
            The array of ramp slopes; the unit is fraction of the
            persistence saturation limit per group.
            xxx Are these the right units?
        """

        # This is based on Michael Regan's trapcapturemodel.pro.

        # xxx This is a workaround, because we don't have such a reference
        # file.  Should we create one?
        persistencesaturationlimit = 90000.

        (_, ngroups, ny, nx) = self.output_obj.shape
        index = (self.output_obj.data[integ, :, :, :] > self.saturation.data)
        satcount = index.sum(axis=0, dtype=np.int32)
        goodindex = ngroups - satcount - 1
        del index
        # If goodindex is -1, the code below will index out of bounds, and
        # if goodindex is 0, we'll divide by zero; so modify it if it's < 1.
        mask = np.where(goodindex < 1)
        if len(mask[0]) > 0:
            log.warning("%d ramps were saturated in the first or "
                        "second group.", len(mask[0]))
            goodindex[mask] = 1
        del mask

        # We need the science data in units of electrons.  (no, maybe not)
        # xxx data = self.output_obj.data * self.gain.data
        data = self.output_obj.data

        (iy, ix) = np.indices((ny, nx), dtype=np.int32)
        accumulatedEl = data[integ, goodindex, iy, ix].astype(np.float64) - \
                        data[integ, 0, :, :].astype(np.float64)
        # persistencesaturationlimit should be a 2-D reference file.
        slope = (accumulatedEl / persistencesaturationlimit) / \
                goodindex.astype(np.float64)
        log.info("xxx compute_slope:  slope.mean = %g",
                 slope.mean())

        return slope


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


    def compute_capture(self, capture_param_k, trap_density, slope, t_int):
        """Compute the number of traps that will be filled in time t_int.

        This is based on Michael Regan's predictrampcapture3.pro.

        Parameters
        ----------
        capture_param_k: tuple
            Three values read from a reference table.  These will be
            from three separate columns but one row; the row corresponds
            to the current trap family.  (The _k in the variable name
            refers to the index of the trap family.)

        trap_density: 2-D ndarray
            Image of the total number of traps per pixel.

        slope: 2-D ndarray
            Image of the slope of the ramp per pixel.  The slope is
            computed from the non-saturated pixel values.  This value is
            not precise; it's based on the first and last (non-saturated)
            values, which includes cosmic-ray hits but treats them as if
            they were not discontinuous jumps but were instead spread
            uniformly along the ramp.

        t_int: float
            The time interval (unit = group count) over which the charge
            capture is to be computed.  (t_int implies integration time,
            but that is not necessarily the actual time interval.)

        Returns
        -------
        2-D ndarray
            The computed traps_filled at the end of the exposure.
        """

        (par0, par1, par2) = capture_param_k

        t1 = t_int**2 / 2.
        t2 = (-t_int / par1) * np.exp(par1 * t_int)
        t3 = (1.0 / par1**2) * np.exp(par1 * t_int)
        t4 = -1.0 / par1**2
        traps_filled = (trap_density * slope**2 *
                        (par0 * (t1 + t2 + t3 + t4) + t_int**2 * par2 / 2.))

        return traps_filled


    def compute_decay(self, traps_filled, decay_param, t_int):
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

        t_int: float
            The time interval (unit = group count) over which the trap
            decay is to be computed.

        Returns
        -------
        decay: 2-D ndarray
            Image of the computed number of trap decays for each pixel,
            for the current trap family.
        """

        log.info("xxx compute_decay:  traps_filled.mean = %g",
                 traps_filled.mean())
        return traps_filled * (1. - np.exp(decay_param * t_int))
