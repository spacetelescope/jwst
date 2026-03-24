import logging

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.lib import pipe_utils
from jwst.refpix import irs2_subtract_reference, reference_pixels
from jwst.stpipe import Step

__all__ = ["RefPixStep"]

log = logging.getLogger(__name__)


class RefPixStep(Step):
    """Use reference pixels to correct bias drifts."""

    class_alias = "refpix"

    spec = """
        odd_even_columns = boolean(default=True) # Compute reference signal separately for even/odd columns
        use_side_ref_pixels = boolean(default=True) # Use side reference pixels for reference signal for each row
        side_smoothing_length = integer(default=11) # Median window smoothing height for side reference signal
        side_gain = float(default=1.0) # Multiplicative factor for side reference signal before subtracting from rows
        odd_even_rows = boolean(default=True) # Compute reference signal separately for even- and odd-numbered rows
        ovr_corr_mitigation_ftr = float(default=3.0) # Factor to avoid overcorrection of bad reference pixels for IRS2
        preserve_irs2_refpix = boolean(default=False) # Preserve reference pixels in output
        irs2_mean_subtraction = boolean(default=False) # Apply a mean offset subtraction before IRS2 correction
        refpix_algorithm = option("median", "sirs", default="median") # NIR full-frame side pixel algorithm
        sigreject = float(default=4.0) # Number of sigmas to reject as outliers
        gaussmooth = float(default=1.0) # Width of Gaussian smoothing kernel to use as a low-pass filter
        halfwidth = integer(default=30) # Half-width of convolution kernel to build
    """  # noqa: E501

    reference_file_types = ["refpix", "sirskernel"]

    def process(self, step_input):
        """
        Execute the reference pixel correction step.

        Parameters
        ----------
        step_input : str or `~stdatamodels.jwst.datamodels.RampModel`
            Input datamodel or file name.

        Returns
        -------
        result : ~stdatamodels.jwst.datamodels.RampModel`
            Result of applying the reference pixel correction step.
        """
        conv_kernel_params = {
            "refpix_algorithm": self.refpix_algorithm,
            "sirs_kernel_model": None,
            "sigreject": self.sigreject,
            "gaussmooth": self.gaussmooth,
            "halfwidth": self.halfwidth,
        }

        # Open the input data model
        result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)

        if pipe_utils.is_irs2(result):
            # Flag bad reference pixels first
            irs2_subtract_reference.flag_bad_refpix(
                result, n_sigma=self.ovr_corr_mitigation_ftr, flag_only=True
            )

            # If desired, do the normal refpix correction before IRS2, without
            # side pixel handling
            if self.irs2_mean_subtraction:
                if self.use_side_ref_pixels:
                    log.info("Turning off side pixel correction for IRS2")
                    self.use_side_ref_pixels = False
                reference_pixels.correct_model(
                    result,
                    self.odd_even_columns,
                    self.use_side_ref_pixels,
                    self.side_smoothing_length,
                    self.side_gain,
                    self.odd_even_rows,
                    conv_kernel_params,
                )

            # Now that values are updated, replace bad reference pixels
            irs2_subtract_reference.flag_bad_refpix(result, replace_only=True)

            # Get the necessary refpix reference file for IRS2 correction
            self.irs2_name = self.get_reference_file(result, "refpix")
            log.info(f"Using refpix reference file: {self.irs2_name}")

            # Check for a valid reference file
            if self.irs2_name == "N/A":
                log.warning("No refpix reference file found")
                log.warning("RefPix step will be skipped")
                result.meta.cal_step.refpix = "SKIPPED"
                return result

            # Load the reference file into a datamodel
            irs2_model = datamodels.IRS2Model(self.irs2_name)

            # Apply the IRS2 correction scheme
            result = irs2_subtract_reference.correct_model(
                result, irs2_model, preserve_refpix=self.preserve_irs2_refpix
            )

            if result.meta.cal_step.refpix != "SKIPPED":
                result.meta.cal_step.refpix = "COMPLETE"
            del irs2_model
            return result

        else:
            # Not an NRS IRS2 exposure. Do the normal refpix correction.

            # Get the reference file from CRDS or use user-supplied one
            # only for NIR full-frame data
            if self.refpix_algorithm == "sirs":
                if result.meta.instrument.name != "MIRI" and "FULL" in result.meta.subarray.name:
                    sirs_ref_filename = self.get_reference_file(result, "sirskernel")
                    if sirs_ref_filename == "N/A":
                        log.warning("No reference file found for the optimized convolution kernel.")
                        log.warning(
                            "REFPIX step will use the running median algorithm for side pixels."
                        )
                    else:
                        log.info(f"Using SIRS reference file: {sirs_ref_filename}")
                        sirs_kernel_model = datamodels.SIRSKernelModel(sirs_ref_filename)
                        conv_kernel_params["sirs_kernel_model"] = sirs_kernel_model
                elif result.meta.instrument.name == "MIRI":
                    log.info(
                        "Simple Improved Reference Subtraction (SIRS) not applied for MIRI data."
                    )
                elif "FULL" not in result.meta.subarray.name:
                    log.info(
                        "Simple Improved Reference Subtraction (SIRS) "
                        "not applied for subarray data."
                    )

            status = reference_pixels.correct_model(
                result,
                self.odd_even_columns,
                self.use_side_ref_pixels,
                self.side_smoothing_length,
                self.side_gain,
                self.odd_even_rows,
                conv_kernel_params,
            )

            if status == reference_pixels.REFPIX_OK:
                result.meta.cal_step.refpix = "COMPLETE"
            elif status == reference_pixels.SUBARRAY_DOESNTFIT:
                log.warning("Subarray doesn't fit in full-sized array")
                result.meta.cal_step.refpix = "SKIPPED"
            elif status == reference_pixels.BAD_REFERENCE_PIXELS:
                log.warning("No valid reference pixels, refpix step skipped")
                result.meta.cal_step.refpix = "SKIPPED"
            elif status == reference_pixels.SUBARRAY_SKIPPED:
                result.meta.cal_step.refpix = "SKIPPED"

            if (
                result.meta.subarray.num_superstripe is not None
                and result.meta.subarray.num_superstripe > 0
            ):
                result = collate_superstripes(result)

            return result


# TODO: reorganize stripe handling utilities into a separate module
def collate_superstripes(input_model):
    """
    Collate superstripes into arrays resembling the full detector/subarray shape.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel`
        The datamodel containing a reference-pixel corrected ramp
        for superstripe data.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The datamodel with superstripes collated into a single frame per set of stripes.
    """
    # First define the parent array shape
    fastaxis = np.abs(input_model.meta.subarray.fastaxis)

    slowsize = 2048
    slowdir = np.sign(input_model.meta.subarray.slowaxis)
    nreads1 = input_model.meta.subarray.multistripe_reads1

    # Generate slowaxis ranges to place stripes into parent frame
    stripe_ranges = pipe_utils.generate_superstripe_ranges(input_model)
    srlist = []
    for key in stripe_ranges:
        if slowdir < 0:
            # If slowaxis is negative, need to read the stripes in reverse order
            # (but not reverse column order within a stripe, as they've been
            # transformed to science frame)
            srlist.append([slowsize - x for x in stripe_ranges[key][0][::-1]])
        else:
            srlist.append(list(stripe_ranges[key][0]))

    # Determine integration/stripe numbers
    nints, ngroups, ny, nx = input_model.data.shape
    n_sstr = input_model.meta.subarray.num_superstripe
    nints_sci = nints // n_sstr

    # Initialize new array shapes
    if fastaxis == 1:
        newdata = np.full((nints_sci, ngroups, slowsize, nx), np.nan, dtype=">f4")
        newgdq = np.full((nints_sci, ngroups, slowsize, nx), 0, dtype="uint8")
        newpdq = np.full(
            (slowsize, nx), datamodels.dqflags.pixel["REFERENCE_PIXEL"], dtype="uint32"
        )
    else:
        newdata = np.full((nints_sci, ngroups, ny, slowsize), np.nan, dtype=">f4")
        newgdq = np.full((nints_sci, ngroups, ny, slowsize), 0, dtype="uint8")
        newpdq = np.full(
            (ny, slowsize), datamodels.dqflags.pixel["REFERENCE_PIXEL"], dtype="uint32"
        )

    # Work through each set of stripes, pushing them into a common frame
    for integ in range(nints_sci):
        for stripe in range(n_sstr):
            """
            Long term, there may be cases where stripes overlap.
            This could be refactored to store each stripe in a separate
            plane, and overlaps could be reduced using a function of choice,
            e.g. np.median.
            """
            # Determine fast orient
            if fastaxis == 1:
                newslice = np.s_[integ, :, srlist[stripe][0] : srlist[stripe][1], :]
                # Determine end of slowread to drop refpix from
                if slowdir < 0:
                    dataslice = np.s_[integ * n_sstr + stripe, :, :-nreads1, :]
                else:
                    dataslice = np.s_[integ * n_sstr + stripe, :, nreads1:, :]
                newdata[newslice] = input_model.data[dataslice]
                newgdq[newslice] = input_model.groupdq[dataslice]
                if integ == 0:
                    newpdq[newslice[-2], newslice[-1]] = input_model.pixeldq[
                        stripe, dataslice[-2], dataslice[-1]
                    ]
            else:
                newslice = np.s_[integ, :, :, srlist[stripe][0] : srlist[stripe][1]]
                if slowdir < 0:
                    dataslice = np.s_[integ * n_sstr + stripe, :, :, :-nreads1]
                else:
                    dataslice = np.s_[integ * n_sstr + stripe, :, :, nreads1:]
                newdata[newslice] = input_model.data[dataslice]
                newgdq[newslice] = input_model.groupdq[dataslice]
                if integ == 0:
                    newpdq[newslice[-2], newslice[-1]] = input_model.pixeldq[
                        stripe, dataslice[-2], dataslice[-1]
                    ]

    new_model = datamodels.RampModel(
        data=newdata,
        groupdq=newgdq,
        pixeldq=newpdq,
        int_times=input_model.int_times,
    )
    new_model.update(input_model)

    new_model = generate_stripe_int_times(new_model)
    new_model = clean_superstripe_metadata(new_model)

    return new_model


def clean_superstripe_metadata(input_model):
    """
    Update model metadata to match changes to arrays.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel`
        The model with updated data array shapes matching a parent
        frame, e.g. full frame or a subarray, which consist of
        multiple superstripe integrations.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The model cleaned of metadata indicating the presence
        of superstripe data.
    """
    input_model.meta.exposure.integration_start = np.ceil(
        input_model.meta.exposure.integration_start / input_model.meta.subarray.num_superstripe
    ).astype(int)
    input_model.meta.exposure.integration_end = np.ceil(
        input_model.meta.exposure.integration_end / input_model.meta.subarray.num_superstripe
    ).astype(int)
    input_model.meta.exposure.nints = np.ceil(
        input_model.meta.exposure.nints / input_model.meta.subarray.num_superstripe
    ).astype(int)

    input_model.meta.subarray.multistripe_reads1 = None
    input_model.meta.subarray.multistripe_reads2 = None
    input_model.meta.subarray.multistripe_skips1 = None
    input_model.meta.subarray.multistripe_skips2 = None
    input_model.meta.subarray.repeat_stripe = None
    input_model.meta.subarray.interleave_reads1 = None
    input_model.meta.subarray.num_superstripe = None
    input_model.meta.subarray.superstripe_step = None
    input_model.meta.subarray.ysize, input_model.meta.subarray.xsize = input_model.data.shape[-2:]

    return input_model


def generate_stripe_int_times(input_model):
    """
    Move input model INT_TIMES to stripe table, then condense ints to parent frame.

    Each output integration in the parent frame will have the start time from
    the first stripe readout and the end time from the last stripe readout included
    in the integration.  The mid time for the integration is assigned to the average
    of the start and end times.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel`
        The model with updated data array shapes matching a parent
        frame, e.g. full frame or a subarray, which consist of
        multiple superstripe integrations.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The model now with two INT_TIMES tables - one reflecting per-stripe
        information, preserved from input, and one reflecting times
        corresponding to the new parent frames.
    """
    if input_model.int_times is None or len(input_model.int_times) == 0:
        return input_model

    nstr = input_model.meta.subarray.num_superstripe
    nints_sci = len(input_model.int_times) // nstr

    otab = np.array(
        list(
            zip(
                np.repeat(np.arange(len(input_model.int_times) // nstr) + 1, nstr),
                np.arange(len(input_model.int_times)) % nstr + 1,
                [integ["int_start_MJD_UTC"] for integ in input_model.int_times],
                [integ["int_mid_MJD_UTC"] for integ in input_model.int_times],
                [integ["int_end_MJD_UTC"] for integ in input_model.int_times],
                [integ["int_start_BJD_TDB"] for integ in input_model.int_times],
                [integ["int_mid_BJD_TDB"] for integ in input_model.int_times],
                [integ["int_end_BJD_TDB"] for integ in input_model.int_times],
                strict=True,
            )
        ),
        dtype=input_model.get_dtype("int_times_stripe"),
    )

    input_model.int_times_stripe = otab

    cds_mjd_beg = np.full(nints_sci, np.nan, dtype="<f8")
    cds_mjd_mid = np.full(nints_sci, np.nan, dtype="<f8")
    cds_mjd_end = np.full(nints_sci, np.nan, dtype="<f8")
    cds_bjd_beg = np.full(nints_sci, np.nan, dtype="<f8")
    cds_bjd_mid = np.full(nints_sci, np.nan, dtype="<f8")
    cds_bjd_end = np.full(nints_sci, np.nan, dtype="<f8")

    input_inttimes = input_model.int_times
    nints_sci = len(input_inttimes) // nstr
    for i in range(nints_sci):
        cds_mjd_beg[i] = input_inttimes[i * nstr]["int_start_MJD_UTC"]
        cds_mjd_end[i] = input_inttimes[(i + 1) * nstr - 1]["int_end_MJD_UTC"]
        cds_mjd_mid[i] = (cds_mjd_end[i] - cds_mjd_beg[i]) / 2.0 + cds_mjd_beg[i]
        cds_bjd_beg[i] = input_inttimes[i * nstr]["int_start_BJD_TDB"]
        cds_bjd_end[i] = input_inttimes[(i + 1) * nstr - 1]["int_end_BJD_TDB"]
        cds_bjd_mid[i] = (cds_bjd_end[i] - cds_bjd_beg[i]) / 2.0 + cds_bjd_beg[i]

    otab2 = np.array(
        list(
            zip(
                np.arange(len(input_model.int_times) // nstr) + 1,
                cds_mjd_beg,
                cds_mjd_mid,
                cds_mjd_end,
                cds_bjd_beg,
                cds_bjd_mid,
                cds_bjd_end,
                strict=True,
            )
        ),
        dtype=input_model.int_times.dtype,
    )

    input_model.int_times = otab2

    return input_model
