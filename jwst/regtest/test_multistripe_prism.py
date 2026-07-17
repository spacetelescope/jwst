import os

import numpy as np
import pytest
from astropy.io import fits

from jwst import datamodels
from jwst.lib import reffile_utils
from jwst.stpipe import Step

parent_root = "jw01118001001_03101_00001-seg001"

params = {}
params["SUB256M2_PRISM"] = {
    "NREADS1": 4,
    "NSKIPS1": 1020,
    "NREADS2": 252,
    "NSKIPS2": 0,
    "NSTRIPES": 2,
}
params["SUB128M4_PRISM"] = {
    "NREADS1": 4,
    "NSKIPS1": 1020,
    "NREADS2": 124,
    "NSKIPS2": 0,
    "NSTRIPES": 4,
}
params["SUB64M8_PRISM"] = {
    "NREADS1": 4,
    "NSKIPS1": 1020,
    "NREADS2": 60,
    "NSKIPS2": 0,
    "NSTRIPES": 8,
}
params["SUB32M16_PRISM"] = {
    "NREADS1": 4,
    "NSKIPS1": 1020,
    "NREADS2": 28,
    "NSKIPS2": 0,
    "NSTRIPES": 16,
}


def create_multistripe(base_datamodel, format="SUB256M2_PRISM"):
    """
    Create a multistripe datamodel for testing.

    Parameters
    ----------
    base_datamodel : datamodels.DataModel
        A base datamodel to modify.
    format : str
        The format of the multistripe datamodel to create.

    Returns
    -------
    multistripe_datamodel : datamodels.DataModel
        A new datamodel with the multistripe format.
    """
    # Create a copy of the base datamodel
    multistripe_datamodel = base_datamodel.copy()
    data_array_science = multistripe_datamodel.data
    fastaxis = multistripe_datamodel.meta.subarray.fastaxis
    slowaxis = multistripe_datamodel.meta.subarray.slowaxis
    nskips1 = params[format]["NSKIPS1"]
    nreads1 = params[format]["NREADS1"]
    nskips2 = params[format]["NSKIPS2"]
    nreads2 = params[format]["NREADS2"]
    nstripes = params[format]["NSTRIPES"]
    # Get the shape of the data
    nints, ngroups, nrows, ncols = data_array_science.shape
    # Output data array has same number of rows but only 1/NSTRIPES the columns in each integration
    # (Reversed in detector format)
    output_data = np.zeros(
        (nints * nstripes, ngroups, nrows, ncols // nstripes), dtype=data_array_science.dtype
    )
    # Transform to detector format
    science_detector = reffile_utils.science_detector_frame_transform(
        data_array_science, fastaxis, slowaxis
    )
    output_detector = reffile_utils.science_detector_frame_transform(
        output_data, fastaxis, slowaxis
    )
    # Loop over integrations
    output_integration = 0
    for i in range(nints):
        input_current_row = 0
        for stripe in range(nstripes):
            # Add the reference pixels, made up values
            output_detector[output_integration, :, 0:4] = output_integration * 0.1
            # Copy the science data
            rowstart = input_current_row
            rowstop = input_current_row + nreads2
            output_detector[output_integration, :, 4 : 4 + nreads2] = science_detector[
                i, :, rowstart:rowstop
            ]
            input_current_row = input_current_row + nreads2
            output_integration += 1
    # Update the datamodel with the new data, transformed to science frame
    multistripe_datamodel.data = reffile_utils.detector_science_frame_transform(
        output_detector, fastaxis, slowaxis
    )
    # Update the metadata to reflect the new shape and parameters
    multistripe_datamodel.meta.subarray.xstart = 1
    if multistripe_datamodel.meta.instrument.detector == "NRS2":
        multistripe_datamodel.meta.subarray.xstart = 2048
    multistripe_datamodel.meta.subarray.multistripe_skips1 = nskips1
    multistripe_datamodel.meta.subarray.multistripe_reads1 = nreads1
    multistripe_datamodel.meta.subarray.multistripe_skips2 = nskips2
    multistripe_datamodel.meta.subarray.multistripe_reads2 = nreads2
    multistripe_datamodel.meta.subarray.repeat_stripe = 1
    multistripe_datamodel.meta.subarray.interleave_reads1 = 1
    multistripe_datamodel.meta.subarray.superstripe_step = nreads2
    multistripe_datamodel.meta.subarray.num_superstripe = nstripes
    multistripe_datamodel.meta.subarray.name = format
    multistripe_datamodel.meta.exposure.nints = nints * nstripes
    multistripe_datamodel.meta.exposure.integration_end = nints * nstripes
    multistripe_datamodel.meta.subarray.xsize = ncols // nstripes
    # Adjust the integration times table
    itable = multistripe_datamodel.int_times
    colnames = itable.columns.names
    formats = itable.columns.formats
    newcoldefs = []
    newcollen = nints * nstripes
    newcoldefs.append(
        fits.Column(
            name="integration_number", format="J", array=np.arange(newcollen, dtype=np.int32)
        )
    )
    for column_number, column in enumerate(colnames):
        if column == "integration_number":
            continue
        oldcol = itable[column]
        newcol = np.zeros((newcollen), dtype=oldcol.dtype)
        diffs = np.diff(oldcol)
        deltas = diffs / nstripes
        # interpolate times between existing table entries
        for i in range(nints - 1):
            for j in range(nstripes):
                newcol[i * nstripes + j] = oldcol[i] + j * deltas[i]
        # Add the entries at the end
        for j in range(nstripes):
            newcol[(nints - 1) * nstripes + j] = oldcol[-1] + j * deltas[-1]
        newcoldefs.append(fits.Column(name=column, format=formats[column_number], array=newcol))
    newtable = fits.BinTableHDU.from_columns(newcoldefs)
    multistripe_datamodel.int_times = newtable.data
    return multistripe_datamodel


@pytest.mark.parametrize("detector", ["nrs1", "nrs2"])
@pytest.mark.parametrize("multistripe_name", list(params.keys()))
def test_run_pipelines(rtdata_module, detector, multistripe_name):
    """
    Test that the multistripe prism pipeline results match those from the parent data.

    Framed as a test rather than a fixture so that parameterization
    can be used.
    """
    rtdata = rtdata_module
    parent = f"{parent_root}_{detector}_uncal.fits"
    input_file = parent
    if not os.path.exists(parent):
        rtdata.get_data(f"nirspec/tso/{parent}")
        input_file = rtdata.input
    # The refpix step will be skipped for the parent data, as it uses the rows of real data as
    # proxy reference pixels.  Since the multistripe data MUST be run through the refpix step,
    # as this is where the data get reformatted from multistripe to normal format, the
    # MASK file is replaced with one where the REFERENCE_PIXELS bit is unset, and this
    # MASK file is used for both parent and multistripe data.
    mask = f"mask_norefpixels_{detector}.fits"
    maskfile = mask
    if not os.path.exists(mask):
        rtdata.get_data(f"nirspec/tso/{mask}")
        maskfile = rtdata.input
    # The SUPERBIAS file is used in both the saturation and superbias steps.  CRDS selection
    # is format-dependent, but the steps will also run with a full-frame superbias file, with
    # the code extracting the appropriate section.  To ensure the same reference is used for both
    # parent and multistripe data, the full-frame superbias file is used for both.
    superbias = f"superbias_full_{detector}.fits"
    superbiasfile = superbias
    if not os.path.exists(superbias):
        rtdata.get_data(f"nirspec/tso/{superbias}")
        superbiasfile = rtdata.input
    # Technically, only the data from after the refpix step and before need to be compared,
    # since the multistripe data are reformatted to normal format in the refpix step and
    # the stripe_utils routines are not used after the refpix step.
    # Since the refpix step is skipped for parent data, there is no saved data from after
    # the refpix step, so the comparison is done after the linearity step.
    # Results are saved after the saturation, superbias and linearity steps.
    # The comparison fails after the jump step, because occasionally there are events
    # in the parent data that aren't in the the multistripe data (since the reformatted
    # multistripe data are truncated in size relative to the parent data).  These events
    # can cause flagging of pixels quite far away from the event, so they spill over
    # into the region shared by the multistripe data, and cause the comparison to fail.
    # Similarly, the clean_flicker_noise step causes differences because the data sizes
    # are different, and the use of Fourier transforms causes non-local effects.
    # The saturation.n_pix_grow_sat parameter is set to 0 to avoid the growth of any
    # saturated pixels in the region that is only in the parent data into the common region.
    step_args = [
        "calwebb_detector1",
        input_file,
        f"--steps.dq_init.override_mask={maskfile}",
        f"--steps.saturation.override_superbias={superbiasfile}",
        "--steps.saturation.n_pix_grow_sat=0",
        f"--steps.superbias.override_superbias={superbiasfile}",
        "--steps.refpix.skip=True",
        "--steps.clean_flicker_noise.skip=True",
        "--steps.linearity.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.superbias.save_results=True",
    ]
    Step.from_cmdline(step_args)
    # Now run the multistripe data through the detector1 pipeline
    parent_model = datamodels.open(input_file)
    multistripe_model = create_multistripe(parent_model, format=multistripe_name)
    multistripe_filename = f"{multistripe_name.lower()}_{detector}_uncal.fits"
    multistripe_model.to_fits(multistripe_filename, overwrite=True)
    step_args = [
        "calwebb_detector1",
        multistripe_filename,
        f"--steps.dq_init.override_mask=mask_norefpixels_{detector}.fits",
        f"--steps.saturation.override_superbias=superbias_full_{detector}.fits",
        "--steps.saturation.n_pix_grow_sat=0",
        f"--steps.superbias.override_superbias=superbias_full_{detector}.fits",
        "--steps.dark_current.skip=True",
        "--steps.clean_flicker_noise.skip=True",
        "--steps.linearity.save_results=True",
        "--steps.saturation.save_results=True",
    ]
    Step.from_cmdline(step_args)
    assert True


@pytest.mark.parametrize("detector", ["nrs1", "nrs2"])
@pytest.mark.parametrize("multistripe_name", list(params.keys()))
def test_uncal(detector, multistripe_name):
    """
    Test that the create_multistripe uncal files match the data in the parent uncal files.
    """
    parent = f"{parent_root}_{detector}_uncal.fits"
    with fits.open(parent) as f1:
        multistripe_filename = f"{multistripe_name.lower()}_{detector}_uncal.fits"
        sci1 = f1["sci"].data
        with fits.open(multistripe_filename) as f2:
            nstripes = f2[0].header["SSTR_NST"]
            columns_per_integration = 512 // nstripes - 4
            sci2 = f2["sci"].data
            for i in range(nstripes):
                if detector == "nrs1":
                    start_col = i * columns_per_integration
                    end_col = start_col + columns_per_integration
                    result = np.where(sci1[..., start_col:end_col] != sci2[i::nstripes, :, :, 4:])
                else:
                    end_col = 512 - i * columns_per_integration
                    start_col = end_col - columns_per_integration
                    result = np.where(sci1[..., start_col:end_col] != sci2[i::nstripes, :, :, :-4])
                assert len(result[0]) == 0, f"Data mismatch for {multistripe_name} at stripe {i}"


@pytest.mark.parametrize("detector", ["nrs1", "nrs2"])
@pytest.mark.parametrize("multistripe_name", list(params.keys()))
def test_saturation(detector, multistripe_name):
    """
    Test that the saturation step output is identical for the multistripe and parent data.

    Basically tests that the routines in lib/stripe_utils.py work correctly
    for the NIRSPEC prism multistripe data
    """
    parent = f"{parent_root}_{detector}_saturation.fits"
    everything_ok = True
    bad_extensions = []
    with fits.open(parent) as f1:
        multistripe_filename = f"{multistripe_name.lower()}_{detector}_saturation.fits"
        sci1 = f1["sci"].data
        with fits.open(multistripe_filename) as f2:
            nstripes = params[multistripe_name]["NSTRIPES"]
            columns_per_integration = 512 // nstripes - 4
            sci2 = f2["sci"].data
            for i in range(nstripes):
                if detector == "nrs1":
                    start_col = i * columns_per_integration
                    end_col = start_col + columns_per_integration
                    result_sci = np.where(
                        sci1[..., start_col:end_col] != sci2[i::nstripes, :, :, 4:]
                    )
                    result_pixeldq = np.where(
                        f1["pixeldq"].data[..., start_col:end_col]
                        != f2["pixeldq"].data[i::nstripes, :, 4:]
                    )
                    result_groupdq = np.where(
                        f1["groupdq"].data[..., start_col:end_col]
                        != f2["groupdq"].data[i::nstripes, :, :, 4:]
                    )
                else:
                    end_col = 512 - i * columns_per_integration
                    start_col = end_col - columns_per_integration
                    result_sci = np.where(
                        sci1[..., start_col:end_col] != sci2[i::nstripes, :, :, :-4]
                    )
                    result_pixeldq = np.where(
                        f1["pixeldq"].data[..., start_col:end_col]
                        != f2["pixeldq"].data[i::nstripes, :, :-4]
                    )
                    result_groupdq = np.where(
                        f1["groupdq"].data[..., start_col:end_col]
                        != f2["groupdq"].data[i::nstripes, :, :, :-4]
                    )
                if len(result_sci[0]) != 0:
                    everything_ok = False
                    bad_extensions.append("sci")
                if len(result_pixeldq[0]) != 0:
                    everything_ok = False
                    bad_extensions.append("pixeldq")
                if len(result_groupdq[0]) != 0:
                    everything_ok = False
                    bad_extensions.append("groupdq")

    assert everything_ok, f"Data mismatch found in extensions: {bad_extensions}"


@pytest.mark.parametrize("detector", ["nrs1", "nrs2"])
@pytest.mark.parametrize("multistripe_name", list(params.keys()))
def test_linearity(multistripe_name, detector):
    """
    Test that the linearity step output in the is identical for the multistripe and parent data.

    Basically tests that the routines in lib/stripe_utils.py work correctly
    for the NIRSPEC prism multistripe data.  After the refpix step, the multistripe
    data are in the same format as the parent data, so testing that the data are
    identical for the linearity step is the last test needed.
    """
    everything_ok = True
    bad_extensions = []
    parent = f"{parent_root}_{detector}_linearity.fits"
    with fits.open(parent) as f1:
        multistripe_filename = f"{multistripe_name.lower()}_{detector}_linearity.fits"
        sci1 = f1["sci"].data
        with fits.open(multistripe_filename) as f2:
            nstripes = params[multistripe_name]["NSTRIPES"]
            columns_per_integration = 512 // nstripes - 4
            total_columns = columns_per_integration * nstripes
            sci2 = f2["sci"].data
            if detector == "nrs1":
                start_col = 0
                end_col = start_col + total_columns
                result_sci = np.where(sci1[..., start_col:end_col] != sci2)
                result_pixeldq = np.where(
                    f1["pixeldq"].data[..., start_col:end_col] != f2["pixeldq"].data
                )
                result_groupdq = np.where(
                    f1["groupdq"].data[..., start_col:end_col] != f2["groupdq"].data
                )
            else:
                end_col = 512
                start_col = end_col - total_columns
                result_sci = np.where(sci1[..., start_col:end_col] != sci2)
                result_pixeldq = np.where(
                    f1["pixeldq"].data[..., start_col:end_col] != f2["pixeldq"].data
                )
                result_groupdq = np.where(
                    f1["groupdq"].data[..., start_col:end_col] != f2["groupdq"].data
                )
            if len(result_sci[0]) != 0:
                everything_ok = False
                bad_extensions.append("sci")
            if len(result_pixeldq[0]) != 0:
                everything_ok = False
                bad_extensions.append("pixeldq")
            if len(result_groupdq[0]) != 0:
                everything_ok = False
                bad_extensions.append("groupdq")

    assert everything_ok, f"Data mismatch found in extensions: {bad_extensions}"
