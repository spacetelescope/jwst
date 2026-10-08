"""Mock NIRCam data in rate format."""

import numpy as np
from astropy.io import fits
from gwcs import wcs
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import nircam
from jwst.lib.stripe_utils import generate_substripe_ranges
from jwst.tests.wcs_helpers import get_reference_files

__all_ = ["nircam_image_model", "nircam_tsgrism_rate_model", "nircam_wfss_rate_model"]

# Default wcs information
DEFAULT_WCS_KW = {
    "wcsaxes": 2,
    "ra_ref": 53.1423683802,
    "dec_ref": -27.8171119969,
    "v2_ref": 86.103458,
    "v3_ref": -493.227512,
    "roll_ref": 45.04234459270135,
    "crpix1": 1024.5,
    "crpix2": 1024.5,
    "crval1": 53.1423683802,
    "crval2": -27.8171119969,
    "cdelt1": 1.74460027777777e-05,
    "cdelt2": 1.75306861111111e-05,
    "ctype1": "RA---TAN",
    "ctype2": "DEC--TAN",
    "pc1_1": -1,
    "pc1_2": 0,
    "pc2_1": 0,
    "pc2_2": 1,
    "cunit1": "deg",
    "cunit2": "deg",
}


# Example keyword values from jw01366002001_04103_00001
TSO_WCS_KW = {
    "wcsaxes": 2,
    "ra_ref": 217.3262354867486,
    "dec_ref": -3.444360390029802,
    "v2_ref": 73.106857,
    "v3_ref": -551.643196,
    "v3i_yang": 0.24336066,
    "vparity": -1,
    "roll_ref": 110.12344444842617,
    "xref_sci": 1581.0,
    "yref_sci": 35.0,
    "cdelt1": 1.76686111111111e-05,
    "cdelt2": 1.78527777777777e-05,
    "ctype1": "RA---TAN",
    "ctype2": "DEC--TAN",
    "pc1_1": -1,
    "pc1_2": 0,
    "pc2_1": 0,
    "pc2_2": 1,
    "cunit1": "deg",
    "cunit2": "deg",
}


# Ordered from lowest to highest subarray row (= order in which stripes are
# packed into the packed subarray by generate_stripe_reference).  The regions
# reference file is built with the same ordering, so the two stay in sync.
NRCA1_DHS_STRIPE_IDS = [10, 9, 8, 7]


def _nircam_rate_hdul(
    detector="NRCALONG",
    channel="LONG",
    module="A",
    filter_name="F444W",
    exptype="NRC_IMAGE",
    pupil="GRISMR",
    subarray="FULL",
    wcskeys=None,
):
    if wcskeys is None:
        wcskeys = DEFAULT_WCS_KW

    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header["TELESCOP"] = "JWST"
    phdu.header["FILENAME"] = "test+" + filter_name
    phdu.header["INSTRUME"] = "NIRCAM"
    phdu.header["CHANNEL"] = channel
    phdu.header["DETECTOR"] = detector
    phdu.header["FILTER"] = filter_name
    phdu.header["PUPIL"] = pupil
    phdu.header["MODULE"] = module
    phdu.header["TIME-OBS"] = "8:59:37"
    phdu.header["DATE-OBS"] = "2023-01-01"
    phdu.header["EXP_TYPE"] = exptype
    phdu.header["SUBARRAY"] = subarray
    scihdu = fits.ImageHDU()
    scihdu.header["EXTNAME"] = "SCI"
    scihdu.header.update(wcskeys)
    hdul.append(phdu)
    hdul.append(scihdu)
    return hdul


def nircam_image_rate_model(with_wcs=True):
    """
    Create a mock NIRCam image rate model.

    The data array is zero-filled with shape (10, 10).

    Parameters
    ----------
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The NIRCam image datamodel.
    """
    model = datamodels.ImageModel(_nircam_rate_hdul())
    model.data = np.zeros((10, 10))

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nircam.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def nircam_tsgrism_rate_model(filter_name="F322W2", with_wcs=True):
    """
    Create a mock NIRCam TSGRISM rateints model.

    The data array is zero-filled with shape (10, 10, 10).

    Parameters
    ----------
    filter_name : str, optional
        Filter name.
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.CubeModel`
        The NIRCam TSGRISM datamodel.
    """
    hdul = _nircam_rate_hdul(
        exptype="NRC_TSGRISM",
        pupil="GRISMR",
        filter_name=filter_name,
        detector="NRCALONG",
        subarray="SUBGRISM256",
        wcskeys=TSO_WCS_KW,
    )
    model = datamodels.CubeModel(hdul)
    model.data = np.zeros((10, 10, 10))
    model.meta.dither.x_offset = 0.0
    model.meta.dither.y_offset = 1.4485  # example from jw01366002001_04103_00001

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nircam.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def nircam_wfss_rate_model(pupil="GRISMR", with_wcs=True):
    """
    Create a mock NIRCam WFSS rateints model.

    The data array is zero-filled with shape (10, 10).

    Parameters
    ----------
    pupil : str, optional
        Pupil name.
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.CubeModel`
        The NIRCam WFSS datamodel.
    """
    hdul = _nircam_rate_hdul(exptype="NRC_WFSS", filter_name="F444W", pupil=pupil)
    model = datamodels.ImageModel(hdul)
    model.data = np.zeros((10, 10))

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nircam.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def _populate_dhs_shared_metadata(model, subarray="SUB260STRIPE4_DHS"):
    """
    Populate metadata that is identical between NRCA1 and NRCALONG DHS modes.

    Updates *model* in place.
    """
    # Aperture
    if subarray == "SUB164STRIPE4_DHS":
        model.meta.aperture.pps_name = "NRCA5_164STRIPE4_DHS_F322W2"
    else:
        model.meta.aperture.pps_name = "NRCA5_260STRIPE4_DHS_F322W2"

    # Exposure
    model.meta.exposure.type = "NRC_TSGRISM"
    model.meta.exposure.ngroups = 5

    # Observation
    model.meta.observation.date = "2026-05-03"
    model.meta.observation.time = "00:00:00.000"

    # Subarray
    model.meta.subarray.name = subarray
    model.meta.subarray.fastaxis = -1
    model.meta.subarray.slowaxis = 2
    model.meta.subarray.num_superstripe = 0
    model.meta.subarray.repeat_stripe = 1
    model.meta.subarray.xstart = 1
    model.meta.subarray.xsize = 2048
    model.meta.subarray.ystart = 1
    if subarray == "SUB164STRIPE4_DHS":
        model.meta.subarray.ysize = 164
    else:
        model.meta.subarray.ysize = 260

    # WCS info
    # Example values from jw04453025001_03103_00001-seg001_nrca1
    model.meta.wcsinfo.ra_ref = 265.74911600957546
    model.meta.wcsinfo.dec_ref = 66.93483954623517
    model.meta.wcsinfo.v2_ref = 120.311122
    model.meta.wcsinfo.v3_ref = -495.782084
    model.meta.wcsinfo.roll_ref = 245.14852360749214
    model.meta.wcsinfo.velosys = 617.4
    model.meta.wcsinfo.v3yangle = -0.42593583
    model.meta.wcsinfo.vparity = -1

    # Velocity aberration
    model.meta.velocity_aberration.scale_factor = 1.0

    # Offset from TA
    model.meta.dither.x_offset = 0.0
    model.meta.dither.y_offset = 1.39


def nircam_dhs_nrca1_sub260_rate_model(with_wcs=True):
    """
    Create a mock DHS NRCA1 rate CubeModel with subarray SUB260STRIPE4_DHS.

    Data array has shape (5, 260, 2048) and is filled with the stripe ID values.

    Parameters
    ----------
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.CubeModel`
        Mock NRCA1 DHS rate model.
    """
    model = datamodels.CubeModel((5, 260, 2048))
    _populate_dhs_shared_metadata(model)

    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.channel = "SHORT"
    model.meta.instrument.detector = "NRCA1"
    model.meta.instrument.filter = "F150W2"
    model.meta.instrument.pupil = "GDHS0"
    model.meta.instrument.module = "A"

    model.meta.subarray.interleave_reads1 = 1
    model.meta.subarray.multistripe_reads1 = 1
    model.meta.subarray.multistripe_reads2 = 64
    model.meta.subarray.multistripe_skips1 = 1514
    model.meta.subarray.multistripe_skips2 = 61

    model.meta.wcsinfo.siaf_xref_sci = 1024.5
    model.meta.wcsinfo.siaf_yref_sci = 24.5

    # Fill each packed stripe row range with the stripe's ID value so tests
    # can verify that the correct detector region ends up in each output slit.
    sub_ranges = generate_substripe_ranges(model, science_frame=True)["subarray"]
    for i, stripe_id in enumerate(NRCA1_DHS_STRIPE_IDS):
        y0, y1 = sub_ranges[i]
        model.data[:, y0:y1, :] = stripe_id

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nircam.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def nircam_dhs_nrca1_sub164_rate_model(with_wcs=True):
    """
    Create a mock DHS NRCA1 rate CubeModel with subarray SUB164STRIPE4_DHS.

    Parameters
    ----------
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.CubeModel`
        Mock NRCA1 DHS rate model.
    """
    model = datamodels.CubeModel((5, 164, 2048))
    _populate_dhs_shared_metadata(model, subarray="SUB164STRIPE4_DHS")

    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.channel = "SHORT"
    model.meta.instrument.detector = "NRCA1"
    model.meta.instrument.filter = "F150W2"
    model.meta.instrument.pupil = "GDHS0"
    model.meta.instrument.module = "A"

    model.meta.subarray.interleave_reads1 = 1
    model.meta.subarray.multistripe_reads1 = 1
    model.meta.subarray.multistripe_reads2 = 40
    model.meta.subarray.multistripe_skips1 = 1526
    model.meta.subarray.multistripe_skips2 = 85

    model.meta.wcsinfo.siaf_xref_sci = 1024.5
    model.meta.wcsinfo.siaf_yref_sci = 24.5

    # Fill each packed stripe row range with the stripe's ID value so tests
    # can verify that the correct detector region ends up in each output slit.
    sub_ranges = generate_substripe_ranges(model, science_frame=True)["subarray"]
    for i, stripe_id in enumerate(NRCA1_DHS_STRIPE_IDS):
        y0, y1 = sub_ranges[i]
        model.data[:, y0:y1, :] = stripe_id

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nircam.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def nircam_dhs_nrcalong_sub260_rate_model(with_wcs=True):
    """
    Create a mock DHS NRCALONG rate CubeModel with subarray SUB260STRIPE4_DHS.

    Parameters
    ----------
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.CubeModel`
        Mock NRCALONG DHS rate model.
    """
    model = datamodels.CubeModel((5, 260, 2048))
    _populate_dhs_shared_metadata(model)

    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.channel = "LONG"
    model.meta.instrument.detector = "NRCALONG"
    model.meta.instrument.filter = "F322W2"
    model.meta.instrument.pupil = "GRISMR"
    model.meta.instrument.module = "A"

    model.meta.subarray.interleave_reads1 = 0
    model.meta.subarray.multistripe_reads1 = 1
    model.meta.subarray.multistripe_reads2 = 64
    model.meta.subarray.multistripe_skips1 = 959
    model.meta.subarray.multistripe_skips2 = 0

    model.meta.wcsinfo.siaf_xref_sci = 1584.0
    model.meta.wcsinfo.siaf_yref_sci = 32.5

    # Fill each packed stripe row range with the same pattern
    # since for NRCALONG all readouts should be from the same physical detector region
    sub_ranges = generate_substripe_ranges(model, science_frame=True)["subarray"]
    for sub_range in sub_ranges.values():
        y0, y1 = sub_range
        model.data[:, y0 + 5 : y1 - 5, :] = 1.0

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nircam.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model
