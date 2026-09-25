import numpy as np
from astropy.table import Table, vstack
from stdatamodels.jwst import datamodels

from jwst.adaptive_trace_model.tests.helpers import miri_mrs_model


def dithered_miri_mrs_model(detector="MIRIFUSHORT", channel="12", band="LONG", patt_num=1):
    """
    Make a dithered MIRI MRS model suitable for PFPC correction.

    Parameters
    ----------
    detector : str, optional
        Detector name.
    channel : str, optional
        Channel name.
    band : str, optional
        Band name.
    patt_num : int
        Dither index.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The MRS datamodel.
    """
    model = miri_mrs_model(detector=detector, channel=channel, band=band)

    # May vary for input
    model.meta.dither.position_number = patt_num

    # Set to match PFPC file
    model.meta.dither.primary_type = "4-POINT"
    model.meta.dither.direction = "NEGATIVE"
    model.meta.dither.optimized_for = "POINT-SOURCE"
    model.meta.dither.primary_channel = "ALL_MRS"

    # Must be point source to be processed
    model.meta.target.source_type = "POINT"

    # combine_1d needs exposure time
    model.meta.exposure.exposure_time = 1.0

    # Set a test filename
    model.meta.filename = "test_mrs_cal.fits"

    return model


def miri_mrs_pfpc_model():
    """
    Make a PFPC reference model for MIRI MRS.

    Includes placeholder corrections for all channels and bands in a
    4 point dither pattern. The value for the correction is set to
    channel + dither index.

    Returns
    -------
    pfpc_model : `~stdatamodels.jwst.datamodels.MirMrsPFPCModel`
        The PFPC datamodel.
    """
    pfpc_tab = None
    ndither = 4
    for channel in range(1, 5):
        for band in ["SHORT", "MEDIUM", "LONG"]:
            new_tab = Table()
            new_tab["channel"] = [str(channel)] * ndither
            new_tab["band"] = [band] * ndither
            new_tab["patttype"] = ["4-POINT"] * ndither
            new_tab["dithdirc"] = ["NEGATIVE"] * ndither
            new_tab["dithopfr"] = ["POINT-SOURCE"] * ndither
            new_tab["mrsprchn"] = ["ALL_MRS"] * ndither
            new_tab["patt_num"] = np.arange(1, ndither + 1)
            new_tab["r_flat"] = ["N/A"] * ndither
            new_tab["r_photom"] = ["N/A"] * ndither

            # Linear wavelength sampling over a large range
            new_tab["wavelength"] = [np.linspace(4, 30, 1000)] * ndither

            # Set the correction to the value of the channel + dither number, for testing
            new_tab["correction"] = [np.full(1000, channel) + i for i in range(ndither)]

            if pfpc_tab is None:
                pfpc_tab = new_tab
            else:
                pfpc_tab = vstack([pfpc_tab, new_tab], join_type="exact")

    pfpc_model = datamodels.MirMrsPFPCModel()
    pfpc_model.pfpc_table = pfpc_tab.as_array()
    return pfpc_model
