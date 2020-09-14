import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def get_dispersion_direction(exposure_type, grating="ANY", filter_wh="ANY",
                             pupil="ANY"):
    """Get the dispersion direction.

    Parameters
    ----------
    exposure_type : str
        The exposure type.

    grating : str
        The name of the optic in the grating wheel.

    filter_wh : str
        The name of the optic in the filter wheel.

    pupil : str
        The name of the optic in the pupil wheel.

    Returns
    -------
    int or None : The dispersion direction

-    None  the dispersion direction is not meaningful or not defined
-       1  the dispersion direction is horizontal ("sky" coordinates)
-       2  the dispersion direction is vertical ("sky" coordinates)
    """

    exposure_type = exposure_type.upper()
    if grating is not None:
        grating = grating.upper()
    if filter_wh is not None:
        filter_wh = filter_wh.upper()
    if pupil is not None:
        pupil = pupil.upper()

    # The keys of the `by_exp_type` dictionary are values of exposure type.
    # If the dispersion direction is uniquely determined by the exposure
    # type, the associated value of by_exp_type[exposure_type] will be an
    # integer, the dispersion direction.  If one or more parameters other
    # than the exposure type are needed, the value of
    # by_exp_type[exposure_type] will be a tuple of several strings; use
    # that tuple as the key for dictionary `second_pass`, and the value of
    # that will be the dispersion direction.
    by_exp_type = {
        # FGS science
        "FGS_DARK": None,
        "FGS_FOCUS": None,
        "FGS_IMAGE": None,
        "FGS_INTFLAT": None,
        "FGS_SKYFLAT": None,
        # FGS guide star
        "FGS_ACQ1": None,
        "FGS_ACQ2": None,
        "FGS_FINEGUIDE": None,
        "FGS_ID-IMAGE": None,
        "FGS_ID-STACK": None,
        "FGS_TRACK": None,
        # MIRI
        "MIR_4QPM": None,
        "MIR_CORONCAL": None,
        "MIR_DARKALL": None,
        "MIR_DARKIMG": None,
        "MIR_DARKMRS": None,
        "MIR_FLATALL": None,
        "MIR_FLATIMAGE": None,
        "MIR_FLATIMAGE-EXT": None,
        "MIR_FLATMRS": 2,
        "MIR_FLATMRS-EXT": 2,
        "MIR_IMAGE": None,
        "MIR_LRS-FIXEDSLIT": 2,
        "MIR_LRS-SLITLESS": 2,
        "MIR_LYOT": None,
        "MIR_MRS": 2,
        "MIR_TACQ": None,
        # NIRISS
        "NIS_AMI": None,
        "NIS_DARK": None,
        "NIS_EXTCAL": None,
        "NIS_FOCUS": None,
        "NIS_IMAGE": None,
        "NIS_LAMP": None,
        "NIS_SOSS": 1,
        "NIS_TACQ": None,
        "NIS_TACONFIRM": None,
        "NIS_WFSS": (exposure_type, "ANY", filter_wh, "ANY"),
        # NIRCam
        "NRC_CORON": None,
        "NRC_DARK": None,
        "NRC_FLAT": None,
        "NRC_FOCUS": None,
        "NRC_GRISM": (exposure_type, "ANY", "ANY", pupil),
        "NRC_IMAGE": None,
        "NRC_WFSS": (exposure_type, "ANY", "ANY", pupil),
        "NRC_LED": None,
        "NRC_WFSC": None,
        "NRC_TACONFIRM": None,
        "NRC_TACQ": None,
        "NRC_TSGRISM": (exposure_type, "ANY", "ANY", pupil),
        "NRC_TSIMAGE": None,
        # NIRSpec
        "NRS_AUTOFLAT": 1,
        "NRS_AUTOWAVE": 1,
        "NRS_BRIGHTOBJ": 1,
        "NRS_CONFIRM": None,
        "NRS_DARK": None,
        "NRS_FIXEDSLIT": 1,
        "NRS_FOCUS": None,
        "NRS_IFU": 1,
        "NRS_IMAGE": None,
        "NRS_LAMP":  1,
        "NRS_MIMF": None,
        "NRS_MSASPEC": 1,
        "NRS_MSATA": None,
        "NRS_TACONFIRM": None,
        "NRS_TACQ": None,
        "NRS_TASLIT": None,
        "NRS_WATA": None,
        # Misc
        "N/A": None,
        "ANY": None
    }

    if exposure_type not in by_exp_type.keys():
        return None

    # The strings in each tuple are exposure_type, grating, filter_wh, pupil.
    second_pass = {
        ("NIS_WFSS", "ANY", "GR150R", "ANY"): 2,
        ("NIS_WFSS", "ANY", "GR150C", "ANY"): 1,

        ("NRC_GRISM", "ANY", "ANY", "GRISMR"): 1,
        ("NRC_GRISM", "ANY", "ANY", "GRISMC"): 2,

        ("NRC_TSGRISM", "ANY", "ANY", "GRISMR"): 1,

        ("NRC_WFSS", "ANY", "ANY", "GRISMR"): 1,
        ("NRC_WFSS", "ANY", "ANY", "GRISMC"): 2
    }

    select = by_exp_type[exposure_type]
    if isinstance(select, int):
        return select
    else:
        if select in second_pass.keys():
            return second_pass[select]
        else:
            log.warning("Error in get_dispersion_direction:  {} not in "
                        "`second_pass`".format(select))
            log.warning("Dispersion direction could not be determined.")
            return None
