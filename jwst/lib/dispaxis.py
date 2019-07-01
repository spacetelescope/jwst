import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def get_dispersion_direction(exposure_type, grating="ANY", filter="ANY",
                             pupil="ANY"):
    """Get the dispersion direction.

    Parameters
    ----------
    exposure_type : str
        The exposure type.

    grating : str
        The name of the optic in the grating wheel.

    filter : str
        The name of the optic in the filter wheel.

    pupil : str
        The name of the optic in the pupil wheel.

    Returns
    -------
    int : The dispersion direction

-       0  the dispersion direction is not meaningful or not defined
-       1  the dispersion direction is horizontal ("sky" coordinates)
-       2  the dispersion direction is vertical ("sky" coordinates)
    """

    exposure_type = exposure_type.upper()
    grating = grating.upper()
    filter = filter.upper()
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
        "FGS_DARK": 0,
        "FGS_FOCUS": 0,
        "FGS_IMAGE": 0,
        "FGS_INTFLAT": 0,
        "FGS_SKYFLAT": 0,
        # FGS guide star
        "FGS_ACQ1": 0,
        "FGS_ACQ2": 0,
        "FGS_FINEGUIDE": 0,
        "FGS_ID-IMAGE": 0,
        "FGS_ID-STACK": 0,
        "FGS_TRACK": 0,
        # MIRI
        "MIR_4QPM": 0,
        "MIR_CORONCAL": 0,
        "MIR_DARKALL": 0,
        "MIR_DARKIMG": 0,
        "MIR_DARKMRS": 0,
        "MIR_FLATALL": 0,
        "MIR_FLATIMAGE": 0,
        "MIR_FLATIMAGE-EXT": 0,
        "MIR_FLATMRS": 2,
        "MIR_FLATMRS-EXT": 2,
        "MIR_IMAGE": 0,
        "MIR_LRS-FIXEDSLIT": 2,
        "MIR_LRS-SLITLESS": 2,
        "MIR_LYOT": 0,
        "MIR_MRS": 2,
        "MIR_TACQ": 0,
        # NIRISS
        "NIS_AMI": 0,
        "NIS_DARK": 0,
        "NIS_EXTCAL": 0,
        "NIS_FOCUS": 0,
        "NIS_IMAGE": 0,
        "NIS_LAMP": 0,
        "NIS_SOSS": 1,
        "NIS_TACQ": 0,
        "NIS_TACONFIRM": 0,
        "NIS_WFSS": (exposure_type, "ANY", filter, "ANY"),
        # NIRCam
        "NRC_CORON": 0,
        "NRC_DARK": 0,
        "NRC_FLAT": 0,
        "NRC_FOCUS": 0,
        "NRC_GRISM": (exposure_type, "ANY", "ANY", pupil),
        "NRC_IMAGE": 0,
        "NRC_WFSS": (exposure_type, "ANY", "ANY", pupil),
        "NRC_LED": 0,
        "NRC_WFSC": 0,
        "NRC_TACONFIRM": 0,
        "NRC_TACQ": 0,
        "NRC_TSGRISM": (exposure_type, "ANY", "ANY", pupil),
        "NRC_TSIMAGE": 0,
        # NIRSpec
        "NRS_AUTOFLAT": 0,
        "NRS_AUTOWAVE": 1,
        "NRS_BRIGHTOBJ": 1,
        "NRS_CONFIRM": 0,
        "NRS_DARK": 0,
        "NRS_FIXEDSLIT": 1,
        "NRS_FOCUS": 0,
        "NRS_IFU": 1,
        "NRS_IMAGE": 0,
        "NRS_LAMP":  1,
        "NRS_MIMF": 0,
        "NRS_MSASPEC": 1,
        "NRS_MSATA": 0,
        "NRS_TACONFIRM": 0,
        "NRS_TACQ": 0,
        "NRS_TASLIT": 0,
        "NRS_WATA": 0,
        # Misc
        "N/A": 0,
        "ANY": 0
    }

    if exposure_type not in by_exp_type.keys():
        return 0

    # The strings in each tuple are exposure_type, grating, filter, pupil.
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
            return 0
