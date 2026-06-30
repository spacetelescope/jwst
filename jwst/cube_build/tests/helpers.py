import numpy as np
from stdatamodels.jwst.datamodels import MiriIFUCubeParsModel


def make_miri_cube_pars():
    """
    Set up the miri cube pars reference file.

    Returns
    -------
    model : MiriIFUCubeParsModel
        The MIRI IFU cube pars model.
    """
    model = MiriIFUCubeParsModel()
    model.meta.reftype = "CUBEPAR"
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "N/A"
    model.meta.exposure.type = "MIR_MRS"

    # make the first extension
    channel = np.array(["1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4"])
    subchannel = np.array(
        [
            "SHORT",
            "MEDIUM",
            "LONG",
            "SHORT",
            "MEDIUM",
            "LONG",
            "SHORT",
            "MEDIUM",
            "LONG",
            "SHORT",
            "MEDIUM",
            "LONG",
        ]
    )

    spsize = np.array([0.13, 0.13, 0.13, 0.17, 0.17, 0.17, 0.2, 0.2, 0.2, 0.35, 0.35, 0.35])
    wsamp = np.array(
        [0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.003, 0.003, 0.003, 0.006, 0.006, 0.006]
    )

    wmin = np.array([4.89, 5.65, 6.52, 7.49, 8.65, 10.00, 11.53, 13.37, 15.44, 17.66, 20.54, 23.95])
    wmax = np.array([5.75, 6.64, 7.66, 8.78, 10.14, 11.7, 13.48, 15.63, 18.05, 20.92, 24.40, 28.45])

    dtype1 = model.get_dtype("ifucubepars_table")
    model.ifucubepars_table = np.array(
        list(zip(channel, subchannel, spsize, wsamp, wmin, wmax, strict=True)), dtype=dtype1
    )

    # make the second extension
    roispat = np.array([0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.20, 0.20, 0.20, 0.40, 0.40, 0.40])
    roispec = np.array(
        [0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.003, 0.003, 0.003, 0.006, 0.006, 0.006]
    )

    power = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
    softrad = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])

    dtype2 = model.get_dtype("ifucubepars_msm_table")
    model.ifucubepars_msm_table = np.array(
        list(zip(channel, subchannel, roispat, roispec, power, softrad, strict=True)), dtype=dtype2
    )

    # need EMSM to be defined for test_invalid_coord_sys because setting to internal_cal
    # will hard-code weighting to emsm
    dtype4 = model.get_dtype("ifucubepars_emsm_table")
    model.ifucubepars_emsm_table = np.array(
        list(zip(channel, subchannel, roispat, roispec, softrad, strict=True)), dtype=dtype4
    )

    # make the third extension
    # Define the multiextension wavelength solution - only use a few number for testing
    finalwave = np.array([5, 10, 15, 20, 25])
    roispat = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    roispec = np.array([0.001, 0.002, 0.003, 0.004, 0.005])
    power = np.array([1, 2, 3, 4, 5])
    softrad = np.array([0.01, 0.02, 0.03, 0.04, 0.05])

    dtype3 = model.get_dtype("ifucubepars_multichannel_msm_wavetable")
    model.ifucubepars_multichannel_msm_wavetable = np.array(
        list(zip(finalwave, roispat, roispec, power, softrad, strict=True)), dtype=dtype3
    )

    # need EMSM to be defined for test_invalid_coord_sys because setting to internal_cal
    # will hard-code weighting to emsm
    dtype5 = model.get_dtype("ifucubepars_multichannel_emsm_wavetable")
    model.ifucubepars_multichannel_emsm_wavetable = np.array(
        list(zip(finalwave, roispat, roispec, softrad, strict=True)), dtype=dtype5
    )

    # set units of all tables to default

    for table_name in [
        "ifucubepars_table",
        "ifucubepars_msm_table",
        "ifucubepars_multichannel_msm_wavetable",
        "ifucubepars_emsm_table",
        "ifucubepars_multichannel_emsm_wavetable",
    ]:
        table = getattr(model, table_name)
        unit_table_name = table_name + "_units"
        for name in table.dtype.names:
            default_unit = getattr(model, unit_table_name).get_default(name)
            setattr(getattr(model, unit_table_name), name, default_unit)

    return model
