"""Module for lists of modes grouped in different ways."""

from jwst.associations.lib.dms_base import (
    ACQ_EXP_TYPES,
    IMAGE2_SCIENCE_EXP_TYPES,
    IMAGE2_NONSCIENCE_EXP_TYPES,
    SPEC2_SCIENCE_EXP_TYPES,
)

IMAGING_TYPES = set(
    tuple(ACQ_EXP_TYPES)
    + tuple(IMAGE2_SCIENCE_EXP_TYPES)
    + tuple(IMAGE2_NONSCIENCE_EXP_TYPES)
    + ("fgs_image", "fgs_focus")
)

SPEC_TYPES = SPEC2_SCIENCE_EXP_TYPES

# FGS guide star exposures
FGS_GUIDE_EXP_TYPES = [
    "fgs_acq1",
    "fgs_acq2",
    "fgs_fineguide",
    "fgs_id-image",
    "fgs_id-stack",
    "fgs_track",
]

# NIRSPEC lamp mode spec types
NRS_LAMP_MODE_SPEC_TYPES = [
    "brightobj",
    "fixedslit",
    "ifu",
    "msaspec",
]


def is_nrs_lamp(datamodel):
    """
    Check if exposure type is nrs_lamp or nrs_autowave.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    exp_type = datamodel.meta.exposure.type.lower()
    return exp_type in ["nrs_lamp", "nrs_autowave"]


def is_nrs_linelamp(datamodel):
    """
    Check if lamp state is lin or ref.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    lamp_state = datamodel.meta.instrument.lamp_state.lower()
    return lamp_state[0:3] in ["lin", "ref"]


def is_nrs_flatlamp(datamodel):
    """
    Check if lamp state is flat.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    lamp_state = datamodel.meta.instrument.lamp_state.lower()
    return lamp_state[0:4] == "flat"


def is_nrs_slit_linelamp(datamodel):
    """
    Check NRS slit line lamp.

    Specifically, check if lamp mode is msaspec, fixedslit, or brightobj
    when lamp state is lin or ref and exposure type is
    nrs_autowave or nrs_lamp.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    exp_type = datamodel.meta.exposure.type.lower()
    return (
        lamp_mode in ["msaspec", "fixedslit", "brightobj"]
        and is_nrs_linelamp(datamodel)
        and exp_type in ["nrs_autowave", "nrs_lamp"]
    )


def is_nrs_ifu_linelamp(datamodel):
    """
    Check if lamp mode is ifu and lamp state is lin or ref.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == "ifu" and is_nrs_linelamp(datamodel)


def is_nrs_ifu_flatlamp(datamodel):
    """
    Check if lamp mode is ifu and lamp state is flat.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == "ifu" and is_nrs_flatlamp(datamodel)


def is_nrs_ifu_lamp(datamodel):
    """
    Check if lamp mode is ifu and exposure type is nrs_lamp or nrs_autowave.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == "ifu" and is_nrs_lamp(datamodel)


def is_nrs_msaspec_lamp(datamodel):
    """
    Check if lamp mode is msaspec and exposure type is nrs_lamp or nrs_autowave.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == "msaspec" and is_nrs_lamp(datamodel)


def is_nrs_msaspec_linelamp(datamodel):
    """
    Check if lamp mode is msaspec and lamp state is lin or ref.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == "msaspec" and is_nrs_linelamp(datamodel)


def is_nrs_msaspec_flatlamp(datamodel):
    """
    Check if lamp mode is msaspec and lamp state is flat.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == "msaspec" and is_nrs_flatlamp(datamodel)


def is_nrs_autoflat(datamodel):
    """
    Check of exposure type is nrs_autoflat.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    exp_type = datamodel.meta.exposure.type
    return exp_type.lower() == "nrs_autoflat"


def is_moving_target(datamodel):
    """
    Determine if a moving target exposure.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`
        JWST data model.

    Returns
    -------
    status : bool
        `True` if it is.
    """
    if (
        hasattr(datamodel.meta.target, "type")
        and datamodel.meta.target.type is not None
        and datamodel.meta.target.type.lower() == "moving"
    ):
        return True
    return False
