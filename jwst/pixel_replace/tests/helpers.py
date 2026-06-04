import astropy.units as u
import gwcs
import numpy as np
from astropy.modeling.models import Const1D, Mapping
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels.dqflags import pixel as flags

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_ifu_file

__all__ = [
    "cal_data",
    "nirspec_tso",
    "nirspec_fs_slitmodel",
    "nirspec_msa_multislit",
    "nirspec_ifu",
]


def cal_data(shape, bad_idx, dispaxis=1, model="slit"):
    """
    Make a basic cal model with a bad pixel.

    Parameters
    ----------
    shape : tuple
        Data shape.
    bad_idx : tuple
        Bad pixel location.
    dispaxis : int
        Dispersion axis.
    model : str, optional
        Model type.  May be "slit", "image", or "ifu".

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`, \
            `~stdatamodels.jwst.datamodels.IFUImageModel`, or \
            `~stdatamodels.jwst.datamodels.SlitModel`
        The data model.
    """
    if model == "image":
        model = datamodels.ImageModel(shape)
    elif model == "ifu":
        model = datamodels.IFUImageModel(shape)
    else:
        model = datamodels.SlitModel(shape)
    model.meta.wcsinfo.dispersion_direction = dispaxis

    # Set the data and error arrays to all 1s except one bad pixel
    # to correct at the middle of the array
    ones = np.ones(shape, dtype=float)
    model.data = ones.copy()
    model.dq = model.get_default("dq")
    model.err = ones.copy()
    model.var_poisson = ones.copy()
    model.var_rnoise = ones.copy()
    model.var_flat = ones.copy()

    bad_flag = flags["DO_NOT_USE"] + flags["OTHER_BAD_PIXEL"]
    model.data[bad_idx] = np.nan
    model.err[bad_idx] = np.nan
    model.var_poisson[bad_idx] = np.nan
    model.var_rnoise[bad_idx] = np.nan
    model.var_flat[bad_idx] = np.nan
    model.dq[bad_idx] = bad_flag

    # Also add a non-science region in one row and one column
    non_science = flags["DO_NOT_USE"] + flags["NON_SCIENCE"]
    model.data[..., 1] = np.nan
    model.err[..., 1] = np.nan
    model.var_poisson[..., 1] = np.nan
    model.var_rnoise[..., 1] = np.nan
    model.var_flat[..., 1] = np.nan
    model.dq[..., 1] = non_science

    model.data[..., 1, :] = np.nan
    model.err[..., 1, :] = np.nan
    model.var_poisson[..., 1, :] = np.nan
    model.var_rnoise[..., 1, :] = np.nan
    model.var_flat[..., 1, :] = np.nan
    model.dq[..., 1, :] = non_science

    return model


def nirspec_tso():
    """
    Make a NIRSpec TSO model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.SlitModel`
        The data model.
    bad_idx : tuple
        The bad pixel location.
    """
    bad_idx = (1, 10, 10)
    model = cal_data(shape=(3, 20, 20), bad_idx=bad_idx, dispaxis=1)
    model.meta.instrument.name = "NIRSPEC"
    model.meta.exposure.type = "NRS_BRIGHTOBJ"
    return model, bad_idx


def nirspec_fs_slitmodel():
    """
    Make a NIRSpec FS model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.SlitModel`
        The data model.
    bad_idx : tuple
        The bad pixel location.
    """
    bad_idx = (10, 10)
    model = cal_data(shape=(20, 20), bad_idx=bad_idx, dispaxis=1)
    model.meta.instrument.name = "NIRSPEC"
    model.meta.exposure.type = "NRS_FIXEDSLIT"
    return model, bad_idx


def nirspec_msa_multislit():
    """
    Make a NIRSpec MSA model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        The data model.
    bad_idx : tuple
        The bad pixel location.
    """
    bad_idx = (10, 10)
    slit_model = cal_data(shape=(20, 20), bad_idx=bad_idx, dispaxis=1)
    model = datamodels.MultiSlitModel()
    model.slits.append(slit_model)
    model.meta.instrument.name = "NIRSPEC"
    model.meta.exposure.type = "NRS_MSASPEC"
    return model, bad_idx


def nirspec_ifu():
    """
    Make a NIRSpec IFU model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The data model.
    bad_idx : tuple
        The bad pixel location.
    """
    shape = (2048, 2048)
    bad_idx = (1424, 690)

    # IFU mode requires WCS information, so make a more realistic model
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    hdul["SCI"].data = np.ones(shape, dtype=float)

    model = datamodels.IFUImageModel(hdul)
    model = AssignWcsStep.call(model)

    test_data = cal_data(shape=shape, bad_idx=bad_idx, dispaxis=1, model="ifu")
    model.data = test_data.data
    model.dq = test_data.dq
    model.err = test_data.err
    model.var_poisson = test_data.var_poisson
    model.var_rnoise = test_data.var_rnoise
    model.var_flat = test_data.var_flat

    test_data.close()

    return model, bad_idx


def miri_lrs():
    """
    Make a MIRI LRS fixed slit model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The data model.
    bad_idx : tuple
        The bad pixel location.
    """
    bad_idx = (10, 10)
    model = cal_data(shape=(20, 20), bad_idx=bad_idx, dispaxis=2, model="image")
    model.meta.instrument.name = "MIRI"
    model.meta.exposure.type = "MIR_LRS-FIXEDSLIT"
    return model, bad_idx


def miri_mrs():
    """
    Make a MIRI MRS model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The data model.
    bad_idx : tuple
        The bad pixel location.
    """
    shape = (20, 20)
    bad_idx = (10, 10)
    model = cal_data(shape=shape, bad_idx=bad_idx, dispaxis=2, model="ifu")
    model.meta.instrument.name = "MIRI"
    model.meta.exposure.type = "MIR_MRS"

    # Mock a wcs that just returns 1 for alpha, beta, lam, region
    transform = Mapping((0, 1, 1), n_inputs=2) | Const1D(1) & Const1D(1) & Const1D(1)
    label_mapper = gwcs.selector.LabelMapperArray(np.full(shape, 1))
    det2alpha_beta = gwcs.selector.RegionsSelector(
        ("x", "y"), ("alpha", "beta", "lam"), label_mapper=label_mapper, selector={1: transform}
    )
    output_frame = gwcs.CompositeFrame(
        [
            gwcs.Frame2D(name="alpha_beta_spatial", axes_order=(0, 1), unit=(u.arcsec, u.arcsec)),
            gwcs.SpectralFrame(name="lam", axes_order=(2,), unit=(u.um,)),
        ],
        name="alpha_beta",
    )
    model.meta.wcs = gwcs.WCS(
        [(gwcs.Frame2D(name="detector"), det2alpha_beta), (output_frame, None)]
    )
    return model, bad_idx
