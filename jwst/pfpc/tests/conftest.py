import pytest

from jwst.pfpc.tests import helpers


@pytest.fixture(scope="module")
def mrs_dith1_ch12_medium():
    """
    Make a MIRI MRS cal model with channel 1,2 and band MEDIUM, dither index 1.

    Yields
    ------
    `~stdatamodels.jwst.datamodels.IFUImageModel`
        The MRS model.
    """
    model = helpers.dithered_miri_mrs_model(detector="MIRIFUSHORT", channel="12", band="MEDIUM")
    yield model
    model.close()


@pytest.fixture(scope="module")
def mrs_dith2_ch12_medium():
    """
    Make a MIRI MRS cal model with channel 1,2 and band MEDIUM, dither index 2.

    Yields
    ------
    `~stdatamodels.jwst.datamodels.IFUImageModel`
        The MRS model.
    """
    model = helpers.dithered_miri_mrs_model(
        detector="MIRIFUSHORT", channel="12", band="MEDIUM", patt_num=2
    )
    yield model
    model.close()


@pytest.fixture(scope="module")
def mrs_dith1_ch34_short():
    """
    Make a MIRI MRS cal model with channel 3,4 and band SHORT.

    Yields
    ------
    `~stdatamodels.jwst.datamodels.IFUImageModel`
        The MRS model.
    """
    model = helpers.dithered_miri_mrs_model(detector="MIRIFULONG", channel="34", band="SHORT")
    yield model
    model.close()


@pytest.fixture(scope="module")
def nrs_ifu_dith1():
    """
    Make a NIRSpec IFU cal model with dither index 1.

    Yields
    ------
    `~stdatamodels.jwst.datamodels.IFUImageModel`
        The IFU model.
    """
    model = helpers.dithered_nirspec_ifu_model()
    yield model
    model.close()


@pytest.fixture(scope="function")
def mrs_pfpc_model():
    """
    Make a PFPC reference model for MIRI MRS.

    Note: this fixture must be function-scoped because ``model.copy()``
    does not perform a deepcopy on the FITS table.

    Yields
    ------
    `~stdatamodels.jwst.datamodels.MirMrsPFPCModel`
        The PFPC datamodel.
    """
    pfpc = helpers.miri_mrs_pfpc_model()
    yield pfpc
    pfpc.close()


@pytest.fixture(scope="function")
def nrs_ifu_pfpc_model():
    """
    Make a PFPC reference model for NIRSpec IFU.

    Note: this fixture must be function-scoped because ``model.copy()``
    does not perform a deepcopy on the FITS table.

    Yields
    ------
    `~stdatamodels.jwst.datamodels.MirMrsPFPCModel`
        The PFPC datamodel. Borrows the MIRI MRS model type for now.
    """
    pfpc = helpers.nirspec_ifu_pfpc_model()
    yield pfpc
    pfpc.close()


@pytest.fixture(scope="module")
def mrs_dith1_ch1_medium_x1d():
    """
    Make a MIRI MRS x1d model with channel 1,2 and band MEDIUM, dither index 1.

    Yields
    ------
    `~stdatamodels.jwst.datamodels.IFUImageModel`
        The MRS model.
    """
    model = helpers.miri_mrs_x1d(detector="MIRIFUSHORT", channel="1", band="MEDIUM")
    yield model
    model.close()
