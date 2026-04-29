from stdatamodels.jwst.datamodels import RampModel, SuperBiasModel

from jwst.saturation.tests.helpers import add_test_refmodel_metadata

__all__ = ["superbias_model", "nrc_full_ramp", "nrc_subarray_ramp"]


def superbias_model():
    """
    Create a NIRCam superbias model.

    Returns
    -------
    bias_model : `~stdatamodels.jwst.datamodels.SuperbiasModel`
        The superbias datamodel.
    """
    bias_model = SuperBiasModel((2048, 2048))
    bias_model.meta.subarray.xstart = 1
    bias_model.meta.subarray.ystart = 1
    bias_model.meta.subarray.xsize = 2048
    bias_model.meta.subarray.ysize = 2048
    bias_model.meta.instrument.name = "NIRCAM"
    add_test_refmodel_metadata(bias_model)
    return bias_model


def nrc_full_ramp(ngroups, nrows, ncols):
    """
    Set up mock NIRCam FULL data to test.

    Returns
    -------
    data_model : `~stdatamodels.jwst.datamodels.RampModel`
        The NIRCam ramp model.
    """
    nints = 1

    # create a JWST datamodel for NIRCam FULL data
    data_model = RampModel((nints, ngroups, nrows, ncols))
    data_model.meta.subarray.xstart = 1
    data_model.meta.subarray.ystart = 1
    data_model.meta.subarray.xsize = ncols
    data_model.meta.subarray.ysize = nrows
    data_model.meta.exposure.ngroups = ngroups
    data_model.meta.instrument.name = "NIRCAM"
    data_model.meta.instrument.detector = "NRCA1"
    data_model.meta.observation.date = "2017-10-01"
    data_model.meta.observation.time = "00:00:00"

    return data_model


def nrc_subarray_ramp(xstart, ystart, ngroups, nrows, ncols):
    """
    Set up mock NIRCam subarray data to test.

    Returns
    -------
    data_model : `~stdatamodels.jwst.datamodels.RampModel`
        The NIRCam ramp model.
    """
    nints = 1

    # create a JWST datamodel for NIRCam SUB320A335R data
    data_model = RampModel((nints, ngroups, nrows, ncols))
    data_model.meta.subarray.name = "SUB320A335R"
    data_model.meta.subarray.xstart = xstart
    data_model.meta.subarray.ystart = ystart
    data_model.meta.subarray.xsize = ncols
    data_model.meta.subarray.ysize = nrows
    data_model.meta.exposure.ngroups = ngroups
    data_model.meta.instrument.name = "NIRCAM"
    data_model.meta.instrument.detector = "NRCALONG"
    data_model.meta.observation.date = "2019-10-14"
    data_model.meta.observation.time = "16:44:12.000"

    return data_model
