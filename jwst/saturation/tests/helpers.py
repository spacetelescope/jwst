import numpy as np
from stdatamodels.jwst.datamodels import RampModel, SaturationModel, SuperBiasModel

from jwst.dq_init.tests.helpers import make_superstripe_model

__all__ = [
    "add_test_refmodel_metadata",
    "setup_nrc_cube",
    "setup_miri_cube",
    "setup_nrs_irs2_cube",
    "setup_nrs_nrs_cube",
    "setup_nis_superstripe_cube",
]


def add_test_refmodel_metadata(refmodel):
    """
    Add minimal mock metadata for reference files.

    Parameters
    ----------
    refmodel : `~stdatamodels.jwst.datamodels.DataModel`
        The reference datamodel, updated in place.
    """
    refmodel.meta.telescope = "JWST"
    refmodel.meta.description = "filler"
    refmodel.meta.reftype = "filler"
    refmodel.meta.author = "Py Test"
    refmodel.meta.pedigree = "Pytest"
    refmodel.meta.useafter = "2015-01-01T01:00:00"


def setup_nrc_cube(ngroups, nrows, ncols):
    """
    Set up mock NIRCam data to test.

    The number of integrations is always 1.

    Parameters
    ----------
    ngroups : int
        Number of groups.
    nrows : int
        Number of rows.
    ncols : int
        Number of columns.

    Returns
    -------
    data_model : `~stdatamodels.jwst.datamodels.RampModel`
        The NIRCam ramp model.
    saturation_model : `~stdatamodels.jwst.datamodels.SaturationModel`
        A corresponding saturation model.
    """
    nints = 1

    data_model = RampModel((nints, ngroups, nrows, ncols))
    data_model.meta.subarray.xstart = 1
    data_model.meta.subarray.ystart = 1
    data_model.meta.subarray.xsize = ncols
    data_model.meta.subarray.ysize = nrows
    data_model.meta.exposure.ngroups = ngroups
    data_model.meta.exposure.nframes = 3
    data_model.meta.instrument.name = "NIRCAM"
    data_model.meta.instrument.detector = "NRCA1"
    data_model.meta.observation.date = "2017-10-01"
    data_model.meta.observation.time = "00:00:00"

    saturation_model = SaturationModel((2048, 2048))
    saturation_model.meta.subarray.xstart = 1
    saturation_model.meta.subarray.ystart = 1
    saturation_model.meta.subarray.xsize = 2048
    saturation_model.meta.subarray.ysize = 2048
    saturation_model.meta.instrument.name = "NIRCAM"
    add_test_refmodel_metadata(saturation_model)

    return data_model, saturation_model


def setup_miri_cube(xstart, ystart, ngroups, nrows, ncols):
    """
    Set up mock MIRI data to test.

    The number of integrations is always 1.

    Parameters
    ----------
    xstart : int
        Starting x-value for the subarray.
    ystart : int
        Starting y-value for the subarray.
    ngroups : int
        Number of groups.
    nrows : int
        Number of rows.
    ncols : int
        Number of columns.

    Returns
    -------
    data_model : `~stdatamodels.jwst.datamodels.RampModel`
        The MIRI ramp model.
    saturation_model : `~stdatamodels.jwst.datamodels.SaturationModel`
        A corresponding saturation model.
    """
    nints = 1

    # create a JWST datamodel for MIRI data
    data_model = RampModel((nints, ngroups, nrows, ncols))
    data_model.data += 1
    data_model.meta.instrument.name = "MIRI"
    data_model.meta.instrument.detector = "MIRIMAGE"
    data_model.meta.instrument.filter = "F1500W"
    data_model.meta.instrument.band = "N/A"
    data_model.meta.observation.date = "2016-06-01"
    data_model.meta.observation.time = "00:00:00"
    data_model.meta.exposure.type = "MIR_IMAGE"
    data_model.meta.subarray.name = "MASK1550"
    data_model.meta.subarray.xstart = xstart
    data_model.meta.subarray.xsize = ncols
    data_model.meta.subarray.ystart = ystart
    data_model.meta.subarray.ysize = nrows
    data_model.meta.exposure.ngroups = ngroups
    data_model.meta.exposure.nframes = 3

    # create a saturation model for the saturation step
    saturation_model = SaturationModel((1032, 1024))
    saturation_model.meta.instrument.name = "MIRI"
    saturation_model.meta.instrument.detector = "MIRIMAGE"
    saturation_model.meta.subarray.xstart = 1
    saturation_model.meta.subarray.xsize = 1024
    saturation_model.meta.subarray.ystart = 1
    saturation_model.meta.subarray.ysize = 1032
    add_test_refmodel_metadata(saturation_model)

    return data_model, saturation_model


def setup_nrs_irs2_cube():
    """
    Set up mock NIRSpec IRS2 data to test.

    The output data shape is (1, 5, 3200, 2048).

    Returns
    -------
    data_model : `~stdatamodels.jwst.datamodels.RampModel`
        The NIRSpec IRS2 ramp model.
    saturation_model : `~stdatamodels.jwst.datamodels.SaturationModel`
        A corresponding saturation model.
    bias_model : `~stdatamodels.jwst.datamodels.SuperbiasModel`
        A corresponding superbias model.
    """
    # create a JWST datamodel for NIRSPEC IRS2 data
    data_model = RampModel((1, 5, 3200, 2048))
    data_model.data = np.ones((1, 5, 3200, 2048))
    data_model.groupdq = np.zeros((1, 5, 3200, 2048))
    data_model.pixeldq = np.zeros((3200, 2048))
    data_model.meta.instrument.name = "NIRSPEC"
    data_model.meta.instrument.detector = "NRS1"
    data_model.meta.instrument.filter = "F100LP"
    data_model.meta.observation.date = "2019-07-19"
    data_model.meta.observation.time = "23:23:30.912"
    data_model.meta.exposure.type = "NRS_LAMP"
    data_model.meta.subarray.name = "FULL"
    data_model.meta.subarray.xstart = 1
    data_model.meta.subarray.xsize = 2048
    data_model.meta.subarray.ystart = 1
    data_model.meta.subarray.ysize = 2048
    data_model.meta.exposure.nrs_normal = 16
    data_model.meta.exposure.nrs_reference = 4
    data_model.meta.exposure.readpatt = "NRSIRS2"
    data_model.meta.exposure.nframes = 5
    data_model.meta.exposure.ngroups = 5

    # create a saturation model for the saturation step
    saturation_model = SaturationModel((2048, 2048))
    saturation_model.data = (
        np.ones((2048, 2048)) * 60000
    )  # saturation limit for every pixel is 60000
    saturation_model.meta.subarray.xstart = 1
    saturation_model.meta.subarray.xsize = 2048
    saturation_model.meta.subarray.ystart = 1
    saturation_model.meta.subarray.ysize = 2048
    add_test_refmodel_metadata(saturation_model)

    # create a bias model for group 2 saturation checking
    bias_model = SuperBiasModel((3200, 2048))
    bias_model.data = np.ones((3200, 2048)) * 15000  # bias for every pixel is 15000
    bias_model.meta.instrument.name = "NIRSPEC"
    bias_model.meta.instrument.detector = "NRS1"
    bias_model.meta.subarray.xstart = 1
    bias_model.meta.subarray.xsize = 2048
    bias_model.meta.subarray.ystart = 1
    bias_model.meta.subarray.ysize = 2048
    bias_model.meta.exposure.readpatt = "NRSIRS2"
    add_test_refmodel_metadata(bias_model)

    return data_model, saturation_model, bias_model


def setup_nrs_nrs_cube():
    """
    Set up mock NIRSpec data to test, with read pattern NRS.

    The output data shape is (1, 5, 3200, 2048).

    Returns
    -------
    data_model : `~stdatamodels.jwst.datamodels.RampModel`
        The NIRSpec IRS2 ramp model.
    saturation_model : `~stdatamodels.jwst.datamodels.SaturationModel`
        A corresponding saturation model.
    """
    # create a JWST datamodel for NIRSPEC NRS-mode data
    data_model = RampModel((1, 5, 3200, 2048))
    data_model.data = np.ones((1, 5, 3200, 2048))
    data_model.groupdq = np.zeros((1, 5, 3200, 2048))
    data_model.pixeldq = np.zeros((3200, 2048))
    data_model.meta.instrument.name = "NIRSPEC"
    data_model.meta.instrument.detector = "NRS1"
    data_model.meta.instrument.filter = "F100LP"
    data_model.meta.observation.date = "2019-07-19"
    data_model.meta.observation.time = "23:23:30.912"
    data_model.meta.exposure.type = "NRS_LAMP"
    data_model.meta.subarray.name = "FULL"
    data_model.meta.subarray.xstart = 1
    data_model.meta.subarray.xsize = 2048
    data_model.meta.subarray.ystart = 1
    data_model.meta.subarray.ysize = 2048
    data_model.meta.exposure.nrs_normal = 16
    data_model.meta.exposure.nrs_reference = 4
    data_model.meta.exposure.readpatt = "NRS"
    data_model.meta.exposure.nframes = 4
    data_model.meta.exposure.ngroups = 5

    # create a saturation model for the saturation step
    saturation_model = SaturationModel((2048, 2048))
    saturation_model.data = (
        np.ones((2048, 2048)) * 60000
    )  # saturation limit for every pixel is 60000
    saturation_model.meta.instrument.name = "NIRSPEC"
    saturation_model.meta.instrument.detector = "NRS1"
    saturation_model.meta.subarray.xstart = 1
    saturation_model.meta.subarray.xsize = 2048
    saturation_model.meta.subarray.ystart = 1
    saturation_model.meta.subarray.ysize = 2048
    add_test_refmodel_metadata(saturation_model)

    return data_model, saturation_model


def setup_nis_superstripe_cube():
    """
    Set up mock NIRISS superstripe data to test.

    The output data shape is (30, 5, 256, 208).

    Returns
    -------
    data_model : `~stdatamodels.jwst.datamodels.RampModel`
        The NIRSpec IRS2 ramp model.
    saturation_model : `~stdatamodels.jwst.datamodels.SaturationModel`
        A corresponding saturation model.
    bias_model : `~stdatamodels.jwst.datamodels.SuperbiasModel`
        A corresponding superbias model.
    """
    data_model = make_superstripe_model()
    num_stripes = data_model.meta.subarray.num_superstripe
    data_model.pixeldq = np.zeros((num_stripes, *data_model.data.shape[-2:]), dtype=np.uint32)

    saturation_model = SaturationModel((2048, 2048))
    saturation_model.meta.subarray.xstart = 1
    saturation_model.meta.subarray.ystart = 1
    saturation_model.meta.subarray.xsize = 2048
    saturation_model.meta.subarray.ysize = 2048
    saturation_model.meta.instrument.name = "NIRISS"
    add_test_refmodel_metadata(saturation_model)

    # create a bias model
    bias_model = SuperBiasModel((2028, 2048))
    bias_model.data = np.ones((2048, 2048)) * 15000  # bias for every pixel is 15000
    bias_model.meta.instrument.name = "NIRISS"
    bias_model.meta.instrument.detector = "NIS"
    bias_model.meta.subarray.xstart = 1
    bias_model.meta.subarray.xsize = 2048
    bias_model.meta.subarray.ystart = 1793
    bias_model.meta.subarray.ysize = 256
    bias_model.meta.subarray.fastaxis = -2
    bias_model.meta.subarray.slowaxis = -1
    bias_model.meta.exposure.readpatt = "ANY"
    add_test_refmodel_metadata(bias_model)

    return data_model, saturation_model, bias_model
