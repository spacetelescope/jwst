import math

import pytest
import numpy as np

from astropy import units as u

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import SpecModel, MultiSpecModel

from jwst.photom import photom

MJSR_TO_UJA2 = (u.megajansky / u.steradian).to(u.microjansky / (u.arcsecond**2))

# Multiply by this to convert from square arcseconds to steradians
A2_TO_SR = (np.pi / (180. * 3600.))**2


def mk_data(shape):
    """Utility function for creating test data.

    Parameters
    ----------
    shape : tuple
        Shape (either 2-D or 3-D) to use for creating data arrays.

    Returns
    -------
    tuple of ndarrays
        data is the science data array (float32).
        dq is the data quality array (uint32).
        err is the array of error estimates (float32).
        var_p is the contribution of Poisson noise to the variance (float32)
        var_r is the contribution of read noise to the variance (float32)
        var_f is the contribution of the flat-field to the variance (float32).
    """

    nelem = 1
    for k in shape:
        nelem *= k
    data = np.arange(1, nelem + 1, dtype=np.float32).reshape(shape)
    dq = np.zeros(shape, dtype=np.uint32)
    err = np.ones(shape, dtype=np.float32)
    var_p = np.ones(shape, dtype=np.float32)            # var_poisson
    var_r = np.ones(shape, dtype=np.float32)            # var_rnoise
    var_f = np.ones(shape, dtype=np.float32)            # var_flat

    return data, dq, err, var_p, var_r, var_f


def mk_wavelength(shape, min_wl, max_wl, dispaxis=1):
    """Create a 2-D array of wavelengths, linearly spaced in one axis.

    Parameters
    ----------
    shape : tuple
        Shape (either 2-D or 3-D) to use for creating data arrays.

    min_wl : float
        The minimum wavelength for the output array.

    max_wl : float
        The maximum wavelength for the output array.

    dispaxis : int
        Dispersion direction:  1 --> horizontal, 2 --> vertical

    Returns
    -------
    wl : 2-D ndarray
        The array of wavelengths.  The values will vary in the dispersion
        direction but be constant in the cross-dispersion direction.
    """

    # The wavelength attribute is always 2-D.
    if len(shape) > 2:
        shape = shape[-2:]
    wl = np.zeros(shape, dtype=np.float32)

    if dispaxis == 1:
        nx = shape[-1]
        x = np.linspace(min_wl, max_wl, nx, dtype=np.float32).reshape(1, nx)
        wl[:] = x.copy()
    elif dispaxis == 2:
        ny = shape[-2]
        y = np.linspace(min_wl, max_wl, ny, dtype=np.float32).reshape(ny, 1)
        wl[:] = y.copy()
    else:
        raise RuntimeError("dispaxis must be either 1 or 2.")

    return wl


def mk_soss_spec(settings, speclen):
    """Create a 2-D array of wavelengths, linearly spaced in one axis.

    Parameters
    ----------
    settings : list
        Each entry in the list should be a dict of settings to match
        to photom reference file entries - filter, pupil, and order.

    speclen : list
        Length of spectrum for each entry in settings; list lengths
        must match

    Returns
    -------
    model : MultiSpecModel
        The simulated output of extract_1d to be calibrated; this is
        the SOSS-specific ordering to be tested.
    """
    model = MultiSpecModel()
    for i, inspec in enumerate(settings):
        # Make number of columns equal to length of SpecModel's spec_table dtype, then assign
        # dtype to each column. Use to initialize SpecModel for entry into output MultiSpecModel
        otab = np.array(list(zip(*([np.linspace(0.6, 4.0, speclen[i])] +
                                   [np.ones(speclen[i]) for _ in range(len(SpecModel().spec_table.dtype) - 1)]))),
                        dtype=SpecModel().spec_table.dtype)
        specmodel = datamodels.SpecModel(spec_table=otab)
        model.meta.instrument.filter = inspec['filter']
        model.meta.instrument.pupil = inspec['pupil']
        specmodel.spectral_order = inspec['order']
        model.spec.append(specmodel)

    return model


def create_input(instrument, detector, exptype,
                 filter=None, pupil=None, grating=None, band=None):
    """Create dummy data (an open model) of the appropriate type.

    Parameters
    ----------
    instrument : str
        The instrument name (all upper case letters), one of:
        'NIRISS', 'NIRSPEC', 'NIRCAM', 'MIRI', 'FGS'.

    detector : str
        Detector name.  This is only used for populating a keyword.

    exptype : str
        Exposure type.  These can be explicitly checked for:
        'MIR_MRS', 'MIR_LRS-FIXEDSLIT',
        'NIS_SOSS', 'NIS_WFSS',
        'NRC_WFSS',
        'NRS_BRIGHTOBJ', 'NRS_FIXEDSLIT', 'NRS_IFU', 'NRS_MSASPEC'.

    filter : str or None
        Name of the element in the filter wheel.  For NIRISS WFSS, this
        is used to determine the dispersion direction.

    pupil : str or None
        Name of the element in the pupil wheel.  For NIRCam WFSS, this
        is used to determine the dispersion direction.

    grating : str or None
        Name of the element in the grating wheel.  This is only used for
        populating a keyword.

    band : str or None
        Band (MIRI only).  This is only used for populating a keyword.

    Returns
    -------
    input_model : `~jwst.datamodels.JwstDataModel`
        An open data model object of the appropriate type.
    """

    data = None  # Not defined for niriss_soss
    if instrument == 'NIRISS':
        if exptype == 'NIS_WFSS':
            nslits = 2
            input_model = datamodels.MultiSlitModel()
            if filter.endswith('R'):
                shape = (69, 5)
                dispaxis = 2                    # vertical
            else:
                shape = (5, 69)
                dispaxis = 1                    # horizontal
            input_model.meta.target.source_type = 'POINT'
            (data, dq, err, var_p, var_r, var_f) = mk_data(shape)
            wl = mk_wavelength(shape, 1.0, 5.0, dispaxis)
            for k in range(nslits):
                slit = datamodels.SlitModel(data=data, dq=dq, err=err,
                                            wavelength=wl)
                slit.var_poisson = var_p
                slit.var_rnoise = var_r
                slit.var_flat = var_f
                slit.meta.wcsinfo.spectral_order = k + 1
                # Not realistic, just something for a default.
                slit.meta.photometry.pixelarea_arcsecsq = 0.0025
                slit.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
                input_model.slits.append(slit.copy())
        elif exptype == 'NIS_SOSS':
            settings = [
                {'filter': filter, 'pupil': pupil, 'order': 1},
                {'filter': filter, 'pupil': pupil, 'order': 2},
            ]
            speclen = [200, 200]
            input_model = mk_soss_spec(settings, speclen)
        else:                                   # NIS_IMAGE
            shape = (96, 128)
            (data, dq, err, var_p, var_r, var_f) = mk_data(shape)
            input_model = datamodels.ImageModel(data=data, dq=dq, err=err)
            input_model.var_poisson = var_p
            input_model.var_rnoise = var_r
            input_model.var_flat = var_f
            input_model.meta.target.source_type = 'POINT'
            input_model.meta.photometry.pixelarea_arcsecsq = 0.0025
            input_model.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
    elif instrument == 'NIRSPEC':
        if exptype == 'NRS_FIXEDSLIT':
            nslits = 5
            input_model = datamodels.MultiSlitModel()
            shape = (5, 69)
            (data, dq, err, var_p, var_r, var_f) = mk_data(shape)
            wl = mk_wavelength(shape, 1.0, 5.0, dispaxis=1)
            input_model.meta.target.source_type = 'POINT'  # output will be MJy
            slitnames = ['S200A1', 'S200A2', 'S400A1', 'S1600A1', 'S200B1']
            srctypes = ['EXTENDED', 'POINT', 'EXTENDED', 'POINT', 'EXTENDED']
            for k in range(nslits):
                slit = datamodels.SlitModel(data=data, dq=dq, err=err,
                                            wavelength=wl)
                slit.name = slitnames[k]
                slit.source_type = srctypes[k]
                slit.var_poisson = var_p
                slit.var_rnoise = var_r
                slit.var_flat = var_f
                slit.meta.photometry.pixelarea_arcsecsq = 0.0025
                slit.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
                input_model.slits.append(slit.copy())
        elif exptype == 'NRS_BRIGHTOBJ':
            shape = (3, 5, 69)
            (data, dq, err, var_p, var_r, var_f) = mk_data(shape)
            wl = mk_wavelength(shape, 1.0, 5.0, dispaxis=1)
            input_model = datamodels.SlitModel(data=data, dq=dq, err=err,
                                               wavelength=wl)
            input_model.meta.target.source_type = 'POINT'  # output will be MJy
            input_model.var_poisson = var_p
            input_model.var_rnoise = var_r
            input_model.var_flat = var_f
            input_model.meta.photometry.pixelarea_arcsecsq = 0.0025
            input_model.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
            input_model.name = 'S1600A1'
        elif exptype == 'NRS_MSASPEC':
            nslits = 3
            input_model = datamodels.MultiSlitModel()
            shape = (5, 69)
            (data, dq, err, var_p, var_r, var_f) = mk_data(shape)
            wl = mk_wavelength(shape, 1.0, 5.0, dispaxis=1)
            for k in range(nslits):
                slit = datamodels.SlitModel(data=data, dq=dq, err=err,
                                            wavelength=wl)
                slit.name = str(k + 1)
                slit.var_poisson = var_p
                slit.var_rnoise = var_r
                slit.var_flat = var_f
                # This may be the only case where the source type is specified
                # in the slit attributes rather than in the primary header.
                slit.source_type = 'POINT'
                slit.meta.photometry.pixelarea_arcsecsq = 0.0025
                slit.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
                input_model.slits.append(slit.copy())
        else:
            # NRS_IFU needs the wcs, so we won't cover this case.  Use a
            # regression test instead.
            raise RuntimeError("exp_type {} is not currently tested"
                               .format(exptype))
    elif instrument == 'NIRCAM':
        if exptype == 'NRC_WFSS':
            nslits = 1
            input_model = datamodels.MultiSlitModel()
            if pupil.endswith('C'):
                shape = (69, 5)
                dispaxis = 2                    # vertical
            else:
                shape = (5, 69)
                dispaxis = 1                    # horizontal
            (data, dq, err, var_p, var_r, var_f) = mk_data(shape)
            wl = mk_wavelength(shape, 2.4, 5.0, dispaxis)
            input_model.meta.target.source_type = 'POINT'
            for k in range(nslits):
                slit = datamodels.SlitModel(data=data, dq=dq, err=err,
                                            wavelength=wl)
                slit.name = str(k + 1)
                slit.var_poisson = var_p
                slit.var_rnoise = var_r
                slit.var_flat = var_f
                slit.meta.wcsinfo.spectral_order = 1
                slit.meta.photometry.pixelarea_arcsecsq = 0.0025
                slit.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
                input_model.slits.append(slit.copy())
        else:                                   # NRC_IMAGE
            (data, dq, err, var_p, var_r, var_f) = mk_data((128, 256))
            input_model = datamodels.ImageModel(data=data, dq=dq, err=err)
            input_model.var_poisson = var_p
            input_model.var_rnoise = var_r
            input_model.var_flat = var_f
            input_model.meta.target.source_type = 'POINT'
            input_model.meta.photometry.pixelarea_arcsecsq = 0.0025
            input_model.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
    elif instrument == 'MIRI':
        if exptype == 'MIR_MRS':
            (data, dq, err, var_p, var_r, var_f) = mk_data((128, 256))
            input_model = datamodels.IFUImageModel(data=data, dq=dq, err=err)
            input_model.var_poisson = var_p
            input_model.var_rnoise = var_r
            input_model.var_flat = var_f
            input_model.meta.target.source_type = 'POINT'
            input_model.meta.photometry.pixelarea_arcsecsq = 0.0025
            input_model.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
        elif exptype == 'MIR_LRS-FIXEDSLIT':
            shape = (120, 100)
            array = np.zeros(shape, dtype=np.float32)
            data = np.arange(69 * 5, dtype=np.float32).reshape(69, 5)
            array[3:72, 15:20] = data
            dq = np.zeros(shape, dtype=np.uint32)
            err = np.ones(shape, dtype=np.float32)
            input_model = datamodels.ImageModel(data=array, dq=dq, err=err)
            # There is no wavelength attribute for ImageModel, but this
            # should work anyway.
            wl = mk_wavelength(shape, 5.0, 12.0, dispaxis=2)
            input_model.wavelength = wl.copy()
            input_model.var_poisson = np.ones(shape, dtype=np.float32)
            input_model.var_rnoise = np.ones(shape, dtype=np.float32)
            input_model.var_flat = np.ones(shape, dtype=np.float32)
            input_model.meta.subarray.name = 'FULL'
            input_model.meta.target.source_type = 'POINT'
            input_model.meta.photometry.pixelarea_arcsecsq = 0.0025
            input_model.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
        else:                                   # MIR_IMAGE
            shape = (128, 256)
            data = np.arange(128 * 256, dtype=np.float32).reshape(shape)
            dq = np.zeros(shape, dtype=np.uint32)
            err = np.ones(shape, dtype=np.float32)
            input_model = datamodels.ImageModel(data=data, dq=dq, err=err)
            input_model.var_poisson = np.ones(shape, dtype=np.float32)
            input_model.var_rnoise = np.ones(shape, dtype=np.float32)
            input_model.var_flat = np.ones(shape, dtype=np.float32)
            input_model.meta.subarray.name = 'SUB256'       # matches 'GENERIC'
            input_model.meta.target.source_type = 'POINT'
            input_model.meta.exposure.mid_time = 60000.0 # Added for new PHOTOM step
            input_model.meta.photometry.pixelarea_arcsecsq = 0.0025
            input_model.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
    elif instrument == 'FGS':
        shape = (64, 64)
        data = np.arange(64 * 64, dtype=np.float32).reshape(shape)
        dq = np.zeros(shape, dtype=np.uint32)
        err = np.ones(shape, dtype=np.float32)
        input_model = datamodels.ImageModel(data=data, dq=dq, err=err)
        input_model.var_poisson = np.ones(shape, dtype=np.float32)
        input_model.var_rnoise = np.ones(shape, dtype=np.float32)
        input_model.var_flat = np.ones(shape, dtype=np.float32)
        input_model.meta.target.source_type = 'POINT'
        input_model.meta.photometry.pixelarea_arcsecsq = 0.0025
        input_model.meta.photometry.pixelarea_steradians = 0.0025 * A2_TO_SR
    else:
        raise RuntimeError("instrument {} is not recognized".format(instrument))

    input_model.meta.instrument.name = instrument
    input_model.meta.instrument.detector = detector
    input_model.meta.exposure.type = exptype
    input_model.meta.subarray.xstart = 1
    input_model.meta.subarray.ystart = 1

    if data is not None:
        input_model.meta.subarray.xsize = data.shape[-1]
        input_model.meta.subarray.ysize = data.shape[-2]
    if filter is not None:
        input_model.meta.instrument.filter = filter
    if pupil is not None:
        input_model.meta.instrument.pupil = pupil
    if grating is not None:
        input_model.meta.instrument.grating = grating
    if band is not None:
        input_model.meta.instrument.band = band

    return input_model


def create_photom_nrs_fs(min_wl=1.0, max_wl=5.0, min_r=8.0, max_r=9.0):
    """Create a photom table for NIRSpec FS.

    Parameters
    ----------
    min_wl : float
        Minimum wavelength to assign when populating an array of wavelengths.

    max_wl : float
        Maximum wavelength to assign when populating an array of wavelengths.

    min_r : float
        Minimum value to assign when populating the relresponse array.

    max_r : float
        Maximum value to assign when populating the relresponse array.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a NIRSpec fixed-slit photom reference file.
    """

    filter = ["F100LP", "F100LP", "F100LP", "F100LP", "F100LP",
              "F100LP", "F100LP", "F100LP", "F100LP", "F100LP",
              "F170LP", "F170LP", "F170LP", "F170LP", "F170LP",
              "F170LP", "F170LP", "F170LP", "F170LP", "F170LP"]
    grating = ["G140M", "G140M", "G140M", "G140M", "G140M",
               "G235M", "G235M", "G235M", "G235M", "G235M",
               "G140M", "G140M", "G140M", "G140M", "G140M",
               "G235M", "G235M", "G235M", "G235M", "G235M"]
    slit = ["S200A1", "S200A2", "S400A1", "S1600A1", "S200B1",
            "S200A1", "S200A2", "S400A1", "S1600A1", "S200B1",
            "S200A1", "S200A2", "S400A1", "S1600A1", "S200B1",
            "S200A1", "S200A2", "S400A1", "S1600A1", "S200B1"]

    nrows = len(filter)
    nx = 3

    # [3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
    #  4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0]
    photmj = np.linspace(3.1, 3.1 + (nrows - 1.) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)
    nelem = np.zeros(nrows, dtype=np.int32) + nx
    x = np.linspace(min_wl, max_wl, nx, dtype=np.float32).reshape(1, nx)
    wavelength = np.zeros((nrows, nx), dtype=np.float32)
    wavelength[:] = x.copy()
    y = np.linspace(min_r, max_r, nx, dtype=np.float32).reshape(1, nx)
    relresponse = np.zeros((nrows, nx), dtype=np.float32)
    relresponse[:] = y.copy()
    reluncertainty = np.ones((nrows, nx), dtype=np.float32)

    nx = wavelength.shape[-1]

    dtype = np.dtype([('filter', 'S12'),
                      ('grating', 'S15'),
                      ('slit', 'S15'),
                      ('photmj', '<f4'),
                      ('uncertainty', '<f4'),
                      ('nelem', '<i2'),
                      ('wavelength', '<f4', (nx,)),
                      ('relresponse', '<f4', (nx,)),
                      ('reluncertainty', '<f4', (nx,))])
    reftab = np.array(list(zip(filter, grating, slit,
                               photmj, uncertainty, nelem,
                               wavelength, relresponse, reluncertainty)),
                      dtype=dtype)
    ftab = datamodels.NrsFsPhotomModel(phot_table=reftab)

    return ftab


def create_photom_nrs_msa(min_wl=1.0, max_wl=5.0, min_r=8.0, max_r=9.0):
    """Create a photom table for NIRSpec MSA.

    Parameters
    ----------
    min_wl : float
        Minimum wavelength to assign when populating an array of wavelengths.

    max_wl : float
        Maximum wavelength to assign when populating an array of wavelengths.

    min_r : float
        Minimum value to assign when populating the relresponse array.

    max_r : float
        Maximum value to assign when populating the relresponse array.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a NIRSpec MSA photom reference file.
    """

    filter = ["F100LP", "F100LP", "F170LP", "F170LP"]
    grating = ["G140M", "G235M", "G140M", "G235M"]

    nrows = len(filter)
    nx = 3

    # [3.1, 3.2, 3.3, 3.4]
    photmj = np.linspace(3.1, 3.1 + (nrows - 1.) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)
    nelem = np.zeros(nrows, dtype=np.int32) + nx
    x = np.linspace(min_wl, max_wl, nx, dtype=np.float32).reshape(1, nx)
    wavelength = np.zeros((nrows, nx), dtype=np.float32)
    wavelength[:] = x.copy()
    y = np.linspace(min_r, max_r, nx, dtype=np.float32).reshape(1, nx)
    relresponse = np.zeros((nrows, nx), dtype=np.float32)
    relresponse[:] = y.copy()
    reluncertainty = np.ones((nrows, nx), dtype=np.float32)

    dtype = np.dtype([('filter', 'S12'),
                      ('grating', 'S15'),
                      ('photmj', '<f4'),
                      ('uncertainty', '<f4'),
                      ('nelem', '<i2'),
                      ('wavelength', '<f4', (nx,)),
                      ('relresponse', '<f4', (nx,)),
                      ('reluncertainty', '<f4', (nx,))])
    reftab = np.array(list(zip(filter, grating,
                               photmj, uncertainty, nelem,
                               wavelength, relresponse, reluncertainty)),
                      dtype=dtype)
    ftab = datamodels.NrsMosPhotomModel(phot_table=reftab)

    return ftab


def create_photom_niriss_wfss(min_wl=1.0, max_wl=5.0, min_r=8.0, max_r=9.0):
    """Create a photom table for NIRISS WFSS.

    Parameters
    ----------
    min_wl : float
        Minimum wavelength to assign when populating an array of wavelengths.

    max_wl : float
        Maximum wavelength to assign when populating an array of wavelengths.

    min_r : float
        Minimum value to assign when populating the relresponse array.

    max_r : float
        Maximum value to assign when populating the relresponse array.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a NIRISS WFSS photom reference file.
    """

    filter = ["GR150C", "GR150C", "GR150C", "GR150C",
              "GR150R", "GR150R", "GR150R", "GR150R"]
    pupil = ["F140M", "F140M", "F200W", "F200W",
             "F140M", "F140M", "F200W", "F200W"]
    order = [1, 2, 1, 2, 1, 2, 1, 2]

    nrows = len(filter)
    nx = 3

    # [3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8]
    photmjsr = np.linspace(3.1, 3.1 + (nrows - 1.) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)
    nelem = np.zeros(nrows, dtype=np.int32) + nx
    x = np.linspace(min_wl, max_wl, nx, dtype=np.float32).reshape(1, nx)
    wavelength = np.zeros((nrows, nx), dtype=np.float32)
    wavelength[:] = x.copy()
    y = np.linspace(min_r, max_r, nx, dtype=np.float32).reshape(1, nx)
    relresponse = np.zeros((nrows, nx), dtype=np.float32)
    relresponse[:] = y.copy()
    reluncertainty = np.ones((nrows, nx), dtype=np.float32)

    dtype = np.dtype([('filter', 'S12'),
                      ('pupil', 'S15'),
                      ('order', '<i2'),
                      ('photmjsr', '<f4'),
                      ('uncertainty', '<f4'),
                      ('nelem', '<i2'),
                      ('wavelength', '<f4', (nx,)),
                      ('relresponse', '<f4', (nx,)),
                      ('reluncertainty', '<f4', (nx,))])
    reftab = np.array(list(zip(filter, pupil, order,
                               photmjsr, uncertainty, nelem,
                               wavelength, relresponse, reluncertainty)),
                      dtype=dtype)
    ftab = datamodels.NisWfssPhotomModel(phot_table=reftab)

    return ftab


def create_photom_niriss_soss(min_r=8.0, max_r=9.0):
    """Create a photom table for NIRISS SOSS.

    Parameters
    ----------
    min_r : float
        Minimum value to assign when populating the relresponse array.

    max_r : float
        Maximum value to assign when populating the relresponse array.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a NIRISS SOSS photom reference file.
    """

    filter = ["CLEAR", "CLEAR"]
    pupil = ["GR700XD", "GR700XD"]
    order = [1, 2]

    nrows = len(filter)
    nx = 3

    # [3.1, 3.2]
    photmj = np.linspace(3.1, 3.1 + (nrows - 1.) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)
    nelem = np.zeros(nrows, dtype=np.int32) + nx
    wavelength = np.zeros((nrows, nx), dtype=np.float32)
    wavelength[0, :] = np.linspace(0.9, 2.8, nx, dtype=np.float32)
    wavelength[1, :] = np.linspace(0.6, 1.4, nx, dtype=np.float32)
    y = np.linspace(min_r, max_r, nx, dtype=np.float32).reshape(1, nx)
    relresponse = np.zeros((nrows, nx), dtype=np.float32)
    relresponse[:] = y.copy()
    reluncertainty = np.ones((nrows, nx), dtype=np.float32)

    dtype = np.dtype([('filter', 'S12'),
                      ('pupil', 'S12'),
                      ('order', '<i2'),
                      ('photmj', '<f4'),
                      ('uncertainty', '<f4'),
                      ('nelem', '<i2'),
                      ('wavelength', '<f4', (nx,)),
                      ('relresponse', '<f4', (nx,)),
                      ('reluncertainty', '<f4', (nx,))])
    reftab = np.array(list(zip(filter, pupil, order,
                               photmj, uncertainty, nelem,
                               wavelength, relresponse, reluncertainty)),
                      dtype=dtype)
    ftab = datamodels.NisSossPhotomModel(phot_table=reftab)

    return ftab


def create_photom_niriss_image(min_r=8.0, max_r=9.0):
    """Create a photom table for NIRISS image.

    Parameters
    ----------
    min_r : float
        Minimum value to assign when populating the relresponse array.

    max_r : float
        Maximum value to assign when populating the relresponse array.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a NIRISS image photom reference file.
    """

    # The middle row should be selected.
    filter = ["F430M", "CLEAR", "CLEAR"]
    pupil = ["F090W", "F140M", "F140W"]

    nrows = len(filter)

    # [3.1, 3.2, 3.3]
    photmjsr = np.linspace(3.1, 3.1 + (nrows - 1.) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)

    dtype = np.dtype([('filter', 'S12'),
                      ('pupil', 'S12'),
                      ('photmjsr', '<f4'),
                      ('uncertainty', '<f4')])
    reftab = np.array(list(zip(filter, pupil, photmjsr, uncertainty)),
                      dtype=dtype)
    ftab = datamodels.NisImgPhotomModel(phot_table=reftab)

    return ftab


def create_photom_miri_mrs(shape, value, pixel_area, photmjsr):
    """Create a photom reference file for MIRI MRS.

    Parameters
    ----------
    shape : tuple
        The shape to use when creating image arrays.

    value : float
        The value to assign to the SCI data array.

    pixel_area : float
        The pixel solid angle in steradians.

    photmjsr : float
        The value to assign to the MJy / sr keyword.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a MIRI MRS photom reference file.
    """

    data = np.zeros(shape, dtype=np.float32) + value
    err = np.ones(shape, dtype=np.float32)
    dq = np.zeros(shape, dtype=np.uint32)
    pixsiz = np.zeros(shape, dtype=np.float32) + pixel_area

    ftab = datamodels.MirMrsPhotomModel(data=data, err=err, dq=dq,
                                        pixsiz=pixsiz)
    ftab.meta.photometry.conversion_megajanskys = photmjsr

    return ftab


def create_photom_miri_lrs(min_wl=5.0, max_wl=10.0, min_r=8.0, max_r=9.0):
    """Create a photom table for MIRI LRS.

    Parameters
    ----------
    min_wl : float
        Minimum wavelength to assign when populating an array of wavelengths.

    max_wl : float
        Maximum wavelength to assign when populating an array of wavelengths.

    min_r : float
        Minimum value to assign when populating the relresponse array.

    max_r : float
        Maximum value to assign when populating the relresponse array.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a MIRI LRS photom reference file.
    """

    filter = ["F560W", "P750L", "F1000W"]
    subarray = ["GENERIC", "FULL", "GENERIC"]

    nrows = len(filter)
    nx = 3

    # [3.1, 3.2, 3.3]
    photmjsr = np.linspace(3.1, 3.1 + (nrows - 1.) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)
    nelem = np.zeros(nrows, dtype=np.int32) + nx
    x = np.linspace(min_wl, max_wl, nx, dtype=np.float32).reshape(1, nx)
    wavelength = np.zeros((nrows, nx), dtype=np.float32)
    wavelength[:] = x.copy()
    y = np.linspace(min_r, max_r, nx, dtype=np.float32).reshape(1, nx)
    relresponse = np.zeros((nrows, nx), dtype=np.float32)
    relresponse[:] = y.copy()
    reluncertainty = np.ones((nrows, nx), dtype=np.float32)

    dtype = np.dtype([('filter', 'S12'),
                      ('subarray', 'S15'),
                      ('photmjsr', '<f4'),
                      ('uncertainty', '<f4'),
                      ('nelem', '<i2'),
                      ('wavelength', '<f4', (nx,)),
                      ('relresponse', '<f4', (nx,)),
                      ('reluncertainty', '<f4', (nx,))])
    reftab = np.array(list(zip(filter, subarray,
                               photmjsr, uncertainty, nelem,
                               wavelength, relresponse, reluncertainty)),
                      dtype=dtype)
    ftab = datamodels.MirLrsPhotomModel(phot_table=reftab)

    return ftab


def create_photom_miri_image(min_wl=16.5, max_wl=19.5,
                             min_r=8.0, max_r=9.0):
    """Create a photom table for MIRI image mode.

    Parameters
    ----------
    min_wl : float
        Minimum wavelength to assign when populating an array of wavelengths.

    max_wl : float
        Maximum wavelength to assign when populating an array of wavelengths.

    min_r : float
        Minimum value to assign when populating the relresponse array.

    max_r : float
        Maximum value to assign when populating the relresponse array.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a MIRI image photom reference file.
    """

    filter = ["F1800W", "F2100W", "F2550W"]
    subarray = ["SUB256", "SUB256", "SUB256"]

    nrows = len(filter)

    # [3.1, 3.2, 3.3]
    photmjsr = np.linspace(3.1, 3.1 + (nrows - 1.) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)

    dtype = np.dtype([('filter', 'S12'),
                      ('subarray', 'S15'),
                      ('photmjsr', '<f4'),
                      ('uncertainty', '<f4')])
    reftab = np.array(list(zip(filter, subarray, photmjsr, uncertainty)),
                      dtype=dtype)
    timecoeff_amp = np.linspace(2.1, 2.1 + (nrows - 1.) * 0.1, nrows)
    timecoeff_tau = np.asarray([145, 145, 145])
    timecoeff_t0 = np.asarray([59720, 59720, 59720])
    dtypec = np.dtype([('amplitude', '<f4'),
                      ('tau', '<f4'),
                      ('t0', '<f4')])
    reftabc = np.array(list(zip(timecoeff_amp, timecoeff_tau, timecoeff_t0)),
                      dtype=dtypec)
    ftab = datamodels.MirImgPhotomModel(phot_table=reftab, timecoeff=reftabc)

    return ftab


def create_photom_nircam_image(min_r=8.0, max_r=9.0):
    """Create a photom table for NIRCam image.

    Parameters
    ----------
    min_r : float
        Minimum value to assign when populating the relresponse array.

    max_r : float
        Maximum value to assign when populating the relresponse array.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a NIRCam image photom reference file.
    """

    filter = ["F090W", "F115W", "F150W", "F200W"]
    pupil = ["F162M", "F164N", "CLEAR", "WLP8"]

    nrows = len(filter)

    # [3.1, 3.2, 3.3, 3.4]
    photmjsr = np.linspace(3.1, 3.1 + (nrows - 1.) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)

    dtype = np.dtype([('filter', 'S12'),
                      ('pupil', 'S12'),
                      ('photmjsr', '<f4'),
                      ('uncertainty', '<f4')])
    reftab = np.array(list(zip(filter, pupil, photmjsr, uncertainty)),
                      dtype=dtype)

    ftab = datamodels.NrcImgPhotomModel(phot_table=reftab)

    return ftab


def create_photom_nircam_wfss(min_wl=2.4, max_wl=5.0, min_r=8.0, max_r=9.0):
    """Create a photom table for NIRCam WFSS.

    Parameters
    ----------
    min_wl : float
        Minimum wavelength to assign when populating an array of wavelengths.

    max_wl : float
        Maximum wavelength to assign when populating an array of wavelengths.

    min_r : float
        Minimum value to assign when populating the relresponse array.

    max_r : float
        Maximum value to assign when populating the relresponse array.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a NIRCam WFSS photom reference file.
    """

    filter = ["F277W", "F322W2", "F356W", "F410M", "F444W"]
    pupil = ["GRISMR", "GRISMR", "GRISMR", "GRISMR", "GRISMR"]
    order = [1, 1, 1, 1, 1]

    nrows = len(filter)
    nx = 3

    # [3.1, 3.2, 3.3, 3.4, 3.5]
    photmjsr = np.linspace(3.1, 3.1 + (nrows - 1.) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)
    nelem = np.zeros(nrows, dtype=np.int32) + nx
    x = np.linspace(min_wl, max_wl, nx, dtype=np.float32).reshape(1, nx)
    wavelength = np.zeros((nrows, nx), dtype=np.float32)
    wavelength[:] = x.copy()
    y = np.linspace(min_r, max_r, nx, dtype=np.float32).reshape(1, nx)
    relresponse = np.zeros((nrows, nx), dtype=np.float32)
    relresponse[:] = y.copy()
    relunc = np.zeros((nrows, nx), dtype=np.float32)

    dtype = np.dtype([('filter', 'S12'),
                      ('pupil', 'S15'),
                      ('order', '<i2'),
                      ('photmjsr', '<f4'),
                      ('uncertainty', '<f4'),
                      ('nelem', '<i2'),
                      ('wavelength', '<f4', (nx,)),
                      ('relresponse', '<f4', (nx,)),
                      ('reluncertainty', '<f4', (nx,))])
    reftab = np.array(list(zip(filter, pupil, order, photmjsr, uncertainty,
                               nelem, wavelength, relresponse, relunc)),
                      dtype=dtype)

    ftab = datamodels.NrcWfssPhotomModel(phot_table=reftab)

    return ftab


def create_photom_fgs_image(value):
    """Create a photom table for FGS.

    Parameters
    ----------
    value : float
        The value to assign to the MJy / sr column.

    Returns
    -------
    ftab : `~jwst.datamodels.JwstDataModel`
        An open data model for a NIRSpec fixed-slit photom reference file.
    """

    photmjsr = [value]
    uncertainty = [0.0]

    dtype = np.dtype([('photmjsr', '<f4'),
                      ('uncertainty', '<f4')])
    reftab = np.array(list(zip(photmjsr, uncertainty)),
                      dtype=dtype)
    ftab = datamodels.FgsImgPhotomModel(phot_table=reftab)

    return ftab


def create_pixel_area_ref(shape, area_ster, area_a2):
    """Create a pixel area (solid angle) reference file.

    Parameters
    ----------
    shape : tuple
        Shape to use when creating the data array.

    area_ster : float
        Pixel area in units of steradians.

    area_a2 : float
        Pixel area in units of arcseconds squared.

    Returns
    -------
    area_ref : `~jwst.datamodels.JwstDataModel`
        An open data model for a pixel area reference file.
    """

    data = np.ones(shape, dtype=np.float32)
    area_ref = datamodels.PixelAreaModel(data=data)
    area_ref.meta.photometry.pixelarea_steradians = area_ster
    area_ref.meta.photometry.pixelarea_arcsecsq = area_a2

    return area_ref


def create_msa_pixel_area_ref(quadrant, shutter_x, shutter_y, pixarea):
    """Create a pixel area (solid angle) reference file for NIRSpec MSA.

    Parameters
    ----------
    quadrant : ndarray, int16
        Array of MOS quadrant indices.

    shutter_x : ndarray, int16
        Array of X shutter locations.

    shutter_y : ndarray, int16
        Array of Y shutter locations.

    pixarea : ndarray, float32
        Array of pixel area values (arcsec^2).

    Returns
    -------
    area_ref : `~jwst.datamodels.JwstDataModel`
        An open data model for a pixel area reference file.
    """

    dtype = np.dtype([('quadrant', '<i2'),
                      ('shutter_x', '<i2'),
                      ('shutter_y', '<i2'),
                      ('pixarea', '<f4')])
    reftab = np.array(list(zip(quadrant, shutter_x, shutter_y, pixarea)),
                      dtype=dtype)
    area_ref = datamodels.NirspecMosAreaModel(area_table=reftab)

    return area_ref


def find_row_in_ftab(input_model, ftab, select, slitname=None, order=None):
    """Find the matching row in the photom reference file.

    Parameters
    ----------
    input_model : `~jwst.datamodels.JwstDataModel`
        input Data Model object

    ftab : `~jwst.datamodels.JwstDataModel`
        This has a `phot_table` attribute, which is a table containing
        photometric information.  This can be any of several data models
        functionally equivalent to _photom_xxxx.fits files in CRDS.

    select : list of str
        The strings in this list can be any of 'filter', 'grating', or
        'pupil'.  For each of these, the value to search for in
        `ftab.phot_table` will be gotten from the metadata for `input_model`.
        `slitname` and `order` are not included in this list because their
        values are not obtained from the metadata.

    slitname : str or None
        The name of the slit.  This will be included in the selection
        criteria if it is not None.

    order : int or None
        The spectral order number.  This will be included in the selection
        criteria if it is not None.

    Returns
    -------
    rownum : int
        The zero-based number of the row in `ftab` that matches the
        selection criteria.

    Raises
    ------
    RuntimeError
        If no matching row is found in `ftab`.
    """

    if 'filter' in select:
        filter = input_model.meta.instrument.filter
        filter_c = ftab.phot_table['filter']
    else:
        filter = None

    if 'grating' in select:
        grating = input_model.meta.instrument.grating
        grating_c = ftab.phot_table['grating']
    else:
        grating = None

    if 'pupil' in select:
        pupil = input_model.meta.instrument.pupil
        pupil_c = ftab.phot_table['pupil']
    else:
        pupil = None

    if slitname is not None:
        slitname_c = ftab.phot_table['slit']

    if order is not None:
        order_c = ftab.phot_table['order']

    nrows = len(ftab.phot_table)
    foundit = False
    for rownum in range(nrows):
        if ((filter is None or filter == filter_c[rownum])
                and (grating is None or grating == grating_c[rownum])
                and (pupil is None or pupil == pupil_c[rownum])
                and (slitname is None or slitname == slitname_c[rownum])
                and (order is None or order == order_c[rownum])):
            foundit = True
            break
    if not foundit:
        raise RuntimeError("Row not found in ftab.")

    return rownum


def test_nirspec_fs():
    """Test calc_nirspec, fixed-slit data"""

    input_model = create_input('NIRSPEC', 'NRS1', 'NRS_FIXEDSLIT',
                               filter='F170LP', grating='G235M')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)

    ftab = create_photom_nrs_fs(min_wl=1.0, max_wl=5.0,
                                min_r=8.0, max_r=9.0)
    ds.calc_nirspec(ftab, 'garbage')    # area_fname isn't used for this mode

    result = []
    for (k, slit) in enumerate(save_input.slits):
        slitname = slit.name
        srctype = slit.source_type
        input = slit.data                       # this is from save_input
        output = ds.input.slits[k].data         # ds.input is the output
        rownum = find_row_in_ftab(save_input, ftab, ['filter', 'grating'],
                                  slitname, order=None)
        photmj = ftab.phot_table['photmj'][rownum]
        nelem = ftab.phot_table['nelem'][rownum]
        wavelength = ftab.phot_table['wavelength'][rownum][0:nelem]
        relresponse = ftab.phot_table['relresponse'][rownum][0:nelem]
        shape = input.shape
        ix = shape[1] // 2
        iy = shape[0] // 2
        wl = slit.wavelength[iy, ix]
        rel_resp = np.interp(wl, wavelength, relresponse,
                             left=np.nan, right=np.nan)
        compare = photmj * rel_resp
        if srctype != 'POINT':
            compare /= slit.meta.photometry.pixelarea_steradians

        # Compare the values at the center pixel.
        ratio = output[iy, ix] / input[iy, ix]
        result.append(np.allclose(ratio, compare, rtol=1.e-7))

        # Check error array and variance arrays.  This doesn't need to be
        # done for every instrument, because the calc_xxx functions all
        # call photom_io, and that's the function that modifies the data,
        # err, and variance arrays.  This does need to be checked for two
        # cases, though:  once for MultiSlitModel and once for everything
        # else, because those two cases are handled in two separate sections
        # of photom_io.
        ratio_err = ds.input.slits[k].err[iy, ix] / slit.err[iy, ix]
        result.append(np.allclose(ratio_err, compare, rtol=1.e-7))
        ratio_var_p = np.sqrt(ds.input.slits[k].var_poisson[iy, ix] /
                              slit.var_poisson[iy, ix])
        result.append(np.allclose(ratio_var_p, compare, rtol=1.e-7))
        ratio_var_r = np.sqrt(ds.input.slits[k].var_rnoise[iy, ix] /
                              slit.var_rnoise[iy, ix])
        result.append(np.allclose(ratio_var_r, compare, rtol=1.e-7))
        ratio_var_f = np.sqrt(ds.input.slits[k].var_flat[iy, ix] /
                              slit.var_flat[iy, ix])
        result.append(np.allclose(ratio_var_f, compare, rtol=1.e-7))
        # The output units are flux density for point sources and
        # surface brightness for extended.
        if srctype == 'POINT':
            result.append(ds.input.slits[k].meta.bunit_data == 'MJy')
            result.append(ds.input.slits[k].meta.bunit_err == 'MJy')
        else:
            result.append(ds.input.slits[k].meta.bunit_data == 'MJy/sr')
            result.append(ds.input.slits[k].meta.bunit_err == 'MJy/sr')

    assert np.all(result)

    ftab.close()


def test_nirspec_bright():
    """Test calc_nirspec, bright-object data"""

    input_model = create_input('NIRSPEC', 'NRS1', 'NRS_BRIGHTOBJ',
                               filter='F170LP', grating='G235M')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)

    # The FS photom table can be used for BRIGHTOBJ as well.
    ftab = create_photom_nrs_fs(min_wl=1.0, max_wl=5.0,
                                min_r=8.0, max_r=9.0)
    ds.calc_nirspec(ftab, 'garbage')

    shape = input_model.data.shape

    slitname = 'S1600A1'                        # for brightobj mode

    input = save_input.data
    output = ds.input.data
    rownum = find_row_in_ftab(save_input, ftab, ['filter', 'grating'],
                              slitname, order=None)
    photmj = ftab.phot_table['photmj'][rownum]
    nelem = ftab.phot_table['nelem'][rownum]
    wavelength = ftab.phot_table['wavelength'][rownum][0:nelem]
    relresponse = ftab.phot_table['relresponse'][rownum][0:nelem]
    ix = shape[-1] // 2
    iy = shape[-2] // 2
    wl = save_input.wavelength[iy, ix]
    rel_resp = np.interp(wl, wavelength, relresponse,
                         left=np.nan, right=np.nan)

    compare = photmj * rel_resp
    ratio = output[:, iy, ix] / input[:, iy, ix]
    result = []
    result.append(np.allclose(ratio, compare, rtol=1.e-7))

    # Check error array and variance arrays.  This doesn't need to be
    # done for every instrument, because the calc_xxx functions all call
    # photom_io, and that's the function that modifies the data, err, and
    # variance arrays.  This does need to be checked for two cases, though:
    # once for MultiSlitModel and once for everything else, because those
    # two cases are handled in two separate sections of photom_io.
    ratio_err = ds.input.err[:, iy, ix] / save_input.err[:, iy, ix]
    result.append(np.allclose(ratio_err, compare, rtol=1.e-7))
    ratio_var_p = np.sqrt(ds.input.var_poisson[:, iy, ix] /
                          save_input.var_poisson[:, iy, ix])
    result.append(np.allclose(ratio_var_p, compare, rtol=1.e-7))
    ratio_var_r = np.sqrt(ds.input.var_rnoise[:, iy, ix] /
                          save_input.var_rnoise[:, iy, ix])
    result.append(np.allclose(ratio_var_r, compare, rtol=1.e-7))
    ratio_var_f = np.sqrt(ds.input.var_flat[:, iy, ix] /
                          save_input.var_flat[:, iy, ix])
    result.append(np.allclose(ratio_var_f, compare, rtol=1.e-7))
    result.append(ds.input.meta.bunit_data == 'MJy')
    result.append(ds.input.meta.bunit_err == 'MJy')

    assert np.all(result)


def test_nirspec_msa():
    """Test calc_nirspec, MSA data"""

    input_model = create_input('NIRSPEC', 'NRS1', 'NRS_MSASPEC',
                               filter='F170LP', grating='G235M')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)

    ftab = create_photom_nrs_msa(min_wl=1.0, max_wl=5.0,
                                 min_r=8.0, max_r=9.0)
    ds.calc_nirspec(ftab, 'garbage')

    # xxx The slit name is currently not used by photom for MSA data, but
    # it probably will be used at some time in the future.
    rownum = find_row_in_ftab(save_input, ftab, ['filter', 'grating'],
                              slitname=None, order=None)
    photmj = ftab.phot_table['photmj'][rownum]
    nelem = ftab.phot_table['nelem'][rownum]
    wavelength = ftab.phot_table['wavelength'][rownum][0:nelem]
    relresponse = ftab.phot_table['relresponse'][rownum][0:nelem]

    result = []
    for (k, slit) in enumerate(save_input.slits):
        input = slit.data                       # this is from save_input
        output = ds.input.slits[k].data         # ds.input is the output

        shape = input.shape
        ix = shape[1] // 2
        iy = shape[0] // 2
        wl = slit.wavelength[iy, ix]
        rel_resp = np.interp(wl, wavelength, relresponse,
                             left=np.nan, right=np.nan)
        compare = photmj * rel_resp

        ratio = output[iy, ix] / input[iy, ix]
        result.append(np.allclose(ratio, compare, rtol=1.e-7))

    assert np.all(result)


""" Skip this test because it would require a realistic wcs.
def test_nirspec_ifu():

    input_model = create_input('NIRSPEC', 'NRS1', 'NRS_IFU',
                               filter='F170LP', grating='G235M')
    ds = photom.DataSet(input_model)
"""


def test_niriss_wfss():
    """Test calc_niriss, WFSS data"""

    input_model = create_input('NIRISS', 'NIS', 'NIS_WFSS',
                               filter='GR150R', pupil='F140M')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)
    ftab = create_photom_niriss_wfss(min_wl=1.0, max_wl=5.0,
                                     min_r=8.0, max_r=9.0)
    ds.calc_niriss(ftab)

    result = []
    for (k, slit) in enumerate(save_input.slits):
        input = slit.data                       # this is from save_input
        output = ds.input.slits[k].data         # ds.input is the output
        sp_order = slit.meta.wcsinfo.spectral_order
        rownum = find_row_in_ftab(save_input, ftab, ['filter', 'pupil'],
                                  slitname=None, order=sp_order)
        photmjsr = ftab.phot_table['photmjsr'][rownum]
        nelem = ftab.phot_table['nelem'][rownum]
        wavelength = ftab.phot_table['wavelength'][rownum][0:nelem]
        relresponse = ftab.phot_table['relresponse'][rownum][0:nelem]
        shape = input.shape
        ix = shape[1] // 2
        iy = shape[0] // 2
        wl = slit.wavelength[iy, ix]
        rel_resp = np.interp(wl, wavelength, relresponse,
                             left=np.nan, right=np.nan)
        compare = photmjsr * rel_resp
        # Compare the values at the center pixel.
        ratio = output[iy, ix] / input[iy, ix]
        result.append(np.allclose(ratio, compare, rtol=1.e-7))

    assert np.all(result)


def test_niriss_soss():
    """Test calc_niriss, SOSS data"""

    input_model = create_input('NIRISS', 'NIS', 'NIS_SOSS',
                               filter='CLEAR', pupil='GR700XD')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)
    ftab = create_photom_niriss_soss(min_r=8.0, max_r=9.0)
    ds.calc_niriss(ftab)

    input = save_input.spec[0].spec_table['FLUX']
    output = ds.input.spec[0].spec_table['FLUX']                      # ds.input is the output
    sp_order = 1                                # to agree with photom.py
    rownum = find_row_in_ftab(save_input, ftab, ['filter', 'pupil'],
                              slitname=None, order=sp_order)
    photmj = ftab.phot_table['photmj'][rownum]
    nelem = ftab.phot_table['nelem'][rownum]
    wavelength = ftab.phot_table['wavelength'][rownum][0:nelem]
    relresponse = ftab.phot_table['relresponse'][rownum][0:nelem]
    test_ind = len(input) // 2
    wl = input_model.spec[0].spec_table['WAVELENGTH'][test_ind]
    rel_resp = np.interp(wl, wavelength, relresponse,
                         left=np.nan, right=np.nan)
    compare = photmj * rel_resp
    # Compare the values at the center pixel.
    ratio = output[test_ind] / input[test_ind]
    assert np.allclose(ratio, compare, rtol=1.e-7)


def test_niriss_image():
    """Test calc_niriss, image data"""

    input_model = create_input('NIRISS', 'NIS', 'NIS_IMAGE',
                               filter='CLEAR', pupil='F140M')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)
    ftab = create_photom_niriss_image(min_r=8.0, max_r=9.0)
    ds.calc_niriss(ftab)

    input = save_input.data
    output = ds.input.data                      # ds.input is the output
    rownum = find_row_in_ftab(save_input, ftab, ['filter', 'pupil'],
                              slitname=None)
    photmjsr = ftab.phot_table['photmjsr'][rownum]
    shape = input.shape
    ix = shape[1] // 2
    iy = shape[0] // 2
    compare = photmjsr
    # Compare the values at the center pixel.
    ratio = output[iy, ix] / input[iy, ix]
    assert np.allclose(ratio, compare, rtol=1.e-7)


def test_miri_mrs():
    """Test calc_miri, MRS data"""

    input_model = create_input('MIRI', 'MIRIFULONG', 'MIR_MRS',
                               filter='F1500W', band='LONG')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)
    value = 1.436
    pixel_area = 0.0436
    photmjsr = 17.3
    shape = save_input.data.shape
    ftab = create_photom_miri_mrs(shape, value=value,
                                  pixel_area=pixel_area, photmjsr=photmjsr)
    ds.calc_miri(ftab)

    input = save_input.data
    output = ds.input.data                      # ds.input is the output
    ix = shape[1] // 2
    iy = shape[0] // 2

    result = []
    # Check the photometry keywords.
    result.append(math.isclose(photmjsr,
                               ds.input.meta.photometry.conversion_megajanskys,
                               rel_tol=1.e-12))
    result.append(math.isclose(photmjsr * MJSR_TO_UJA2,
                               ds.input.meta.photometry.conversion_microjanskys,
                               rel_tol=1.e-12))
    # Check the data values.
    compare = value
    ratio = output[iy, ix] / input[iy, ix]
    result.append(math.isclose(ratio, compare, rel_tol=1.e-7))
    assert np.all(result)


def test_miri_lrs():
    """Test calc_miri, LRS data"""

    input_model = create_input('MIRI', 'MIRIMAGE', 'MIR_LRS-FIXEDSLIT',
                               filter='P750L')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)
    ftab = create_photom_miri_lrs(min_wl=5.0, max_wl=10.0,
                                  min_r=8.0, max_r=9.0)
    ds.calc_miri(ftab)

    input = save_input.data
    output = ds.input.data                      # ds.input is the output
    # Actual row selection can also require a match with SUBARRAY.
    rownum = find_row_in_ftab(save_input, ftab, ['filter'],
                              slitname=None, order=None)
    photmjsr = ftab.phot_table['photmjsr'][rownum]
    nelem = ftab.phot_table['nelem'][rownum]
    wavelength = ftab.phot_table['wavelength'][rownum][0:nelem]
    relresponse = ftab.phot_table['relresponse'][rownum][0:nelem]
    # see `array[3:72, 15:20] = data` in exptype == 'MIR_LRS-FIXEDSLIT'
    ix = 17
    iy = 37
    wl = input_model.wavelength[iy, ix]
    rel_resp = np.interp(wl, wavelength, relresponse,
                         left=np.nan, right=np.nan)
    compare = photmjsr * rel_resp
    # Compare the values at the center pixel.
    ratio = output[iy, ix] / input[iy, ix]
    assert np.allclose(ratio, compare, rtol=1.e-7)


def test_miri_image():
    """Test calc_miri, image data"""

    input_model = create_input('MIRI', 'MIRIMAGE', 'MIR_IMAGE',
                               filter='F1800W')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)
    ftab = create_photom_miri_image(min_wl=16.5, max_wl=19.5,
                                    min_r=8.0, max_r=9.0)
    ds.calc_miri(ftab)

    input = save_input.data
    output = ds.input.data                      # ds.input is the output
    rownum = find_row_in_ftab(save_input, ftab, ['filter'],
                              slitname=None, order=None)
    photmjsr = ftab.phot_table['photmjsr'][rownum]
    amplitude = ftab.timecoeff['amplitude'][rownum]
    tau = ftab.timecoeff['tau'][rownum]
    t0 = ftab.timecoeff['t0'][rownum]
    shape = input.shape
    ix = shape[1] // 2
    iy = shape[0] // 2
    compare = photmjsr + amplitude * np.exp(-(60000 - t0)/tau) # Added for new PHOTOM step
    # Compare the values at the center pixel.
    ratio = output[iy, ix] / input[iy, ix]
    assert np.allclose(ratio, compare, rtol=1.e-7)


def test_nircam_image():
    """Test calc_nircam, image data"""

    input_model = create_input('NIRCAM', 'NRCA3', 'NRC_IMAGE',
                               filter='F150W', pupil='CLEAR')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)
    ftab = create_photom_nircam_image(min_r=8.0, max_r=9.0)
    ds.calc_nircam(ftab)

    input = save_input.data
    output = ds.input.data                      # ds.input is the output
    rownum = find_row_in_ftab(save_input, ftab, ['filter', 'pupil'],
                              slitname=None, order=None)
    photmjsr = ftab.phot_table['photmjsr'][rownum]
    shape = input.shape
    ix = shape[1] // 2
    iy = shape[0] // 2
    compare = photmjsr
    # Compare the values at the center pixel.
    ratio = output[iy, ix] / input[iy, ix]
    assert np.allclose(ratio, compare, rtol=1.e-7)


def test_nircam_spec():
    """Test calc_nircam, WFSS data"""

    input_model = create_input('NIRCAM', 'NRCALONG', 'NRC_WFSS',
                               filter='F356W', pupil='GRISMR')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)
    ftab = create_photom_nircam_wfss(min_wl=2.4, max_wl=5.0,
                                     min_r=8.0, max_r=9.0)
    ds.calc_nircam(ftab)

    for (k, slit) in enumerate(save_input.slits):

        input = slit.data
        output = ds.input.slits[k].data         # ds.input is the output
        rownum = find_row_in_ftab(save_input, ftab, ['filter', 'pupil'],
                                  slitname=None, order=None)
        photmjsr = ftab.phot_table['photmjsr'][rownum]
        shape = input.shape
        ix = shape[1] // 2
        iy = shape[0] // 2
        nelem = ftab.phot_table['nelem'][rownum]
        wavelength = ftab.phot_table['wavelength'][rownum][0:nelem]
        relresponse = ftab.phot_table['relresponse'][rownum][0:nelem]
        shape = input.shape
        ix = shape[1] // 2
        iy = shape[0] // 2
        wl = slit.wavelength[iy, ix]
        rel_resp = np.interp(wl, wavelength, relresponse,
                             left=np.nan, right=np.nan)
        compare = photmjsr * rel_resp
        # Compare the values at the center pixel.
        ratio = output[iy, ix] / input[iy, ix]
        assert np.allclose(ratio, compare, rtol=1.e-7)


def test_fgs():
    """Test calc_fgs"""

    input_model = create_input('FGS', 'GUIDER1', 'FGS_IMAGE')
    save_input = input_model.copy()
    ds = photom.DataSet(input_model)
    value = 3.9
    ftab = create_photom_fgs_image(value)
    ds.calc_fgs(ftab)

    input = save_input.data
    output = ds.input.data                      # ds.input is the output
    # The FGS reference file has only one row, and there is no selection
    # criterion.
    rownum = 0
    photmjsr = ftab.phot_table['photmjsr'][rownum]
    shape = input.shape
    ix = shape[1] // 2
    iy = shape[0] // 2
    compare = photmjsr
    # Compare the values at the center pixel.
    ratio = output[iy, ix] / input[iy, ix]
    assert np.allclose(ratio, compare, rtol=1.e-7)


def test_apply_photom_1():
    """Test apply_photom"""

    # apply_photom() calls calc_niriss, etc., depending on EXP_TYPE.  We've
    # already tested each of these above.  The unique test in this function
    # is checking that the pixel area keywords are populated correctly.

    input_model = create_input('NIRCAM', 'NRCA3', 'NRC_IMAGE',
                               filter='F150W', pupil='CLEAR')
    ds = photom.DataSet(input_model)
    ftab = create_photom_nircam_image(min_r=8.0, max_r=9.0)

    area_ster = 2.31307642258977E-14
    area_a2 = 0.000984102303070964
    ftab.meta.photometry.pixelarea_steradians = area_ster
    ftab.meta.photometry.pixelarea_arcsecsq = area_a2

    shape = input_model.data.shape
    area_ref = create_pixel_area_ref(shape, area_ster, area_a2)

    # `apply_photom` expects its arguments to both be strings.  We're passing
    # open data models instead, to avoid the need for actual files on disk.
    # We can get away with this only because datamodels.NircamPhotomModel
    # (and other photom models) can take either an open model or the name of
    # a file as input.
    output_model = ds.apply_photom(ftab, area_ref)
    assert (math.isclose(output_model.meta.photometry.pixelarea_steradians,
                         area_ster, rel_tol=1.e-7))
    assert (math.isclose(output_model.meta.photometry.pixelarea_arcsecsq,
                         area_a2, rel_tol=1.e-7))

    # check the x,y size of area array in output is the same size as the sci data
    sh_area = output_model.area.shape
    assert shape[-1] == sh_area[-1]
    assert shape[-2] == sh_area[-2]


@pytest.mark.parametrize('srctype', ['POINT', 'EXTENDED'])
def test_apply_photom_2(srctype):
    """Test apply_photom"""

    # Check that for NIRSpec data and for an extended source, the conversion
    # factor is divided by the pixel area.

    # `apply_photom` expects its arguments to both be strings.  We're passing
    # open data models instead, to avoid the need for actual files on disk.
    # We can get away with this only because datamodels.NircamPhotomModel
    # (and other photom models) can take either an open model or the name of
    # a file as input.

    input_model = create_input('NIRSPEC', 'NRS1', 'NRS_MSASPEC',
                               filter='F170LP', grating='G235M')
    # For MSA data, the source type can vary from slit to slit, so srctype
    # is an attribute of the slit.  For other types of data, srctype is a
    # primary header keyword.
    q = 2
    xc = 7
    yc = 13
    for (k, slit) in enumerate(input_model.slits):
        # These values would normally be different for each slit,
        # but this is just a test.
        slit.source_type = srctype
        slit.quadrant = q
        slit.xcen = xc
        slit.ycen = yc

    save_input = input_model.copy()

    ds = photom.DataSet(input_model)

    ftab = create_photom_nrs_msa(min_wl=1.0, max_wl=5.0,
                                 min_r=8.0, max_r=9.0)

    # A pixel area of 2 sr is not realistic, but it's different from 1,
    # i.e. so it will make a difference when it's applied.

    quadrant = np.array([1, 2, q, 3, 4], dtype=np.int16)
    shutter_x = np.array([4, 2, xc, 11, 5], dtype=np.int16)
    shutter_y = np.array([9, 1, yc, 17, 3], dtype=np.int16)
    # I want these to be the values in steradians, but the column in the
    # pixel area reference file has values in arcsec^2.
    pixarea = np.array([991., 992., 2., 994., 995.], dtype=np.float32)
    pixarea /= A2_TO_SR                 # convert from sr to arcsec^2

    area_ref = create_msa_pixel_area_ref(quadrant, shutter_x, shutter_y,
                                         pixarea)

    output_model = ds.apply_photom(ftab, area_ref)

    rownum = find_row_in_ftab(save_input, ftab, ['filter', 'grating'],
                              slitname=None, order=None)
    photmj = ftab.phot_table['photmj'][rownum]
    nelem = ftab.phot_table['nelem'][rownum]
    wavelength = ftab.phot_table['wavelength'][rownum][0:nelem]
    relresponse = ftab.phot_table['relresponse'][rownum][0:nelem]

    result = []
    for (k, slit) in enumerate(save_input.slits):
        input = slit.data                       # this is from save_input
        output = output_model.slits[k].data

        shape = input.shape
        ix = shape[1] // 2
        iy = shape[0] // 2
        wl = slit.wavelength[iy, ix]
        rel_resp = np.interp(wl, wavelength, relresponse,
                             left=np.nan, right=np.nan)
        if slit.source_type != 'POINT':
            rel_resp /= output_model.slits[k].meta.photometry.pixelarea_steradians
        compare = photmj * rel_resp

        ratio = output[iy, ix] / input[iy, ix]
        result.append(np.allclose(ratio, compare, rtol=1.e-7))

    assert np.all(result)


def test_find_row():
    ftab = create_photom_nircam_wfss(min_wl=2.4, max_wl=5.0,
                                     min_r=8.0, max_r=9.0)
    ind = photom.find_row(ftab.phot_table, {'filter': 'F444W', 'pupil': 'GRISMR', 'order': 1})
    assert ind == 4

    # Use lower case
    ftab.phot_table[-1][1] = 'grismr'
    ind = photom.find_row(ftab.phot_table, {'filter': 'F444W', 'pupil': 'GRISMR', 'order': 1})
    assert ind == 4

    ind = photom.find_row(ftab.phot_table, {'filter': 'F444W', 'pupil': 'GRISMR', 'order': 2})
    assert ind is None
