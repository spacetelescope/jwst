import gwcs
import numpy as np
import stdatamodels.jwst.datamodels as dm
from astropy.modeling.models import Identity, Mapping, Scale, Shift

from jwst.assign_wcs.util import wcs_bbox_from_shape
from jwst.datamodels.utils.tso_multispec import make_tso_specmodel
from jwst.exp_to_source import multislit_to_container

__all__ = [
    "simple_wcs_func",
    "simple_wcs_transpose_func",
    "simple_wcs_ifu_func",
    "mock_nirspec_fs_one_slit_func",
    "mock_nirspec_mos_func",
    "mock_nirspec_bots_func",
    "mock_miri_lrs_fs_func",
    "mock_miri_ifu_func",
    "mock_miri_wfss_l2",
    "mock_nis_wfss_l2",
    "mock_nis_wfss_l3",
    "mock_niriss_soss_func",
    "mock_niriss_soss_256_func",
    "mock_niriss_soss_96_func",
    "make_spec_model",
    "make_tso_specmodel",
]


def simple_wcs_func():
    """
    Mock a simple horizontal dispersion WCS.

    Returns
    -------
    wcs : `~gwcs.wcs.WCS`
        The WCS object.
    """
    shape = (50, 50)
    xcenter = shape[1] // 2.0

    cdelt1 = 0.1
    cdelt2 = 1e-7
    cdelt3 = 0.01

    crval1 = 45.0
    crval2 = 45.1
    crval3 = 7.5

    ra = Shift(1 - xcenter) | Scale(cdelt1) | Shift(crval1)
    dec = Scale(cdelt2) | Shift(crval2)
    wave = Scale(cdelt3) | Shift(crval3)
    det2slit = Mapping((1, 1, 0), n_inputs=2) | Identity(2) & wave
    slit2world = ra & dec & Identity(1)

    input_frame = gwcs.Frame2D(name="detector")
    slit_frame = gwcs.Frame2D(name="slit_frame")
    output_frame = gwcs.Frame2D(name="world")
    pipeline = [(input_frame, det2slit), (slit_frame, slit2world), (output_frame, None)]
    wcs = gwcs.WCS(pipeline)

    # Add a bounding box
    wcs.bounding_box = wcs_bbox_from_shape(shape)

    return wcs


def simple_wcs_transpose_func():
    """
    Mock a simple vertical dispersion WCS.

    Returns
    -------
    wcs : `~gwcs.wcs.WCS`
        The WCS object.
    """
    shape = (50, 50)
    xcenter = shape[1] // 2.0

    cdelt1 = 0.1
    cdelt2 = 1e-7
    cdelt3 = -0.01

    crval1 = 45.0
    crval2 = 45.1
    crval3 = 7.5

    ra = Shift(1 - xcenter) | Scale(cdelt1) | Shift(crval1)
    dec = Scale(cdelt2) | Shift(crval2)
    wave = Scale(cdelt3) | Shift(crval3)
    mapping = Mapping((0, 0, 1), n_inputs=2)
    transform = mapping | ra & dec & wave

    input_frame = gwcs.Frame2D(name="detector")
    output_frame = gwcs.Frame2D(name="world")
    wcs = gwcs.WCS(transform, input_frame=input_frame, output_frame=output_frame)

    # Add a bounding box
    wcs.bounding_box = wcs_bbox_from_shape(shape)

    return wcs


def simple_wcs_ifu_func():
    """
    Mock a simple IFU WCS.

    Returns
    -------
    wcs : `~gwcs.wcs.WCS`
        The WCS object.
    """
    shape = (10, 50, 50)
    xcenter = shape[1] // 2.0

    cdelt1 = 0.1
    cdelt2 = 1e-7
    cdelt3 = -0.01

    crval1 = 45.0
    crval2 = 45.1
    crval3 = 7.5

    ra = Shift(1 - xcenter) | Scale(cdelt1) | Shift(crval1)
    dec = Scale(cdelt2) | Shift(crval2)
    wave = Scale(cdelt3) | Shift(crval3)
    mapping = Mapping((0, 1, 2), n_inputs=3)
    transform = mapping | ra & dec & wave

    input_frame = gwcs.CoordinateFrame(
        name="detector", naxes=3, axes_order=(0, 1, 2), axes_type=["SPATIAL", "SPATIAL", "SPECTRAL"]
    )
    output_frame = gwcs.CoordinateFrame(
        name="world", naxes=3, axes_order=(0, 1, 2), axes_type=["SPATIAL", "SPATIAL", "SPECTRAL"]
    )
    wcs = gwcs.WCS(transform, input_frame=input_frame, output_frame=output_frame)

    return wcs


def mock_nirspec_fs_one_slit_func():
    """
    Mock one slit in NIRSpec FS mode.

    Returns
    -------
    SlitModel
        The mock model.
    """
    model = dm.SlitModel()
    model.meta.instrument.name = "NIRSPEC"
    model.meta.instrument.detector = "NRS1"
    model.meta.instrument.filter = "F290LP"
    model.meta.instrument.grating = "G395H"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.instrument.fixed_slit = "S200A1"
    model.meta.exposure.nints = 1
    model.meta.exposure.type = "NRS_FIXEDSLIT"
    model.meta.subarray.name = "ALLSLITS"
    model.source_type = "EXTENDED"

    model.meta.wcsinfo.dispersion_direction = 1
    model.meta.wcs = simple_wcs_func()

    model.data = np.arange(50 * 50, dtype=float).reshape((50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05
    return model


def mock_nirspec_mos_func():
    """
    Mock three slits in NIRSpec MOS mode.

    Returns
    -------
    MultiSlitModel
        The mock model.
    """
    model = dm.MultiSlitModel()
    model.meta.instrument.name = "NIRSPEC"
    model.meta.instrument.detector = "NRS1"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.exposure.type = "NRS_MSASPEC"
    model.meta.exposure.nints = 1

    nslit = 3
    for i in range(nslit):
        slit = mock_nirspec_fs_one_slit_func()
        slit.name = str(i + 1)
        model.slits.append(slit)

    return model


def mock_nirspec_bots_func():
    """
    Mock a single slit with 10 integrations in NIRSpec BOTS mode.

    Returns
    -------
    CubeModel
        The mock model.
    """
    model = dm.CubeModel()
    model.meta.instrument.name = "NIRSPEC"
    model.meta.instrument.detector = "NRS1"
    model.meta.instrument.filter = "F290LP"
    model.meta.instrument.grating = "G395H"
    model.meta.observation.date = "2022-05-30"
    model.meta.observation.time = "01:03:16.369"
    model.meta.instrument.fixed_slit = "S1600A1"
    model.meta.exposure.type = "NRS_BRIGHTOBJ"
    model.meta.subarray.name = "SUB2048"
    model.meta.exposure.nints = 10
    model.meta.exposure.segment_number = 2
    model.meta.visit.tsovisit = True

    model.name = "S1600A1"
    model.meta.target.source_type = "POINT"
    model.meta.wcsinfo.dispersion_direction = 1
    model.meta.wcs = simple_wcs_func()
    model.meta.aperture.position_angle = 150.0

    model.data = np.arange(10 * 50 * 50, dtype=float).reshape((10, 50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05

    # Add an int_times table
    integrations = [
        (
            1,
            59729.04367729,
            59729.04378181,
            59729.04388632,
            59729.04731706,
            59729.04742158,
            59729.04752609,
        ),
        (
            2,
            59729.04389677,
            59729.04400128,
            59729.04410579,
            59729.04753654,
            59729.04764105,
            59729.04774557,
        ),
        (
            3,
            59729.04411625,
            59729.04422076,
            59729.04432527,
            59729.04775602,
            59729.04786053,
            59729.04796504,
        ),
        (
            4,
            59729.04433572,
            59729.04444023,
            59729.04454475,
            59729.04797549,
            59729.04808001,
            59729.04818452,
        ),
        (
            5,
            59729.0445552,
            59729.04465971,
            59729.04476422,
            59729.04819497,
            59729.04829948,
            59729.048404,
        ),
        (
            6,
            59729.04477467,
            59729.04487918,
            59729.0449837,
            59729.04841445,
            59729.04851896,
            59729.04862347,
        ),
        (
            7,
            59729.04499415,
            59729.04509866,
            59729.04520317,
            59729.04863392,
            59729.04873844,
            59729.04884295,
        ),
        (
            8,
            59729.04521362,
            59729.04531813,
            59729.04542265,
            59729.0488534,
            59729.04895791,
            59729.04906242,
        ),
        (
            9,
            59729.0454331,
            59729.04553761,
            59729.04564212,
            59729.04907288,
            59729.04917739,
            59729.0492819,
        ),
        (
            10,
            59729.04565257,
            59729.04575709,
            59729.0458616,
            59729.04929235,
            59729.04939686,
            59729.04950138,
        ),
    ]

    integration_table = np.array(
        integrations,
        dtype=[
            ("integration_number", "i4"),
            ("int_start_MJD_UTC", "f8"),
            ("int_mid_MJD_UTC", "f8"),
            ("int_end_MJD_UTC", "f8"),
            ("int_start_BJD_TDB", "f8"),
            ("int_mid_BJD_TDB", "f8"),
            ("int_end_BJD_TDB", "f8"),
        ],
    )
    model.int_times = integration_table

    return model


def mock_miri_lrs_fs_func():
    """
    Mock a spectral image in MIRI LRS FS mode.

    Returns
    -------
    ImageModel
        The mock model.
    """
    model = dm.ImageModel()
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "MIRIMAGE"
    model.meta.instrument.filter = "P750L"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.exposure.nints = 1
    model.meta.exposure.type = "MIR_LRS-FIXEDSLIT"
    model.meta.subarray.name = "FULL"
    model.meta.target.source_type = "EXTENDED"
    model.meta.dither.dithered_ra = 45.0
    model.meta.dither.dithered_dec = 45.0

    model.meta.wcsinfo.dispersion_direction = 2
    model.meta.wcs = simple_wcs_transpose_func()

    model.data = np.arange(50 * 50, dtype=float).reshape((50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05
    return model


def mock_miri_wfss_l2():
    """
    Mock 3 slits in MIRI WFSS mode, level 2 style.

    The slits correspond to a single exposure, with one slit per extracted source.

    Returns
    -------
    MultiSlitModel
        The mock model.
    """
    model = dm.MultiSlitModel()
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "MIRIMAGE"
    model.meta.instrument.filter = "P750L"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.observation.program_number = "1"
    model.meta.observation.observation_number = "1"
    model.meta.observation.visit_number = "1"
    model.meta.observation.visit_group = "1"
    model.meta.observation.sequence_id = "1"
    model.meta.observation.activity_id = "1"
    model.meta.observation.exposure_number = "5"
    model.meta.exposure.type = "MIR_WFSS"
    model.meta.wcsinfo = {}
    model.meta.wcsinfo.s_region = (
        "POLYGON ICRS  247.883569817 30.206692493 247.901987783 "
        "30.174116268 247.864126916 30.158804440 247.846405241 30.190721550"
    )

    slit0 = mock_miri_lrs_fs_func()

    nslit = 3
    for i in range(nslit):
        slit = slit0.copy()
        slit.name = str(i + 1)
        slit.meta.exposure.type = "MIR_WFSS"
        model.slits.append(slit)

    return model


def mock_miri_ifu_func():
    """
    Mock an IFU cube in MIRI MRS mode.

    Returns
    -------
    IFUCubeModel
        The mock model.
    """
    model = dm.IFUCubeModel()
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "MIRIFULONG"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.exposure.type = "MIR_MRS"

    model.meta.wcsinfo.dispersion_direction = 2
    model.meta.photometry.pixelarea_steradians = 1.0
    model.meta.wcs = simple_wcs_ifu_func()
    model.meta.aperture.position_angle = 150.0

    model.data = np.arange(10 * 50 * 50, dtype=float).reshape((10, 50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05
    model.weightmap = np.full_like(model.data, 1.0)
    return model


def mock_nis_wfss_l2():
    """
    Mock 3 slits in NIRISS WFSS mode, level 2 style.

    The slits correspond to a single exposure, with one slit per extracted source.

    Returns
    -------
    MultiSlitModel
        The mock model.
    """
    model = dm.MultiSlitModel()
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.instrument.filter = "GR150R"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.observation.program_number = "1"
    model.meta.observation.observation_number = "1"
    model.meta.observation.visit_number = "1"
    model.meta.observation.visit_group = "1"
    model.meta.observation.sequence_id = "1"
    model.meta.observation.activity_id = "1"
    model.meta.observation.exposure_number = "5"
    model.meta.exposure.type = "NIS_WFSS"
    model.meta.wcsinfo = {}
    model.meta.wcsinfo.s_region = (
        "POLYGON ICRS  247.883569817 30.206692493 247.901987783 "
        "30.174116268 247.864126916 30.158804440 247.846405241 30.190721550"
    )

    slit0 = mock_nirspec_fs_one_slit_func()

    nslit = 3
    for i in range(nslit):
        slit = slit0.copy()
        slit.name = str(i + 1)
        slit.meta.exposure.type = "NIS_WFSS"
        model.slits.append(slit)

    return model


def mock_nis_wfss_l3():
    """
    Mock 3 slits in NIRISS WFSS mode, level 3 style.

    Here the container has one MultiSlitModel per source, and each model has one
    slit per exposure.

    Returns
    -------
    list of SourceModelContainer
        The mock models.
    """
    model = mock_nis_wfss_l2()
    for i, slit in enumerate(model.slits):
        slit.meta.filename = f"test{i}_s2d.fits"
    container_dict = multislit_to_container([model])
    sources = list(container_dict.values())
    model.close()
    return sources


def mock_niriss_soss_func():
    """
    Mock a multi-integration cube with metadata for NIRISS SOSS mode.

    Returns
    -------
    CubeModel
        The mock model.
    """
    model = dm.CubeModel()
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.instrument.filter = "CLEAR"
    model.meta.instrument.pupil_position = 245.79
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.exposure.type = "NIS_SOSS"
    model.meta.exposure.nints = 3

    model.meta.target.source_type = "POINT"
    model.meta.wcsinfo.dispersion_direction = 1
    model.meta.wcs = simple_wcs_func()

    return model


def mock_niriss_soss_256_func():
    """
    Mock 3 integrations in NIRISS SOSS mode, subarray SUBSTRIP256.

    Returns
    -------
    CubeModel
        The mock model.
    """
    model = mock_niriss_soss_func()
    model.meta.subarray.name = "SUBSTRIP256"

    shape = (3, 256, 2048)
    model.data = np.ones(shape, dtype=np.float32)
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.err = model.data * 0.02
    model.var_poisson = model.data * 0.001
    model.var_rnoise = model.data * 0.001
    model.var_flat = model.data * 0.001
    return model


def mock_niriss_soss_96_func():
    """
    Mock 3 integrations in NIRISS SOSS mode, subarray SUBSTRIP96.

    Returns
    -------
    CubeModel
        The mock model.
    """
    model = mock_niriss_soss_func()
    model.meta.subarray.name = "SUBSTRIP96"

    shape = (3, 96, 2048)
    model.data = np.ones(shape, dtype=np.float32)

    # add a little noise to the data so fitting is more robust
    rng = np.random.default_rng(seed=77)
    noise = rng.standard_normal(shape) * 1e-6
    model.data += noise

    model.dq = np.zeros(shape, dtype=np.uint32)
    model.err = model.data * 0.02
    model.var_poisson = model.data * 0.001
    model.var_rnoise = model.data * 0.001
    model.var_flat = model.data * 0.001

    return model


def mock_niriss_soss_full_func():
    """
    Mock 3 integrations in NIRISS SOSS mode, subarray FULL.

    Returns
    -------
    CubeModel
        The mock model.
    """
    model = mock_niriss_soss_func()
    model.meta.subarray.name = "FULL"

    shape = (3, 2048, 2048)
    model.data = np.ones(shape, dtype=np.float32)
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.err = model.data * 0.02
    model.var_poisson = model.data * 0.001
    model.var_rnoise = model.data * 0.001
    model.var_flat = model.data * 0.001

    return model


def mock_niriss_soss_f277w_func():
    """
    Mock 3 integrations in NIRISS SOSS mode, filter F277W.

    Returns
    -------
    CubeModel
        The mock model.
    """
    model = dm.CubeModel((3, 3, 3))
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.instrument.filter = "F277W"
    model.meta.exposure.type = "NIS_SOSS"
    model.meta.subarray.name = "FULL"
    model.data = np.arange(27).reshape((3, 3, 3))

    return model


def make_spec_model(name="slit1", value=1.0):
    """
    Make a simple spectrum.

    Returns
    -------
    SpecModel
        The mock model.
    """
    wavelength = np.arange(20, dtype=np.float32)
    flux = np.full(10, value)
    error = 0.05 * flux
    f_var_poisson = error**2
    f_var_rnoise = np.zeros_like(flux)
    f_var_flat = np.zeros_like(flux)
    surf_bright = flux / 10
    sb_error = error / 10
    sb_var_poisson = f_var_poisson / 10**2
    sb_var_rnoise = f_var_rnoise / 10**2
    sb_var_flat = f_var_flat / 10**2
    dq = np.zeros(20, dtype=np.uint32)
    background = np.zeros_like(flux)
    berror = np.zeros_like(flux)
    b_var_poisson = np.zeros_like(flux)
    b_var_rnoise = np.zeros_like(flux)
    b_var_flat = np.zeros_like(flux)
    npixels = np.full(20, 10)

    spec_dtype = dm.SpecModel().spec_table.dtype
    otab = np.array(
        list(
            zip(
                wavelength,
                flux,
                error,
                f_var_poisson,
                f_var_rnoise,
                f_var_flat,
                surf_bright,
                sb_error,
                sb_var_poisson,
                sb_var_rnoise,
                sb_var_flat,
                dq,
                background,
                berror,
                b_var_poisson,
                b_var_rnoise,
                b_var_flat,
                npixels,
                strict=False,
            ),
        ),
        dtype=spec_dtype,
    )

    spec_model = dm.SpecModel(spec_table=otab)
    spec_model.name = name

    return spec_model


def make_tso_spec_model(n_spectra=10):
    """
    Make a multi-integration, multi-spec model.

    Parameters
    ----------
    n_spectra : int, optional
        The number of integrations to generate spectra for.

    Returns
    -------
    TSOMultiSpecModel
        The composite spectral model.
    """
    spec_list = []
    for i in range(n_spectra):
        spec_model = make_spec_model(name=f"slit{i + 1}", value=i + 1)
        spec_list.append(spec_model)

    tso_spec = make_tso_specmodel(spec_list)
    model = dm.TSOMultiSpecModel()
    model.spec.append(tso_spec)
    return model
