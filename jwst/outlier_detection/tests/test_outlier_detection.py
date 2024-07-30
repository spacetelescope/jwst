import pytest
import numpy as np
from scipy.ndimage import gaussian_filter
from glob import glob
import os

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.outlier_detection import OutlierDetectionStep
from jwst.outlier_detection.utils import flag_resampled_model_crs
from jwst.outlier_detection.outlier_detection_step import (
    IMAGE_MODES,
    TSO_SPEC_MODES,
    TSO_IMAGE_MODES,
    CORON_IMAGE_MODES,
)
from jwst.assign_wcs.pointing import create_fitswcs

OUTLIER_DO_NOT_USE = np.bitwise_or(
    datamodels.dqflags.pixel["DO_NOT_USE"], datamodels.dqflags.pixel["OUTLIER"]
)

# TSO types to test
exptypes_tso = [(exptype, True) for exptype in TSO_SPEC_MODES + TSO_IMAGE_MODES]
exptypes_tso.append(("MIR_IMAGE", True))
# CORON types to test
exptypes_coron = [(exptype, False) for exptype in CORON_IMAGE_MODES]


@pytest.fixture
def sci_blot_image_pair():
    """Provide a science and blotted ImageModel pair."""
    shape = (20, 20)
    sigma = 0.02
    background = 3

    sci = datamodels.ImageModel(shape)

    # Populate keywords
    sci.meta.exposure.exposure_time = 1
    sci.meta.background.subtracted = False
    sci.meta.background.level = background

    rng = np.random.default_rng(720)
    sci.data = rng.normal(loc=background, size=shape, scale=sigma)
    sci.err = np.zeros(shape) + sigma
    sci.var_rnoise += 0

    # Add a source in the center
    signal = 20 * sigma
    sci.data[10, 10] += signal
    # update the noise for this source to include the photon/measurement noise
    sci.err[10, 10] = np.sqrt(sigma ** 2 + signal)

    # The blot image is just a smoothed version of the science image that has
    # its background subtracted
    blot = sci.copy()
    blot.data = gaussian_filter(blot.data, sigma=3)
    blot.data -= background

    return sci, blot


def test_flag_cr(sci_blot_image_pair):
    """Test the flag_cr function.  Test logic, not the actual noise model."""
    sci, blot = sci_blot_image_pair
    assert (sci.dq == 0).all()

    # Drop some CRs on the science array
    sci.data[3, 3] += 100
    sci.data[3, 7] += 1e3
    sci.data[7, 3] += 1e4
    sci.data[7, 7] += 1e5

    # run flag_cr() which updates in-place.  Copy sci first.
    data_copy = sci.data.copy()
    flag_resampled_model_crs(
        sci,
        blot.data,
        5.0,
        4.0,
        1.2,
        0.7,
        0,
    )

    # Make sure science data array is unchanged after flag_cr()
    np.testing.assert_allclose(sci.data, data_copy)

    # Verify that both DQ flags are set in the DQ array for all outliers
    assert sci.dq[3, 3] == OUTLIER_DO_NOT_USE
    assert sci.dq[3, 7] == OUTLIER_DO_NOT_USE
    assert sci.dq[7, 3] == OUTLIER_DO_NOT_USE
    assert sci.dq[7, 7] == OUTLIER_DO_NOT_USE

    # Verify the source wasn't flagged
    assert sci.dq[10, 10] == datamodels.dqflags.pixel["GOOD"]


# not a fixture - now has options
def we_many_sci(
    numsci=3, sigma=0.02, background=1.5, signal=7, exptype="MIR_IMAGE", tsovisit=False
):
    """Provide numsci science images with different noise but identical source
    and same background level"""
    shape = (20, 20)

    sci1 = datamodels.ImageModel(shape)

    # Populate keywords
    sci1.meta.instrument.name = "MIRI"
    sci1.meta.instrument.detector = "MIRIMAGE"
    sci1.meta.exposure.type = exptype
    sci1.meta.visit.tsovisit = tsovisit
    sci1.meta.observation.date = "2020-01-01"
    sci1.meta.observation.time = "00:00:00"
    sci1.meta.telescope = "JWST"
    sci1.meta.exposure.exposure_time = 1
    sci1.meta.wcsinfo.wcsaxes = 2
    sci1.meta.wcsinfo.ctype1 = "RA---TAN"
    sci1.meta.wcsinfo.ctype2 = "DEC--TAN"
    sci1.meta.wcsinfo.cdelt1 = 3e-6
    sci1.meta.wcsinfo.cdelt2 = 3e-6
    sci1.meta.wcsinfo.roll_ref = 0
    sci1.meta.wcsinfo.ra_ref = 1.5e-5
    sci1.meta.wcsinfo.dec_ref = 1.5e-5
    sci1.meta.wcsinfo.v3yangle = 0
    sci1.meta.wcsinfo.vparity = -1
    sci1.meta.wcsinfo.pc1_1 = 1
    sci1.meta.wcsinfo.pc1_2 = 0
    sci1.meta.wcsinfo.pc2_1 = 0
    sci1.meta.wcsinfo.pc2_2 = 1
    sci1.meta.wcsinfo.crpix1 = 5
    sci1.meta.wcsinfo.crpix2 = 5
    sci1.meta.wcsinfo.crval1 = 0
    sci1.meta.wcsinfo.crval2 = 0
    sci1.meta.wcsinfo.cunit1 = "deg"
    sci1.meta.wcsinfo.cunit2 = "deg"
    sci1.meta.background.subtracted = False
    sci1.meta.background.level = background

    # Replace the FITS-type WCS with an Identity WCS
    sci1.meta.wcs = create_fitswcs(sci1)
    rng = np.random.default_rng(720)
    sci1.data = rng.normal(loc=background, size=shape, scale=sigma)
    sci1.err = np.zeros(shape) + sigma
    sci1.data[7, 7] += signal
    # update the noise for this source to include the photon/measurement noise
    sci1.err[7, 7] = np.sqrt(sigma ** 2 + signal)
    sci1.var_rnoise = np.zeros(shape) + 1.0
    sci1.meta.filename = "foo1_cal.fits"

    # add pixel areas
    sci1.meta.photometry.pixelarea_steradians = 1.0
    sci1.meta.photometry.pixelarea_arcsecsq = 1.0

    # Make copies with different noise
    all_sci = [sci1]
    for i in range(numsci - 1):
        tsci = sci1.copy()
        tsci.data = rng.normal(loc=background, size=shape, scale=sigma)
        # Add a source in the center
        tsci.data[7, 7] += signal
        tsci.meta.filename = f"foo{i + 2}_cal.fits"
        all_sci.append(tsci)

    return all_sci


def container_to_cube(container):
    """
    Utility function to convert the test container to a cube
    for some tests
    """
    cube_data = np.array([i.data for i in container])
    cube_err = np.array([i.err for i in container])
    cube_dq = np.array([i.dq for i in container])
    cube_var_noise = np.array([i.var_rnoise for i in container])
    cube = datamodels.CubeModel(data=cube_data, err=cube_err, dq=cube_dq, var_noise=cube_var_noise)

    # update metadata of cube to match the first image
    cube.meta = container[0].meta
    return cube


@pytest.fixture
def we_three_sci():
    """Provide 3 science images with different noise but identical source
    and same background level"""
    return we_many_sci(numsci=3)


def test_outlier_step_no_outliers(we_three_sci, tmp_cwd):
    """Test whole step, no outliers"""
    container = ModelContainer(list(we_three_sci))
    pristine = ModelContainer([m.copy() for m in container])
    OutlierDetectionStep.call(container, in_memory=True)

    # Make sure nothing changed in SCI and DQ arrays
    for image, uncorrected in zip(pristine, container):
        np.testing.assert_allclose(image.data, uncorrected.data)
        np.testing.assert_allclose(image.dq, uncorrected.dq)


def test_outlier_step_base(we_three_sci, tmp_cwd):
    """Test whole step with an outlier including saving intermediate and results files"""
    container = ModelLibrary(list(we_three_sci))

    # Drop a CR on the science array
    with container:
        zeroth = container.borrow(0)
        zeroth.data[12, 12] += 1
        container.shelve(zeroth)

    # Verify that intermediary files are removed
    OutlierDetectionStep.call(container, in_memory=True)
    i2d_files = glob(os.path.join(tmp_cwd, '*i2d.fits'))
    median_files = glob(os.path.join(tmp_cwd, '*median.fits'))
    assert len(i2d_files) == 0
    assert len(median_files) == 0

    # Save all the data into a separate array before passing into step
    data_as_cube = []
    with container:
        for model in container:
            data_as_cube.append(model.data)
            container.shelve(model, modify=False)

    result = OutlierDetectionStep.call(
        container, save_results=True, save_intermediate_results=True, in_memory=True
    )

    # Make sure nothing changed in SCI array
    with result:
        for i, corrected in enumerate(result):
            np.testing.assert_allclose(data_as_cube[i], corrected.data)
            result.shelve(corrected, modify=False)

    # Verify source is not flagged
    with result:
        for r in result:
            assert r.dq[7, 7] == datamodels.dqflags.pixel["GOOD"]
            result.shelve(r, modify=False)

    # Verify CR is flagged
    with result:
        zeroth = result.borrow(0)
        assert zeroth.dq[12, 12] == OUTLIER_DO_NOT_USE
        result.shelve(zeroth, modify=False)

    # Verify that intermediary files are saved at the specified location
    i2d_files = glob(os.path.join(tmp_cwd, '*i2d.fits'))
    median_files = glob(os.path.join(tmp_cwd, '*median.fits'))
    assert len(i2d_files) != 0
    assert len(median_files) != 0


def test_outlier_step_on_disk(we_three_sci, tmp_cwd):
    """Test whole step with an outlier including saving intermediate and results files"""

    for model in we_three_sci:
        model.save(model.meta.filename)
    filenames = [model.meta.filename for model in we_three_sci]
    # Drop a CR on the science array
    with datamodels.open(filenames[0]) as dm0:
        dm0.data[12, 12] += 1
        dm0.write(dm0.meta.filename)

    # Initialize inputs for the test based on filenames only
    # needs to be an asn for ModelLibrary to load it in on_disk mode
    asn = {
    'asn_type': 'test',
    'asn_id': 'o001',
    'products': [
        {
            'name': 'product_a',
            'members': [
                {'expname': filenames[0], 'exptype': 'science'},
                {'expname': filenames[1], 'exptype': 'science'},
                {'expname': filenames[2], 'exptype': 'science'},
            ]
        },
    ]
}
    container = ModelLibrary(asn, on_disk=True)

    # Save all the data into a separate array before passing into step
    data_as_cube = []
    with container:
        for model in container:
            data_as_cube.append(model.data)
            container.shelve(model, modify=False)

    result = OutlierDetectionStep.call(
        container, save_results=True, save_intermediate_results=True, in_memory=False
    )

    # Make sure nothing changed in SCI array
    with result:
        for i, corrected in enumerate(result):
            np.testing.assert_allclose(data_as_cube[i], corrected.data)
            result.shelve(corrected, modify=False)

    # Verify source is not flagged
    with result:
        for r in result:
            assert r.dq[7, 7] == datamodels.dqflags.pixel["GOOD"]
            result.shelve(r, modify=False)

    # Verify CR is flagged
    with result:
        zeroth = result.borrow(0)
        assert zeroth.dq[12, 12] == OUTLIER_DO_NOT_USE
        result.shelve(zeroth, modify=False)


def test_outlier_step_square_source_no_outliers(we_three_sci, tmp_cwd):
    """Test whole step with square source with sharp edges, no outliers"""
    container = ModelLibrary(list(we_three_sci))

    # put a square source in all three exposures
    with container:
        for ccont in container:
            ccont.data[5:15, 5:15] += 1e3
            container.shelve(ccont)

    # Save all the data into a separate array before passing into step
    data_as_cube = []
    dq_as_cube = []
    with container:
        for model in container:
            data_as_cube.append(model.data)
            dq_as_cube.append(model.dq)
            container.shelve(model, modify=False)

    result = OutlierDetectionStep.call(container, in_memory=True)

    # Make sure nothing changed in SCI and DQ arrays
    with container:
        for i, image in enumerate(container):
            np.testing.assert_allclose(image.data, data_as_cube[i])
            np.testing.assert_allclose(image.dq, dq_as_cube[i])
            container.shelve(image, modify=False)

    # Make sure nothing changed in SCI and DQ arrays
    with result:
        for i, corrected in enumerate(result):
            np.testing.assert_allclose(data_as_cube[i], corrected.data)
            np.testing.assert_allclose(dq_as_cube[i], corrected.dq)
            result.shelve(corrected, modify=False)


@pytest.mark.parametrize("exptype", IMAGE_MODES)
def test_outlier_step_image_weak_CR_dither(exptype, tmp_cwd):
    """Test whole step with an outlier for imaging modes"""
    bkg = 1.5
    sig = 0.02
    container = ModelLibrary(
        we_many_sci(background=bkg, sigma=sig, signal=7.0, exptype=exptype)
    )

    # Drop a weak CR on the science array
    # no noise so it should always be above the default threshold of 5
    with container:
        zeroth = container.borrow(0)
        zeroth.data[12, 12] = bkg + sig * 10
        container.shelve(zeroth)

    # Save all the data into a separate array before passing into step
    data_as_cube = []
    with container:
        for model in container:
            data_as_cube.append(model.data)
            container.shelve(model, modify=False)

    result = OutlierDetectionStep.call(container, in_memory=True)

    # Make sure nothing changed in SCI array
    with result:
        for i, corrected in enumerate(result):
            np.testing.assert_allclose(data_as_cube[i], corrected.data)
            result.shelve(corrected, modify=False)

    # Verify source is not flagged
    with result:
        for r in result:
            assert r.dq[7, 7] == datamodels.dqflags.pixel["GOOD"]
            result.shelve(r, modify=False)

    # Verify CR is flagged
    with result:
        example = result.borrow(0)
        assert example.dq[12, 12] == OUTLIER_DO_NOT_USE
        result.shelve(example, modify=False)


@pytest.mark.parametrize("exptype, tsovisit", exptypes_coron)
def test_outlier_step_image_weak_CR_coron(exptype, tsovisit, tmp_cwd):
    """Test whole step with an outlier for coronagraphic modes"""
    bkg = 1.5
    sig = 0.02
    container = ModelContainer(
        we_many_sci(
            background=bkg, sigma=sig, signal=7.0, exptype=exptype, tsovisit=tsovisit
        )
    )

    # Drop a weak CR on the science array
    # no noise so it should always be above the default threshold of 5
    container[0].data[12, 12] = bkg + sig * 10

    # coron3 will provide a CubeModel so convert the container to a cube
    cube = container_to_cube(container)

    result = OutlierDetectionStep.call(cube)

    # Make sure nothing changed in SCI array
    for i, image in enumerate(container):
        np.testing.assert_allclose(image.data, result.data[i])

    # Verify source is not flagged
    assert np.all(result.dq[:, 7, 7] == datamodels.dqflags.pixel["GOOD"])

    # Verify CR is flagged
    assert np.all(result.dq[0, 12, 12] == OUTLIER_DO_NOT_USE)


@pytest.mark.parametrize("exptype, tsovisit", exptypes_tso)
def test_outlier_step_weak_cr_tso(exptype, tsovisit):
    '''Test outlier detection with rolling median on time-varying source
    This test fails if rolling_window_width is set to 100, i.e., take simple median
    '''
    bkg = 1.5
    sig = 0.02
    rolling_window_width = 7
    numsci = 50
    signal = 7.0
    im = we_many_sci(
        numsci=numsci, background=bkg, sigma=sig, signal=signal, exptype=exptype, tsovisit=tsovisit
    )

    # Drop a weak CR on the science array
    cr_timestep = 5
    im[cr_timestep].data[12, 12] = bkg + sig * 10

    # make time variability that has larger total amplitude than
    # the CR signal but deviations frame-by-frame are smaller
    real_time_variability = signal * np.cos(np.linspace(0, np.pi, numsci))
    for i, model in enumerate(im):
        model.data[7, 7] += real_time_variability[i]
        model.err[7, 7] = np.sqrt(sig ** 2 + model.data[7, 7])

    cube = container_to_cube(im)

    result = OutlierDetectionStep.call(cube, rolling_window_width=rolling_window_width)

    # Make sure nothing changed in SCI array
    for i, model in enumerate(im):
        np.testing.assert_allclose(model.data, result.data[i])

    # Verify source is not flagged
    assert np.all(result.dq[:, 7, 7] == datamodels.dqflags.pixel["GOOD"])

    # Verify CR is flagged
    assert result.dq[cr_timestep, 12, 12] == OUTLIER_DO_NOT_USE
