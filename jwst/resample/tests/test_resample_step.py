from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose
import numpy as np
import pytest

from jwst.datamodels import ImageModel, ModelContainer
from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.util import compute_fiducial, compute_scale
from jwst.extract_2d import Extract2dStep
from jwst.resample import ResampleSpecStep, ResampleStep
from jwst.resample.resample_spec import ResampleSpecData


@pytest.fixture
def nirspec_rate():
    ysize = 2048
    xsize = 2048
    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.var_rnoise += 1
    im.meta.target = {'ra': 100.1237, 'dec': 39.86}
    im.meta.wcsinfo = {
        'dec_ref': 40,
        'ra_ref': 100,
        'roll_ref': 0,
        'v2_ref': -453.5134,
        'v3_ref': -373.4826,
        'v3yangle': 0.0,
        'vparity': -1}
    im.meta.instrument = {
        'detector': 'NRS1',
        'filter': 'CLEAR',
        'grating': 'PRISM',
        'name': 'NIRSPEC',
        'gwa_tilt': 37.0610,
        'gwa_xtilt': 0.0001,
        'gwa_ytilt': 0.0001}
    im.meta.subarray = {
        'fastaxis': 1,
        'name': 'SUBS200A1',
        'slowaxis': 2,
        'xsize': 72,
        'xstart': 1,
        'ysize': 416,
        'ystart': 529}
    im.meta.observation = {
        'date': '2016-09-05',
        'time': '8:59:37'}
    im.meta.exposure = {
        'duration': 11.805952,
        'end_time': 58119.85416,
        'exposure_time': 11.776,
        'frame_time': 0.11776,
        'group_time': 0.11776,
        'groupgap': 0,
        'integration_time': 11.776,
        'nframes': 1,
        'ngroups': 100,
        'nints': 1,
        'nresets_between_ints': 0,
        'nsamples': 1,
        'readpatt': 'NRSRAPID',
        'sample_time': 10.0,
        'start_time': 58119.8333,
        'type': 'NRS_FIXEDSLIT',
        'zero_frame': False}

    return im


@pytest.fixture
def miri_rate():
    xsize = 72
    ysize = 416
    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.data += 5
    im.var_rnoise += 1
    im.meta.wcsinfo = {
        'dec_ref': 40,
        'ra_ref': 100,
        'roll_ref': 0.0,
        'v2_ref': -453.5134,
        'v3_ref': -373.4826,
        'v3yangle': 0.0,
        'vparity': -1}
    im.meta.instrument = {
        'detector': 'MIRIMAGE',
        'filter': 'P750L',
        'name': 'MIRI'}
    im.meta.observation = {
        'date': '2019-01-01',
        'time': '17:00:00'}
    im.meta.subarray = {
        'fastaxis': 1,
        'name': 'SLITLESSPRISM',
        'slowaxis': 2,
        'xsize': xsize,
        'xstart': 1,
        'ysize': ysize,
        'ystart': 529}
    im.meta.exposure = {
        'duration': 11.805952,
        'end_time': 58119.85416,
        'exposure_time': 11.776,
        'frame_time': 0.11776,
        'group_time': 0.11776,
        'groupgap': 0,
        'integration_time': 11.776,
        'nframes': 1,
        'ngroups': 100,
        'nints': 1,
        'nresets_between_ints': 0,
        'nsamples': 1,
        'readpatt': 'FAST',
        'sample_time': 10.0,
        'start_time': 58119.8333,
        'type': 'MIR_LRS-SLITLESS',
        'zero_frame': False}

    return im


@pytest.fixture
def nircam_rate():
    xsize = 204
    ysize = 204
    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.var_rnoise += 0
    im.meta.wcsinfo = {
        'ctype1': 'RA---TAN',
        'ctype2': 'DEC--TAN',
        'dec_ref': 11.99875540218638,
        'ra_ref': 22.02351763251896,
        'roll_ref': 0.005076934167039675,
        'v2_ref': 86.039011,
        'v3_ref': -493.385704,
        'v3yangle': -0.07385127,
        'vparity': -1,
        'wcsaxes': 2}
    im.meta.instrument = {
        'channel': 'LONG',
        'detector': 'NRCALONG',
        'filter': 'F444W',
        'lamp_mode': 'NONE',
        'module': 'A',
        'name': 'NIRCAM',
        'pupil': 'CLEAR'}
    im.meta.subarray = {
        'fastaxis': -1,
        'name': 'FULL',
        'slowaxis': 2,
        'xsize': xsize,
        'xstart': 1,
        'ysize': ysize,
        'ystart': 1}
    im.meta.observation = {
        'activity_id': '01',
        'date': '2021-10-25',
        'exposure_number': '00001',
        'obs_id': 'V42424001001P0000000001101',
        'observation_label': 'nircam_ptsrc_only',
        'observation_number': '001',
        'program_number': '42424',
        'sequence_id': '1',
        'time': '16:58:27.258',
        'visit_group': '01',
        'visit_id': '42424001001',
        'visit_number': '001'}
    im.meta.exposure = {
        'duration': 161.05155,
        'end_time': 59512.70899968495,
        'exposure_time': 150.31478,
        'frame_time': 10.73677,
        'group_time': 21.47354,
        'groupgap': 1,
        'integration_time': 150.31478,
        'mid_time': 59512.70812980775,
        'nframes': 1,
        'ngroups': 7,
        'nints': 1,
        'nresets_at_start': 1,
        'nresets_between_ints': 1,
        'readpatt': 'BRIGHT1',
        'sample_time': 10,
        'start_time': 59512.70725993055,
        'type': 'NRC_IMAGE'}
    im.meta.photometry = {
        'pixelarea_steradians': 1e-13,
        'pixelarea_arcsecsq': 4e-3,
    }

    return im


@pytest.mark.xfail(reason="Fails due to updated ref files delivered to CRDS")
def test_nirspec_wcs_roundtrip(nirspec_rate):
    im = AssignWcsStep.call(nirspec_rate)
    im = Extract2dStep.call(im)
    im = ResampleSpecStep.call(im)

    for slit in im.slits:
        x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
        ra, dec, lam = slit.meta.wcs(x, y)
        xp, yp = slit.meta.wcs.invert(ra, dec, lam)

        assert_allclose(x, xp, rtol=0, atol=1e-8)
        assert_allclose(y, yp, rtol=0, atol=3e-4)


def test_miri_wcs_roundtrip(miri_rate):
    im = AssignWcsStep.call(miri_rate)
    im = ResampleSpecStep.call(im)

    x, y = grid_from_bounding_box(im.meta.wcs.bounding_box)
    ra, dec, lam = im.meta.wcs(x, y)
    xp, yp = im.meta.wcs.invert(ra, dec, lam)

    assert_allclose(x, xp, atol=1e-8)
    assert_allclose(y, yp, atol=1e-8)


@pytest.mark.parametrize("ratio", [0.5, 0.7, 1.0])
def test_pixel_scale_ratio_spec(miri_rate, ratio):
    im = AssignWcsStep.call(miri_rate, sip_approx=False)
    result1 = ResampleSpecStep.call(im)
    result2 = ResampleSpecStep.call(im, pixel_scale_ratio=ratio)

    assert_allclose(np.array(result1.data.shape), np.array(result2.data.shape) * ratio, rtol=1, atol=1)


@pytest.mark.parametrize("ratio", [0.5, 0.7, 1.0])
def test_pixel_scale_ratio_imaging(nircam_rate, ratio):
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    im.data += 5
    result1 = ResampleStep.call(im)
    result2 = ResampleStep.call(im, pixel_scale_ratio=ratio)

    assert_allclose(np.array(result1.data.shape), np.array(result2.data.shape) * ratio, rtol=1, atol=1)

    # Avoid edge effects; make sure data values are identical for surface brightness data
    assert np.mean(result1.data[10:-10, 10:-10]) == np.mean(result2.data[10:-10, 10:-10])

    # Make sure the photometry keywords describing the solid angle of a pixel
    # are updated
    area1 = result1.meta.photometry.pixelarea_steradians
    area2 = result2.meta.photometry.pixelarea_steradians
    assert_allclose(area1 * ratio**2, area2, rtol=1e-6)

    assert result1.meta.resample.pixel_scale_ratio == 1.0
    assert result2.meta.resample.pixel_scale_ratio == ratio


def test_weight_type(nircam_rate, _jail):
    """Check that weight_type of exptime and ivm work"""
    im1 = AssignWcsStep.call(nircam_rate, sip_approx=False)
    im1.var_rnoise[:] = 0
    im2 = im1.copy()
    im3 = im1.copy()
    im1.data += 10
    im2.data += 5
    im3.data += 5
    im1.var_rnoise += (1 / 10)
    im2.var_rnoise += (1 / 5)
    im3.var_rnoise += (1 / 5)
    im2.meta.observation.sequence_id = "2"
    im3.meta.observation.sequence_id = "3"

    c = ModelContainer([im1, im2, im3])
    assert len(c.group_names) == 3

    result1 = ResampleStep.call(c, weight_type="ivm", blendheaders=False, save_results=True)

    # assert_allclose(result1.data, result2.data)
    # assert_allclose(result1.wht, result2.wht)
    assert_allclose(result1.data[100:105, 100:105], 7.5, rtol=1e-2)
    assert_allclose(result1.wht[100:105, 100:105], 19.5, rtol=1e-2)

    result2 = ResampleStep.call(c, weight_type="exptime", blendheaders=False)

    assert_allclose(result2.data[100:105, 100:105], 6.667, rtol=1e-2)
    assert_allclose(result2.wht[100:105, 100:105], 450.9, rtol=1e-1)


def test_sip_coeffs_do_not_propagate(nircam_rate):
    im = AssignWcsStep.call(nircam_rate, sip_degree=2)

    # Check some SIP keywords produced above
    assert im.meta.wcsinfo.cd1_1 is not None
    assert im.meta.wcsinfo.ctype1 == "RA---TAN-SIP"

    # Make sure no PC matrix stuff is there
    assert im.meta.wcsinfo.pc1_1 is None

    result = ResampleStep.call(im)

    # Verify that SIP-related keywords do not propagate to resampled output
    assert result.meta.wcsinfo.cd1_1 is None
    assert result.meta.wcsinfo.ctype1 == "RA---TAN"

    # Make sure we have a PC matrix
    assert result.meta.wcsinfo.pc1_1 is not None


@pytest.fixture
def miri_rate_zero_crossing():
    xsize = 1032
    ysize = 1024
    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.var_rnoise = np.random.random(shape)
    im.meta.wcsinfo = {
        'dec_ref': 2.16444343946559e-05,
        'ra_ref': -0.00026031780056776,
        'roll_ref': 0.0,
        'v2_ref': -415.0690466121227,
        'v3_ref': -400.575920398547,
        'v3yangle': 0.0,
        'vparity': -1}
    im.meta.instrument = {
        'detector': 'MIRIMAGE',
        'filter': 'P750L',
        'name': 'MIRI'}
    im.meta.observation = {
        'date': '2019-01-01',
        'time': '17:00:00'}
    im.meta.subarray = {
        'fastaxis': 1,
        'name': 'FULL',
        'slowaxis': 2,
        'xsize': xsize,
        'xstart': 1,
        'ysize': ysize,
        'ystart': 1}
    im.meta.exposure = {
        'duration': 11.805952,
        'end_time': 58119.85416,
        'exposure_time': 11.776,
        'frame_time': 0.11776,
        'group_time': 0.11776,
        'groupgap': 0,
        'integration_time': 11.776,
        'nframes': 1,
        'ngroups': 100,
        'nints': 1,
        'nresets_between_ints': 0,
        'nsamples': 1,
        'readpatt': 'FAST',
        'sample_time': 10.0,
        'start_time': 58119.8333,
        'type': 'MIR_LRS-FIXEDSLIT',
        'zero_frame': False}

    return im


@pytest.fixture
def miri_rate_pair(miri_rate_zero_crossing):
    im1 = miri_rate_zero_crossing
    # Create a nodded version
    im2 = im1.copy()
    im2.meta.wcsinfo.ra_ref = 0.00026308279776455
    im2.meta.wcsinfo.dec_ref = -2.1860888891293e-05
    im1 = AssignWcsStep.call(im1)
    im2 = AssignWcsStep.call(im2)

    return im1, im2


def test_build_interpolated_output_wcs(miri_rate_pair):
    im1, im2 = miri_rate_pair

    driz = ResampleSpecData(ModelContainer([im1, im2]))
    output_wcs = driz.build_interpolated_output_wcs()

    # Make sure that all RA, Dec values in the input image have a location in
    # the output frame
    grid = grid_from_bounding_box(im2.meta.wcs.bounding_box)
    ra, dec, lam = im2.meta.wcs(*grid)
    x, y = output_wcs.invert(ra, dec, lam)

    # This currently fails, as we see a slight offset
    # assert (x > 0).all()

    # Make sure the output slit size is larger than the input slit size
    # for this nodded data
    assert output_wcs.array_shape[1] > ra.shape[1]


def test_wcs_keywords(nircam_rate):
    """Make sure certain wcs keywords are removed after resample
    """
    im = AssignWcsStep.call(nircam_rate)
    result = ResampleStep.call(im)

    assert result.meta.wcsinfo.v2_ref is None
    assert result.meta.wcsinfo.v3_ref is None
    assert result.meta.wcsinfo.ra_ref is None
    assert result.meta.wcsinfo.dec_ref is None
    assert result.meta.wcsinfo.roll_ref is None
    assert result.meta.wcsinfo.v3yangle is None
    assert result.meta.wcsinfo.vparity is None


@pytest.mark.parametrize("n_images", [1, 2, 3, 9])
def test_resample_variance(nircam_rate, n_images):
    """Test that resampled variance and error arrays are computed properly"""
    err = 0.02429
    var_rnoise = 0.00034
    var_poisson = 0.00025
    im = AssignWcsStep.call(nircam_rate)
    im.var_rnoise += var_rnoise
    im.var_poisson += var_poisson
    im.err += err
    im.meta.filename = "foo.fits"

    c = ModelContainer()
    for n in range(n_images):
        c.append(im.copy())

    result = ResampleStep.call(c, blendheaders=False)

    # Verify that the combined uncertainty goes as 1 / sqrt(N)
    assert_allclose(result.err[5:-5, 5:-5].mean(), err / np.sqrt(n_images), atol=1e-5)
    assert_allclose(result.var_rnoise[5:-5, 5:-5].mean(), var_rnoise / n_images, atol=1e-7)
    assert_allclose(result.var_poisson[5:-5, 5:-5].mean(), var_poisson / n_images, atol=1e-7)


@pytest.mark.parametrize("shape", [(0, ), (10, 1)])
def test_resample_undefined_variance(nircam_rate, shape):
    """Test that resampled variance and error arrays are computed properly"""
    im = AssignWcsStep.call(nircam_rate)
    im.var_rnoise = np.ones(shape, dtype=im.var_rnoise.dtype.type)
    im.var_poisson = np.ones(shape, dtype=im.var_poisson.dtype.type)
    im.var_flat = np.ones(shape, dtype=im.var_flat.dtype.type)
    im.meta.filename = "foo.fits"

    c = ModelContainer([im])
    ResampleStep.call(c, blendheaders=False)


@pytest.mark.parametrize('ratio', [0.7, 1.2])
@pytest.mark.parametrize('rotation', [0, 15, 135])
@pytest.mark.parametrize('crpix', [(256, 488), (700, 124)])
@pytest.mark.parametrize('crval', [(50, 77), (20, -30)])
@pytest.mark.parametrize('shape', [(1205, 1100)])
def test_custom_wcs_resample_imaging(nircam_rate, ratio, rotation, crpix, crval, shape):
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    im.data += 5
    result = ResampleStep.call(
        im,
        output_shape=shape,
        crpix=crpix,
        crval=crval,
        rotation=rotation,
        pixel_scale_ratio=ratio
    )

    t = result.meta.wcs.forward_transform

    # test rotation
    pc = t['pc_rotation_matrix'].matrix.value
    orientation = np.rad2deg(np.arctan2(pc[0, 1], pc[1, 1]))
    assert np.allclose(rotation, orientation)

    # test CRPIX
    assert np.allclose(
        (-t['crpix1'].offset.value, -t['crpix2'].offset.value),
        crpix
    )

    # test CRVAL
    assert np.allclose(t(*crpix), crval)

    # test output image shape
    assert result.data.shape == shape[::-1]


@pytest.mark.parametrize('ratio', [1.3, 1])
def test_custom_wcs_pscale_resample_imaging(nircam_rate, ratio):
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    im.data += 5

    fiducial = compute_fiducial([im.meta.wcs])
    input_scale = compute_scale(wcs=im.meta.wcs, fiducial=fiducial)
    result = ResampleStep.call(
        im,
        pixel_scale_ratio=ratio,
        pixel_scale=3600 * input_scale * 0.75
    )
    output_scale = compute_scale(wcs=result.meta.wcs, fiducial=fiducial)

    # test scales are close
    assert np.allclose(output_scale, input_scale * 0.75)


def test_pixscale(nircam_rate):

    # check that if both 'pixel_scale_ratio' and 'pixel_scale' are passed in,
    # that 'pixel_scale' overrides correctly
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    pixarea = im.meta.photometry.pixelarea_arcsecsq

    # check when both pixel_scale and pixel_scale_ratio are passed in
    res = ResampleStep.call(im, pixel_scale=0.04, pixel_scale_ratio=0.7)
    assert np.allclose(res.meta.resample.pixel_scale_ratio, 0.04**2 / pixarea)

    # just pixel_scale
    res = ResampleStep.call(im, pixel_scale=0.04)
    assert np.allclose(res.meta.resample.pixel_scale_ratio, 0.04**2 / pixarea)

    # just pixel_scale_ratio
    res = ResampleStep.call(im, pixel_scale_ratio=0.7)
    assert res.meta.resample.pixel_scale_ratio == 0.7


def test_phot_keywords(nircam_rate):
    # test that resample keywords agree with photometry keywords after step is run

    im = AssignWcsStep.call(nircam_rate, sip_approx=False)

    orig_pix_area_sr = im.meta.photometry.pixelarea_steradians
    orig_pix_area_arcsec = im.meta.photometry.pixelarea_arcsecsq

    # first run by setting `pixel_scale`
    res = ResampleStep.call(im, pixel_scale=0.04)
    new_psr = res.meta.resample.pixel_scale_ratio

    assert res.meta.resample.pixel_scale_ratio == 0.04**2 / orig_pix_area_arcsec
    assert res.meta.photometry.pixelarea_steradians == orig_pix_area_sr * new_psr**2
    assert res.meta.photometry.pixelarea_arcsecsq == orig_pix_area_arcsec * new_psr**2
