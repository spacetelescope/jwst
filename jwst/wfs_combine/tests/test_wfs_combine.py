import pytest

from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.associations.asn_from_list import asn_from_list
from jwst.wfs_combine import WfsCombineStep
from jwst.wfs_combine import wfs_combine
import numpy as np


GOOD = datamodels.dqflags.pixel["GOOD"]
DO_NOT_USE = datamodels.dqflags.pixel["DO_NOT_USE"]
SATURATED = datamodels.dqflags.pixel["SATURATED"]


def gaussian(x, mu, sig):
    return 1. / (np.sqrt(2 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2) / 2)


def gaussian2d(x, y, mux, muy, sigx, sigy):
    return gaussian(x, mux, sigx) * gaussian(y, muy, sigy)


def add_point_source(inarray, scale, xcen, ycen, sigx, sigy):
    outarray = np.zeros_like(inarray)
    for y in range(inarray.shape[0]):
        for x in range(inarray.shape[1]):
            outarray[y, x] = scale * gaussian2d(x, y, xcen, ycen, sigx, sigy)
#            outarray[y, x] = modeling.functional_models.Gaussian2D(amplitude=scale,
#                                                                   x_mean=xcen, y_mean=ycen,
#                                                                   x_stddev=sigx, y_stddev=sigy)
    return outarray


@pytest.fixture(scope="module")
def wfs_association(tmp_path_factory):
    imsize = 10
    tmp_path = tmp_path_factory.mktemp("wfs")
    im1 = datamodels.ImageModel((imsize, imsize))
    im1.meta.wcsinfo = {
        'dec_ref': 11.99875540218638,
        'ra_ref': 22.02351763251896,
        'roll_ref': 0.005076934167039675,
        'v2_ref': 86.039011,
        'v3_ref': -493.385704,
        'v3yangle': -0.07385127,
        'vparity': -1,
        'wcsaxes': 2}
    im1.meta.instrument = {
        'channel': 'SHORT',
        'detector': 'NRCA4',
        'filter': 'F212N',
        'module': 'A',
        'name': 'NIRCAM',
        'pupil': 'CLEAR'}
    im1.meta.observation = {
        'exposure_number': '1',
        'date': '2021-10-25',
        'time': '16:58:27.258'}
    im1.meta.exposure = {
        'type': 'NRC_IMAGE'}

    im1 = AssignWcsStep.call(im1, sip_approx=False)

    im2 = im1.copy()
    im2.meta.observation = {
        'exposure_number': '2',
        'date': '2021-10-25',
        'time': '17:58:27.258'}
    path1 = str(tmp_path / "image1_cal.fits")
    path2 = str(tmp_path / "image2_cal.fits")
    im1.save(path1)
    im2.save(path2)

    asn = asn_from_list([path1, path2],
                        product_name='jw00024-a3001_t001_nircam_nrca4_{suffix}')
    asn.data["program"] = "00024"
    asn.data["asn_type"] = "wfs-image2"
    asn.sequence = 1
    asn_name, serialized = asn.dump(format="json")
    path_asn = tmp_path / asn_name
    with open(path_asn, "w") as f:
        f.write(serialized)

    return path_asn, path1, path2


@pytest.mark.parametrize(
    "data1, data2, dq1, dq2, result_data, result_dq",
    [
        (1, 3, GOOD, GOOD, 2, GOOD),
        (1, 3, GOOD, DO_NOT_USE, 1, GOOD),
        (1, 3, DO_NOT_USE, GOOD, 3, GOOD),
        (1, 3, DO_NOT_USE, DO_NOT_USE, 0, DO_NOT_USE),
        (1, 3, GOOD, SATURATED, 2, GOOD),
        (1, 3, SATURATED, DO_NOT_USE, 1, SATURATED),
        (1, 3, GOOD, DO_NOT_USE + SATURATED, 1, GOOD),
        (1, 3, DO_NOT_USE, DO_NOT_USE + SATURATED, 0, DO_NOT_USE),
        (1, 3, SATURATED, SATURATED, 2, SATURATED),
    ]
)
def test_create_combined(_jail, wfs_association,
                         data1, data2, dq1, dq2, result_data, result_dq):
    path_asn, path1, path2 = wfs_association

    # Populate data and DQ array at [5, 5] with what is to be tested
    with datamodels.open(path1) as im1, datamodels.open(path2) as im2:
        # Populate some data and DQ flags for the test
        im1.data[5, 5] = data1
        im2.data[5, 5] = data2
        im1.dq[5, 5] = dq1
        im2.dq[5, 5] = dq2

        im1.save(path1)
        im2.save(path2)

    result = WfsCombineStep.call(path_asn, do_refine=False, flip_dithers=False, psf_size=50,
                                 blur_size=10, n_size=2)

    # Check that results are as expected
    assert result[0].data[5, 5] == result_data
    assert result[0].dq[5, 5] == result_dq


def test_shift_order_no_refine_no_flip(wfs_association):
    path_asn, path1, path2 = wfs_association
    nircam_pixel_size = 0.031 / 3600.
    delta_pixel = 5
    with datamodels.open(path2) as im2:
        im2.meta.wcsinfo = {
            'dec_ref': 11.99875540218638,
            'ra_ref': 22.02351763251896 + delta_pixel * nircam_pixel_size,
            'roll_ref': 0.005076934167039675, 'v2_ref': 86.039011, 'v3_ref': -493.385704,
            'v3yangle': -0.07385127, 'vparity': -1, 'wcsaxes': 2}
        im2 = AssignWcsStep.call(im2, sip_approx=False)

        im2.save(path2)
    wfs = wfs_combine.DataSet(path1, path2, 'outfile.fits', do_refine=False, flip_dithers=False, psf_size=50,
                              blur_size=10, n_size=2)
    wfs.do_all()
    assert wfs.input_1.meta.observation.exposure_number == '1'
    assert wfs.input_2.meta.observation.exposure_number == '2'
    assert wfs.off_x == abs(delta_pixel)
    assert wfs.off_y == 0
    wfs.input_1.close()
    wfs.input_2.close()


def test_shift_order_no_refine_with_flip(wfs_association):
    # change the sign of the the dither and check that we get the same delta pix but with switched
    # order of images.
    path_asn, path1, path2 = wfs_association
    nircam_pixel_size = 0.031 / 3600.
    delta_pixel = -5
    with datamodels.open(path2) as im2:
        im2.meta.wcsinfo = {
            'dec_ref': 11.99875540218638,
            'ra_ref': 22.02351763251896 + delta_pixel * nircam_pixel_size,
            'roll_ref': 0.005076934167039675, 'v2_ref': 86.039011, 'v3_ref': -493.385704,
            'v3yangle': -0.07385127, 'vparity': -1, 'wcsaxes': 2}
        im2 = AssignWcsStep.call(im2, sip_approx=False)
        im2.save(path2)
    wfs = wfs_combine.DataSet(path1, path2, 'outfile.fits', do_refine=False, flip_dithers=True, psf_size=50,
                              blur_size=10, n_size=2)
    wfs.do_all()
    wfs.input_1.close()
    wfs.input_2.close()
#    assert wfs.input_1.meta.observation.exposure_number == '2'
#    assert wfs.input_2.meta.observation.exposure_number == '1'
    assert wfs.off_x == -1 * delta_pixel


@pytest.mark.parametrize(
    "xshift, yshift, xerror, yerror, flip_dithers",
    [
        (5, 0, 0, 0, True),
        (-5, 0, 0, 0, True),
        (0, 6, 0, 0, True),
        (0, -6, 0, 0, True),
        (5, 5, 0, 0, True),
        (5, 5, 1, 1, True),
        (5, 5, -1, 1, True),
        (6, -5, 0, 0, True),
        (-5, 6, 0, 0, True),
        (5, -5, 1, -1, True),
        (-5, -5, 1, -1, True),
        (-5, 5, -1, 1, True),
        (5, 0, 0, 0, False),
        (-5, 0, 0, 0, False),
        (0, 6, 0, 0, False),
        (0, -6, 0, 0, False),
        (5, 5, 0, 0, False),
        (5, 5, 1, 1, False),
        (5, 5, -1, 1, False),
        (6, -5, 0, 0, False),
        (-5, 6, 0, 0, False),
        (5, -5, 1, -1, False),
        (-5, -5, 1, -1, False),
        (-5, 5, -1, 1, False),
    ]
)
def test_refine_no_error(wfs_association, xshift, yshift, xerror, yerror, flip_dithers):
    data_size = 201
    path_asn, path1, path2 = wfs_association
    nircam_pixel_size = 0.031 / 3600.
    delta_x_pixel = xshift
    # positive y pixel changes are negative declination changes
    delta_y_pixel = -1.0 * yshift
    with datamodels.open(path1) as im1:
        im1.data = np.zeros(shape=(data_size, data_size), dtype=np.float32)
        im1.dq = np.zeros(shape=(data_size, data_size), dtype=np.int32)
        im1.data = add_point_source(im1.data, 200, 100, 100, 4, 4)
        im1.save(path1)
    with datamodels.open(path2) as im2:
        im2.data = np.zeros(shape=(data_size, data_size), dtype=np.float32)
        im2.dq = np.zeros(shape=(data_size, data_size), dtype=np.int32)
        im2.meta.wcsinfo = {
            'dec_ref': 11.99875540218638 + delta_y_pixel * nircam_pixel_size,
            'ra_ref': 22.02351763251896 + delta_x_pixel * nircam_pixel_size,
            'roll_ref': 0.005076934167039675, 'v2_ref': 86.039011, 'v3_ref': -493.385704,
            'v3yangle': -0.07385127, 'vparity': -1, 'wcsaxes': 2}
        im2 = AssignWcsStep.call(im2, sip_approx=False)
        # shift the actual location of the 2nd image including the error in the WCS position
        im2.data = add_point_source(im2.data, 200, 100 + delta_x_pixel + xerror, 100 + delta_y_pixel - yerror, 4, 4)
        im2.save(path2)
    wfs = wfs_combine.DataSet(path1, path2, 'outfile.fits', do_refine=True, flip_dithers=flip_dithers, psf_size=50,
                              blur_size=10, n_size=2)
    wfs.do_all()
    assert np.max(abs(wfs.diff)) < 0.001
    if delta_x_pixel + xerror >= 0 or not flip_dithers:
        assert wfs.off_x == delta_x_pixel + xerror
        # Positive y pixel errors are negative in declination
        assert wfs.off_y == delta_y_pixel - yerror
    else:  # Order of exposures and sign of shifts have switched to always align the dithers
        assert wfs.off_x == abs(delta_x_pixel + xerror)
        # Positive y pixel errors are negative in declination
        assert wfs.off_y == -1.0 * (delta_y_pixel - yerror)
    wfs.input_1.close()
    wfs.input_2.close()


def test_refine_with_error(wfs_association):
    data_size = 201
    path_asn, path1, path2 = wfs_association
    nircam_pixel_size = 0.031 / 3600.
    shift_error = 2
    delta_pixel = 5
    with datamodels.open(path1) as im1:
        im1.data = np.zeros(shape=(data_size, data_size), dtype=np.float32)
        im1.data = add_point_source(im1.data, 200, 100, 100, 4, 4)
        im1.save(path1)
    with datamodels.open(path2) as im2:
        im2.meta.wcsinfo = {
            'dec_ref': 11.99875540218638,
            'ra_ref': 22.02351763251896 + delta_pixel * nircam_pixel_size,
            'roll_ref': 0.005076934167039675, 'v2_ref': 86.039011, 'v3_ref': -493.385704,
            'v3yangle': -0.07385127, 'vparity': -1, 'wcsaxes': 2}
        im2 = AssignWcsStep.call(im2, sip_approx=False)
        im2.data = np.zeros(shape=(data_size, data_size), dtype=np.float32)
        im2.data = add_point_source(im2.data, 200, 100 + delta_pixel + shift_error, 100, 4, 4)
        im2.save(path2)
    wfs = wfs_combine.DataSet(path1, path2, 'outfile.fits', do_refine=True, flip_dithers=True, psf_size=50,
                              blur_size=10, n_size=2)
    wfs.do_all()
    assert wfs.input_1.meta.observation.exposure_number == '1'
    assert wfs.input_2.meta.observation.exposure_number == '2'
    assert wfs.off_x == delta_pixel + shift_error
    assert wfs.off_y == 0
    wfs.input_1.close()
    wfs.input_2.close()
