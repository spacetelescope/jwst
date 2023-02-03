import pytest

from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.associations import asn_from_list
from jwst.wfs_combine import WfsCombineStep


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

    asn = asn_from_list.asn_from_list([path1, path2], product_name="combined",
                                      rule=DMS_Level3_Base)
    asn.data["program"] = "00024"
    asn.data["asn_type"] = "wfs-image2"
    asn.sequence = 1
    asn_json, serialized = asn.dump(format="json")
    path_asn = tmp_path / asn_json
    with open(path_asn, "w") as f:
        f.write(serialized)

    return path_asn, path1, path2


def test_step_pos_shift_no_refine_no_flip(wfs_association):
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

    wfs = WfsCombineStep.call(path_asn, do_refine=False, flip_dithers=False, psf_size=50,
                              blur_size=10, n_size=2)
    assert wfs[0].meta.wcsinfo.ra_ref == 22.02351763251896


def test_step_neg_shift_no_refine_no_flip(wfs_association):
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

    wfs = WfsCombineStep.call(path_asn, do_refine=False, flip_dithers=False, psf_size=50,
                              blur_size=10, n_size=2)
    assert wfs[0].meta.wcsinfo.ra_ref == 22.02351763251896


def test_step_neg_order_no_refine_with_flip(wfs_association):
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

    wfs = WfsCombineStep.call(path_asn, do_refine=False, flip_dithers=True, psf_size=50,
                              blur_size=10, n_size=2)
    assert wfs[0].meta.wcsinfo.ra_ref == 22.02351763251896 + delta_pixel * nircam_pixel_size


def test_step_pos_order_no_refine_with_flip(wfs_association):
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

    wfs = WfsCombineStep.call(path_asn, do_refine=False, flip_dithers=True, psf_size=50,
                              blur_size=10, n_size=2)
    assert wfs[0].meta.wcsinfo.ra_ref == 22.02351763251896
