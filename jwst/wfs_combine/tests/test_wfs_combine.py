import pytest

from jwst import datamodels
from jwst.assign_wcs import AssignWcsStep
from jwst.associations.asn_from_list import asn_from_list
from jwst.wfs_combine import WfsCombineStep


GOOD = datamodels.dqflags.pixel["GOOD"]
DO_NOT_USE = datamodels.dqflags.pixel["DO_NOT_USE"]
SATURATED = datamodels.dqflags.pixel["SATURATED"]


@pytest.fixture(scope="module")
def wfs_association(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("wfs")
    im1 = datamodels.ImageModel((10, 10))

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
        'date': '2021-10-25',
        'time': '16:58:27.258'}
    im1.meta.exposure = {
        'type': 'NRC_IMAGE'}

    im1 = AssignWcsStep.call(im1, sip_approx=False)

    im2 = im1.copy()

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

    result = WfsCombineStep.call(path_asn, do_refine=False)

    # Check that results are as expected
    assert result.data[5, 5] == result_data
    assert result.dq[5, 5] == result_dq
