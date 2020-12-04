import pytest

from jwst import datamodels
from jwst.associations.asn_from_list import asn_from_list
from jwst.wfs_combine.wfs_combine import DataSet
from jwst.wfs_combine import WfsCombineStep


@pytest.fixture
def wfs_association(tmp_path):
    im1 = datamodels.ImageModel((10, 10))
    im2 = datamodels.ImageModel((10, 10))
    im1.data += 10
    im2.data += 20

    path1 = str(tmp_path / "image1.fits")
    path2 = str(tmp_path / "image2.fits")
    im1.save(path1)
    im2.save(path2)

    asn = asn_from_list([path1, path2], product_name='wfs_combine_output')
    asn.data["program"] = "00024"
    asn.data["asn_type"] = "wfs-image2"
    from pprint import pprint
    print("")
    for product in asn.data["products"]:
        for member in product["members"]:
            pprint(member)
    asn.sequence = 1
    asn_name, serialized = asn.dump(format="json")
    path_asn = tmp_path / asn_name
    with open(path_asn, "w") as f:
        f.write(serialized)

    return path_asn


def test_create_combined(wfs_association):
    result = WfsCombineStep.call(wfs_association)
