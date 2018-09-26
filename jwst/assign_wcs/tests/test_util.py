from jwst import datamodels
from jwst.assign_wcs.util import bounding_box_from_shape


def test_bounding_box_from_shape_2d():
    model = datamodels.ImageModel((512, 2048))
    bb = bounding_box_from_shape(model.data.shape)
    assert bb == ((-0.5, 2047.5), (-0.5, 511.5))


def test_bounding_box_from_shape_3d():
    model = datamodels.CubeModel((3, 32, 2048))
    bb = bounding_box_from_shape(model.data.shape)
    assert bb == ((-0.5, 2047.5), (-0.5, 31.5))

    model = datamodels.IFUCubeModel((750, 45, 50))
    bb = bounding_box_from_shape(model.data.shape)
    assert bb == ((-0.5, 49.5), (-0.5, 44.5))