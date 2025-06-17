"""
Test the utility functions
"""

from astropy.modeling.models import Shift, Identity
from astropy.table import QTable
from astropy.utils.data import get_pkg_data_filename
from stdatamodels.jwst import datamodels

from jwst.lib.catalog_utils import SkyObject

from jwst.assign_wcs.util import (
    get_object_info,
    wcs_bbox_from_shape,
    subarray_transform,
    bounding_box_from_subarray,
    transform_bbox_from_shape,
)


def test_transform_bbox_from_shape_2d():
    model = datamodels.ImageModel((512, 2048))
    bb = transform_bbox_from_shape(model.data.shape)
    assert bb == ((-0.5, 511.5), (-0.5, 2047.5))


def test_transform_bbox_from_shape_3d():
    model = datamodels.CubeModel((3, 32, 2048))
    bb = transform_bbox_from_shape(model.data.shape)
    assert bb == ((-0.5, 31.5), (-0.5, 2047.5))


def test_wcs_bbox_from_shape_2d():
    model = datamodels.ImageModel((512, 2048))
    bb = wcs_bbox_from_shape(model.data.shape)
    assert bb == ((-0.5, 2047.5), (-0.5, 511.5))


def test_wcs_bbox_from_shape_3d():
    model = datamodels.CubeModel((3, 32, 2048))
    bb = wcs_bbox_from_shape(model.data.shape)
    assert bb == ((-0.5, 2047.5), (-0.5, 31.5))

    model = datamodels.IFUCubeModel((750, 45, 50))
    bb = wcs_bbox_from_shape(model.data.shape)
    assert bb == ((-0.5, 49.5), (-0.5, 44.5))


def read_catalog(catalogname):
    return get_object_info(catalogname)


def test_create_grism_objects():
    source_catalog = get_pkg_data_filename(
        "data/step_SourceCatalogStep_cat.ecsv", package="jwst.assign_wcs.tests"
    )

    # create from test ascii file
    grism_objects = read_catalog(source_catalog)
    assert isinstance(grism_objects, list), "return grism objects were not a list"

    required_fields = list(SkyObject()._fields)
    go_fields = grism_objects[0]._fields
    assert all([a == b for a, b in zip(required_fields, go_fields)]), (
        "Required fields mismatch for SkyObject and GrismObject"
    )

    # create from QTable object
    tempcat = QTable.read(source_catalog, format="ascii.ecsv")
    grism_object_from_table = read_catalog(tempcat)
    assert isinstance(grism_object_from_table, list), "return grism objects were not a list"


def test_subarray_transform():
    im = datamodels.ImageModel()
    assert subarray_transform(im) is None

    im.meta.subarray.xstart = 3
    transform = subarray_transform(im)
    assert isinstance(transform[0], Shift) and transform[0].offset == 2
    assert isinstance(transform[1], Identity)

    im.meta.subarray.ystart = 5
    transform = subarray_transform(im)
    assert isinstance(transform[0], Shift) and transform[0].offset == 2
    assert isinstance(transform[1], Shift) and transform[1].offset == 4

    im = datamodels.ImageModel()
    im.meta.subarray.ystart = 5
    transform = subarray_transform(im)
    assert isinstance(transform[0], Identity)
    assert isinstance(transform[1], Shift) and transform[1].offset == 4


def test_bounding_box_from_subarray():
    im = datamodels.ImageModel()
    im.meta.subarray.xstart = 4
    im.meta.subarray.ystart = 6
    im.meta.subarray.xsize = 400
    im.meta.subarray.ysize = 600
    assert bounding_box_from_subarray(im) == ((-0.5, 599.5), (-0.5, 399.5))
