"""Test various utility functions"""

from astropy import coordinates as coord
from astropy import wcs as fitswcs
from astropy.modeling import models as astmodels
from gwcs import coordinate_frames as cf
from gwcs.wcstools import wcs_from_fiducial
import numpy as np
import pytest

from stdatamodels.jwst.datamodels import SlitModel, dqflags

from jwst.resample.resample_spec import find_dispersion_axis


DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
GOOD = dqflags.pixel["GOOD"]

DQ = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
BITVALUES = 2**0 + 2**2
BITVALUES_STR = f"{2**0}, {2**2}"
BITVALUES_INV_STR = f"~{2**0}, {2**2}"
JWST_NAMES = "DO_NOT_USE, JUMP_DET"
JWST_NAMES_INV = "~" + JWST_NAMES


@pytest.fixture(scope="module")
def wcs_gwcs():
    crval = (150.0, 2.0)
    crpix = (500.0, 500.0)
    shape = (1000, 1000)
    pscale = 0.06 / 3600

    prj = astmodels.Pix2Sky_TAN()
    fiducial = np.array(crval)

    pc = np.array([[-1.0, 0.0], [0.0, 1.0]])
    pc_matrix = astmodels.AffineTransformation2D(pc, name="pc_rotation_matrix")
    scale = astmodels.Scale(pscale, name="cdelt1") & astmodels.Scale(pscale, name="cdelt2")
    transform = pc_matrix | scale

    out_frame = cf.CelestialFrame(
        name="world", axes_names=("lon", "lat"), reference_frame=coord.ICRS()
    )
    input_frame = cf.Frame2D(name="detector")
    wnew = wcs_from_fiducial(
        fiducial,
        coordinate_frame=out_frame,
        projection=prj,
        transform=transform,
        input_frame=input_frame,
    )

    output_bounding_box = ((0.0, float(shape[1])), (0.0, float(shape[0])))
    offset1, offset2 = crpix
    offsets = astmodels.Shift(-offset1, name="crpix1") & astmodels.Shift(-offset2, name="crpix2")

    wnew.insert_transform("detector", offsets, after=True)
    wnew.bounding_box = output_bounding_box

    tr = wnew.pipeline[0].transform
    pix_area = np.deg2rad(tr["cdelt1"].factor.value) * np.deg2rad(tr["cdelt2"].factor.value)

    wnew.pixel_area = pix_area
    wnew.pixel_shape = shape[::-1]
    wnew.array_shape = shape
    return wnew


@pytest.fixture(scope="module")
def wcs_fitswcs(wcs_gwcs):
    fits_wcs = fitswcs.WCS(wcs_gwcs.to_fits_sip())
    return fits_wcs


@pytest.fixture(scope="module")
def wcs_slicedwcs(wcs_gwcs):
    xmin, xmax = 100, 500
    slices = (slice(xmin, xmax), slice(xmin, xmax))
    sliced_wcs = fitswcs.wcsapi.SlicedLowLevelWCS(wcs_gwcs, slices)
    return sliced_wcs


def test_find_dispersion_axis():
    """
    Test the find_dispersion_axis() function
    """
    dm = SlitModel()

    dm.meta.wcsinfo.dispersion_direction = 1  # horizontal
    assert find_dispersion_axis(dm) == 0  # X axis for wcs functions

    dm.meta.wcsinfo.dispersion_direction = 2  # vertical
    assert find_dispersion_axis(dm) == 1  # Y axis for wcs functions
