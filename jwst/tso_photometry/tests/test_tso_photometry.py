import math

import astropy.units as u
import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.lib import reffile_utils
from jwst.tso_photometry import tso_photometry as tp
from jwst.tso_photometry.tso_photometry_step import TSOPhotometryStep

# Default values for mock data
XCENTER = 75.0
YCENTER = 50.0
VALUE = 17.0
BACKGROUND = 0.8
RADIUS = 5.0
RADIUS_INNER = 8.0
RADIUS_OUTER = 11.0


def get_gain_2d(datamodel):
    """Get the 2D gain values from the gain reference file

    Parameters
    ----------
    datamodel : `CubeModel`
        The input data model for which the gain reference file is needed.

    Returns
    -------
    gain_2d : ndarray
        The 2D gain values from the gain reference file.
    """
    # Get the gain reference file
    gain_filename = TSOPhotometryStep().get_reference_file(datamodel, "gain")
    gain_m = datamodels.GainModel(gain_filename)

    # Get the relevant 2D gain values from the model
    if reffile_utils.ref_matches_sci(datamodel, gain_m):
        gain_2d = gain_m.data
    else:
        # If the reference file does not match the science data, extract
        # the subarray model that matches the science data.
        gain_2d = reffile_utils.get_subarray_model(datamodel, gain_m).data

    return gain_2d


def mk_data_array(shape, value, background, xcenter, ycenter, radius):
    """Create a data array with a flat circular source."""
    data = np.zeros(shape, dtype=np.float32) + background

    yidx, xidx = np.mgrid[: shape[1], : shape[2]]
    source_mask = ((xidx - xcenter) ** 2 + (yidx - ycenter) ** 2) < radius**2
    data[:, source_mask] += value

    return data


def make_mock_wcs(xcenter, ycenter):
    """Placeholder WCS"""

    def mock_wcs(x, y):
        if np.isscalar(x):
            x = [x]
            y = [y]
        crpix1 = xcenter
        crpix2 = ycenter
        cdelt1 = 0.031
        cdelt2 = 0.031
        crval1 = 15.0 * (((34.4147 / 60.0 + 48.0) / 60.0) + 15.0)  # 15h 48m 34.4147s
        crval2 = ((24.295 / 60.0 + 9.0) / 60.0) + 28.0  # 28d 09' 24.295"

        ra, dec = [], []
        for xval, yval in zip(x, y):
            decval = (yval + 1.0 - crpix2) * cdelt2 + crval2
            dec.append(decval)
            cosdec = math.cos(decval * math.pi / 180.0)
            ra.append((xval + 1.0 - crpix1) * cdelt1 / cosdec + crval1)

        if np.isscalar(x):
            return ra[0], dec[0]
        else:
            return ra, dec

    return mock_wcs


def set_meta(datamodel, xcenter, ycenter, sub64p=False):
    """Assign some metadata"""

    datamodel.meta.exposure.nints = datamodel.data.shape[0]
    datamodel.meta.exposure.integration_time = 10.7
    datamodel.meta.exposure.start_time = 58704.62847
    datamodel.meta.observation.date = "2019-08-09"
    datamodel.meta.observation.time = "15:04:59.808"

    datamodel.meta.exposure.integration_start = 5
    datamodel.meta.exposure.integration_end = (
        datamodel.meta.exposure.integration_start + datamodel.meta.exposure.nints - 1
    )

    datamodel.meta.instrument.name = "NIRCAM"
    datamodel.meta.instrument.detector = "NRCA2"
    datamodel.meta.instrument.channel = "SHORT"
    datamodel.meta.instrument.filter = "F210M"
    datamodel.meta.target.catalog_name = "R Coronae Borealis"
    if sub64p:
        datamodel.meta.instrument.pupil = "WLP8"
        datamodel.meta.subarray.name = "SUB64P"
        datamodel.meta.subarray.xstart = 1
        datamodel.meta.subarray.xsize = datamodel.data.shape[2]
        datamodel.meta.subarray.ystart = 1
        datamodel.meta.subarray.ysize = datamodel.data.shape[1]
    else:
        datamodel.meta.instrument.pupil = "CLEAR"
        datamodel.meta.subarray.name = "FULL"
        datamodel.meta.subarray.xstart = 1
        datamodel.meta.subarray.xsize = datamodel.data.shape[2]
        datamodel.meta.subarray.ystart = 1
        datamodel.meta.subarray.ysize = datamodel.data.shape[1]

    datamodel.meta.bunit_data = "MJy/sr"
    datamodel.meta.bunit_err = "MJy/sr"
    # NOTE: this is a placeholder value that leaves the mock test data values
    # unchanged during the unit conversion in tso_photometry.
    datamodel.meta.photometry.pixelarea_steradians = 1.0e-6

    datamodel.meta.wcs = make_mock_wcs(xcenter, ycenter)
    datamodel.meta.wcsinfo.siaf_xref_sci = xcenter + 1
    datamodel.meta.wcsinfo.siaf_yref_sci = ycenter + 1


def include_int_times(datamodel):
    """Create an int_times table and copy it into the data model."""

    int_start = datamodel.meta.exposure.integration_start
    shape = datamodel.data.shape
    nrows = int_start + shape[0] + 2  # create a few extra rows
    # integration_number and time_arr are one_indexed
    integration_number = np.arange(1, nrows + 1, dtype=np.float32)
    time_arr = np.arange(1, nrows + 1, dtype=np.float64) + 58700.0
    bjd_arr = time_arr + 0.5  # placeholder BJD TDB values
    mock_data = np.arange(nrows, dtype=np.float64)

    it_dtype = [
        ("integration_number", "<i4"),
        ("int_start_MJD_UTC", "<f8"),
        ("int_mid_MJD_UTC", "<f8"),
        ("int_end_MJD_UTC", "<f8"),
        ("int_start_BJD_TDB", "<f8"),
        ("int_mid_BJD_TDB", "<f8"),
        ("int_end_BJD_TDB", "<f8"),
    ]

    otab = np.array(
        list(
            zip(integration_number, mock_data, time_arr, mock_data, mock_data, bjd_arr, mock_data)
        ),
        dtype=it_dtype,
    )
    datamodel.int_times = otab.copy()

    return time_arr, bjd_arr


def mock_nircam_image(
    shape=(7, 100, 150),
    xcenter=XCENTER,
    ycenter=YCENTER,
    value=VALUE,
    background=BACKGROUND,
    radius=RADIUS,
    sub64p=False,
    convert_units=False,
):
    data = mk_data_array(shape, value, background, xcenter, ycenter, radius)
    datamodel = datamodels.CubeModel(data)
    set_meta(datamodel, xcenter, ycenter, sub64p=sub64p)
    if convert_units:
        tp.convert_data_units(datamodel)
    return datamodel


def test_tso_phot():
    datamodel = mock_nircam_image(convert_units=True)

    # add some NaNs to one of the integrations
    # make sure at least one is in the annulus and one in the aperture
    data = datamodel.data
    data[-1, 50, 75] = np.nan
    data[-1, 50 + 2, 75 - 2] = np.nan
    data[-1, 50 - 9, 75] = np.nan
    datamodel.data = data

    # Use a larger radius than was used for creating the data.
    catalog = tp.tso_aperture_photometry(
        datamodel,
        XCENTER,
        YCENTER,
        RADIUS + 1.0,
        RADIUS_INNER,
        RADIUS_OUTER,
    )

    assert catalog.meta["instrument"] == datamodel.meta.instrument.name
    assert catalog.meta["filter"] == datamodel.meta.instrument.filter
    assert catalog.meta["subarray"] == datamodel.meta.subarray.name
    assert catalog.meta["detector"] == datamodel.meta.instrument.detector
    assert catalog.meta["channel"] == datamodel.meta.instrument.channel
    assert catalog.meta["target_name"] == datamodel.meta.target.catalog_name

    assert np.allclose(catalog["aperture_x"].value, XCENTER, atol=0.01)
    assert np.allclose(catalog["aperture_y"].value, YCENTER, atol=0.01)

    assert np.allclose(catalog["aperture_sum"].value[:-1], 1263.4778, rtol=1.0e-7)
    assert np.allclose(catalog["net_aperture_sum"].value[:-1], 1173.0, rtol=1.0e-7)
    assert np.allclose(catalog["annulus_sum"].value[:-1], 143.256627, rtol=1.0e-7)
    # check that NaNs made sums smaller
    assert np.isclose(catalog["aperture_sum"].value[-1], 1227.8778178321238, rtol=1.0e-7)
    assert np.isclose(catalog["net_aperture_sum"].value[-1], 1138.9999480843544, rtol=1.0e-7)
    assert np.isclose(catalog["annulus_sum"].value[-1], 142.45662712646373, rtol=1.0e-7)

    # mean of background annulus should be the same even with NaNs because dividing by non-NaN area
    assert np.allclose(catalog["annulus_mean"].value, BACKGROUND, rtol=1.0e-7)
    assert np.allclose(catalog["aperture_sum_err"].value, 0.0, atol=1.0e-7)
    assert np.allclose(catalog["annulus_sum_err"].value, 0.0, atol=1.0e-7)
    assert np.allclose(catalog["annulus_mean_err"].value, 0.0, rtol=1.0e-7)
    assert np.allclose(catalog["net_aperture_sum_err"].value, 0.0, atol=1.0e-7)


def test_tso_phot_sub64p():
    """Test with pupil = 'WLP8' and subarray = 'SUB64P'."""

    shape = (7, 64, 64)
    xcenter = 31.0
    ycenter = 31.0
    radius = 50.0
    radius_inner = None
    radius_outer = None

    datamodel = mock_nircam_image(
        shape=shape,
        xcenter=xcenter,
        ycenter=ycenter,
        radius=radius,
        sub64p=True,
        convert_units=True,
    )

    catalog = tp.tso_aperture_photometry(
        datamodel, xcenter, ycenter, radius, radius_inner, radius_outer
    )

    assert catalog.meta["instrument"] == datamodel.meta.instrument.name
    assert catalog.meta["filter"] == datamodel.meta.instrument.filter
    assert catalog.meta["subarray"] == datamodel.meta.subarray.name
    assert catalog.meta["pupil"] == datamodel.meta.instrument.pupil

    assert np.allclose(catalog["aperture_x"].value, xcenter, atol=0.01)
    assert np.allclose(catalog["aperture_y"].value, ycenter, atol=0.01)
    assert np.allclose(
        catalog["aperture_sum"].value, (VALUE + BACKGROUND) * shape[-1] * shape[-2], rtol=1.0e-6
    )
    assert np.allclose(catalog["aperture_sum_err"].value, 0.0, atol=1.0e-7)


def test_tso_phot_with_int_times():
    datamodel = mock_nircam_image()

    # Add integration times to the model
    int_times, bjd_times = include_int_times(datamodel)

    catalog = tp.tso_aperture_photometry(
        datamodel, XCENTER, YCENTER, RADIUS + 1.0, RADIUS_INNER, RADIUS_OUTER
    )

    offset = datamodel.meta.exposure.integration_start - 1
    slc = slice(offset, offset + datamodel.data.shape[0])
    assert np.allclose(catalog["MJD"], int_times[slc], atol=1.0e-8)
    assert np.allclose(catalog["BJD_TDB"], bjd_times[slc], atol=1.0e-8)


def test_tso_phot_int_times_out_of_range():
    datamodel = mock_nircam_image()

    # Add integration times
    include_int_times(datamodel)

    # Modify the column of integration numbers so that they extend outside
    # the range of integration numbers in the science data.  This shouldn't
    # happen, i.e. this is to test an edge case.
    int_start = datamodel.meta.exposure.integration_start
    datamodel.int_times["integration_number"] += 2 * int_start

    catalog = tp.tso_aperture_photometry(
        datamodel, XCENTER, YCENTER, RADIUS + 1.0, RADIUS_INNER, RADIUS_OUTER
    )

    int_times = np.array(
        [
            58704.62853,
            58704.628655,
            58704.62878,
            58704.62890,
            58704.6290,
            58704.6291511,
            58704.629275,
        ]
    )
    assert np.allclose(catalog["MJD"], int_times, rtol=1.0e-8)
    assert np.all(np.isnan(catalog["BJD_TDB"]))


def test_tso_phot_missing_int_start():
    datamodel = mock_nircam_image()

    # Add integration times
    int_times, bjd_times = include_int_times(datamodel)

    # Remove the integration start: it sill be assumed to be 1
    datamodel.meta.exposure.integration_start = None

    catalog = tp.tso_aperture_photometry(
        datamodel, XCENTER, YCENTER, RADIUS + 1.0, RADIUS_INNER, RADIUS_OUTER
    )

    offset = 0
    slc = slice(offset, offset + datamodel.data.shape[0])
    assert np.allclose(catalog["MJD"], int_times[slc], atol=1.0e-8)
    assert np.allclose(catalog["BJD_TDB"], bjd_times[slc], atol=1.0e-8)


def test_tso_phot_uncalibrated():
    # pupil = 'WLP8' and subarray = 'SUB64P'
    shape = (7, 64, 64)
    xcenter = 31.0
    ycenter = 31.0
    radius = 50.0
    radius_inner = None
    radius_outer = None
    datamodel = mock_nircam_image(
        shape=shape,
        xcenter=xcenter,
        ycenter=ycenter,
        radius=radius,
        sub64p=True,
    )

    # Set uncalibrated data units
    datamodel.meta.bunit_data = "DN/s"
    datamodel.meta.bunit_err = "DN/s"

    # Get the gain reference file
    gain_2d = get_gain_2d(datamodel)

    # Convert the data units
    tp.convert_data_units(datamodel, gain_2d)

    catalog = tp.tso_aperture_photometry(
        datamodel, xcenter, ycenter, radius, radius_inner, radius_outer
    )

    assert catalog.meta["instrument"] == datamodel.meta.instrument.name
    assert catalog.meta["filter"] == datamodel.meta.instrument.filter
    assert catalog.meta["subarray"] == datamodel.meta.subarray.name
    assert catalog.meta["pupil"] == datamodel.meta.instrument.pupil
    assert np.allclose(catalog["aperture_x"].value, xcenter, atol=0.01)
    assert np.allclose(catalog["aperture_y"].value, ycenter, atol=0.01)
    assert catalog["aperture_sum"][0].unit == u.Unit("electron")


def test_tso_phot_unexpected_units():
    datamodel = mock_nircam_image()

    # Set unknown data units
    datamodel.meta.bunit_data = "DN"
    datamodel.meta.bunit_err = "DN"

    # Get the gain reference file
    gain_2d = get_gain_2d(datamodel)
    tp.convert_data_units(datamodel, gain_2d)

    catalog = tp.tso_aperture_photometry(
        datamodel, XCENTER, YCENTER, RADIUS + 1.0, RADIUS_INNER, RADIUS_OUTER
    )

    # Unexpected units are left alone
    assert np.allclose(catalog["aperture_sum"].value, 1263.4778, rtol=1.0e-7)
    assert catalog["aperture_sum"][0].unit == u.Unit("DN")


def test_tso_phot_wrong_model():
    model = datamodels.ImageModel()
    with pytest.raises(TypeError, match="must be a CubeModel"):
        tp.tso_aperture_photometry(
            model, XCENTER, YCENTER, RADIUS + 1.0, RADIUS_INNER, RADIUS_OUTER, None
        )


def test_tso_phot_multiple_center_values():
    datamodel = mock_nircam_image()

    # Use an array of xcenter, ycenter values and provide psf values
    nint = datamodel.data.shape[0]
    xcenter = np.full(nint, XCENTER)
    ycenter = np.full(nint, YCENTER)
    centroid_x = np.arange(nint)
    centroid_y = np.arange(nint)
    psf_width_x = np.full(nint, 1.0)
    psf_width_y = np.full(nint, 2.0)
    psf_flux = np.full(nint, 3.0)

    catalog = tp.tso_aperture_photometry(
        datamodel,
        xcenter,
        ycenter,
        RADIUS + 1.0,
        RADIUS_INNER,
        RADIUS_OUTER,
        centroid_x=centroid_x,
        centroid_y=centroid_y,
        psf_width_x=psf_width_x,
        psf_width_y=psf_width_y,
        psf_flux=psf_flux,
    )

    assert np.allclose(catalog["aperture_x"].value, XCENTER, atol=0.01)
    assert np.allclose(catalog["aperture_y"].value, YCENTER, atol=0.01)
    assert np.allclose(catalog["aperture_sum"].value, 1263.4778, rtol=1.0e-7)
    assert np.allclose(catalog["aperture_x"].value, XCENTER, atol=0.01)
    assert np.allclose(catalog["aperture_y"].value, YCENTER, atol=0.01)
    assert np.allclose(catalog["centroid_x"].value, np.arange(nint))
    assert np.allclose(catalog["centroid_y"].value, np.arange(nint))
    assert np.allclose(catalog["psf_width_x"].value, 1.0)
    assert np.allclose(catalog["psf_width_y"].value, 2.0)
    assert np.allclose(catalog["psf_flux"].value, 3.0)


@pytest.mark.parametrize("fit_psf", [True, False])
def test_fit_source_fail(monkeypatch, fit_psf):
    datamodel = mock_nircam_image()
    mask = np.full(datamodel.data.shape, False)
    box_size = int(RADIUS * 2 + 1)
    xcenter, ycenter = XCENTER, YCENTER

    def mock_centroid(*args, **kwargs):
        raise ValueError("test fail")

    monkeypatch.setattr(tp, "centroid_sources", mock_centroid)

    # Failure in centroid just returns NaNs for all values
    result = tp._fit_source(
        datamodel.data, mask, mask[0], xcenter, ycenter, box_size, fit_psf=fit_psf
    )
    if fit_psf:
        assert len(result) == 5
        assert np.all(np.isnan(result))
    else:
        assert len(result) == 2
        assert np.all(np.isnan(result))


@pytest.mark.parametrize("centroid_values", [([-1], [-1]), ([0], [-1]), ([-1], [0])])
@pytest.mark.parametrize("fit_psf", [True, False])
def test_fit_source_centroid_out_of_bounds(monkeypatch, fit_psf, centroid_values):
    datamodel = mock_nircam_image()
    mask = np.full(datamodel.data.shape, False)
    box_size = int(RADIUS * 2 + 1)
    xcenter, ycenter = XCENTER, YCENTER

    def mock_centroid(*args, **kwargs):
        return centroid_values

    monkeypatch.setattr(tp, "centroid_sources", mock_centroid)

    # Failure in centroid just returns NaNs for all values
    result = tp._fit_source(
        datamodel.data, mask, mask[0], xcenter, ycenter, box_size, fit_psf=fit_psf
    )
    if fit_psf:
        assert len(result) == 5
        assert np.all(np.isnan(result))
    else:
        assert len(result) == 2
        assert np.all(np.isnan(result))


def test_fit_source_psf_fail(monkeypatch):
    datamodel = mock_nircam_image()
    mask = np.full(datamodel.data.shape, False)
    box_size = int(RADIUS * 2 + 1)
    xcenter, ycenter = XCENTER, YCENTER

    def mock_psf(*args, **kwargs):
        raise ValueError("test fail")

    monkeypatch.setattr(tp, "_psf_fit_gaussian_prf", mock_psf)

    # Failure in fit returns valid centroids, but NaNs for all psf values
    result = tp._fit_source(datamodel.data, mask, mask[0], xcenter, ycenter, box_size, fit_psf=True)
    assert len(result) == 5
    assert not np.all(np.isnan(result[:3]))
    assert np.all(np.isnan(result[3:]))


def test_psf_fit():
    datamodel = mock_nircam_image()
    data = datamodel.data[0]
    mask = np.full(data.shape, False)
    fit_box_width = int(RADIUS * 2 + 1)
    xcenter, ycenter = XCENTER, YCENTER

    # Fit values vary for the flat synthetic source, but they should
    # be close to the initial estimates
    x_width, y_width, flux = tp._psf_fit_gaussian_prf(data, mask, fit_box_width, xcenter, ycenter)
    assert np.allclose([x_width, y_width], RADIUS / 2, rtol=0.3)
    assert np.allclose(flux, 1263.4778, rtol=0.3)  # expected aperture sum value


def test_psf_fit_gaussian_prf_fail(monkeypatch):
    datamodel = mock_nircam_image()
    data = datamodel.data[0]
    mask = np.full(data.shape, False)
    fit_box_width = int(RADIUS * 2 + 1)
    xcenter, ycenter = XCENTER, YCENTER

    def mock_phot_call(*args, **kwargs):
        raise ValueError("test fail")

    monkeypatch.setattr(tp.PSFPhotometry, "__call__", mock_phot_call)

    # Failure in psf fit just returns NaNs
    result = tp._psf_fit_gaussian_prf(data, mask, fit_box_width, xcenter, ycenter)
    assert len(result) == 3
    assert np.all(np.isnan(result))


def test_tso_source_centroid_fail(monkeypatch):
    datamodel = mock_nircam_image()
    xcenter, ycenter = XCENTER, YCENTER

    def mock_fit_source(*args, **kwargs):
        nan_array = np.full(datamodel.shape[0], np.nan)
        return nan_array, nan_array

    monkeypatch.setattr(tp, "_fit_source", mock_fit_source)

    # Failure in fit just returns NaNs for centroid,
    # None for PSF values
    result = tp.tso_source_centroid(datamodel, xcenter, ycenter)
    assert len(result) == 5
    assert np.all(np.isnan(result[:2]))
    assert result[2] is None
    assert result[3] is None
    assert result[4] is None
