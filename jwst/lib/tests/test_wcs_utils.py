import numpy as np
from numpy.testing import assert_allclose

from astropy import units as u
from astropy import coordinates as coord
from astropy.modeling.models import Mapping, Identity, Shift, Scale
from gwcs import wcstools, wcs
from gwcs import coordinate_frames as cf

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.transforms.models import NirissSOSSModel
from jwst.lib.wcs_utils import get_wavelengths
from jwst.assign_wcs import util


def create_model():
    det = cf.Frame2D(name="detector", axes_order=(0, 1))

    sky = cf.CelestialFrame(name="sky", axes_order=(0, 1), reference_frame=coord.ICRS())
    slit_spatial = cf.Frame2D(
        name="slit_spatial",
        axes_order=(0, 1),
        unit=("", ""),
        axes_names=("x_slit", "y_slit"),
    )

    spec = cf.SpectralFrame(
        name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
    )
    slit_frame = cf.CompositeFrame([slit_spatial, spec], name="slit_frame")
    world = cf.CompositeFrame([sky, spec], name="world")

    det2slit = Mapping((0, 1, 1)) | (Identity(2) & (Scale(0.5) | Shift(0.5)))
    slit2sky = Identity(3)

    slit_wcs = wcs.WCS([(det, det2slit), (slit_frame, slit2sky), (world, None)])

    # compute wavelengths

    data = np.full((10, 10), fill_value=5.0)

    bounding_box = util.wcs_bbox_from_shape(data.shape)

    x, y = wcstools.grid_from_bounding_box(bounding_box, step=(1, 1))
    _, _, lam = slit_wcs(x, y)
    lam = lam.astype(np.float32)
    model = datamodels.SlitModel(data=data, wavelength=lam)
    model.meta.wcs = slit_wcs

    return model


def create_mock_wl():
    wl = np.arange(10.0)
    wl = wl[:, np.newaxis]
    wl = np.repeat(wl, 10, axis=1)
    wl = (wl * 0.5) + 0.5
    return wl


def test_get_wavelengths():
    # create a mock SlitModel
    model = create_model()

    # calculate what the wavelength array should be
    wl_og = create_mock_wl()

    # Test that the get wavelengths returns the wavelength grid
    wl = get_wavelengths(model)
    assert_allclose(wl, wl_og)

    del model.wavelength

    # Check that wavelengths can be generated from wcs when the
    # wavelength attribute is unavailable
    wl = get_wavelengths(model)
    assert_allclose(wl, wl_og)

    # Check that wavelengths are generated correctly when given a WFSS exp_type
    wl = get_wavelengths(model, exp_type="NRC_TSGRISM")
    assert_allclose(wl, wl_og)


def test_get_wavelengths_soss():
    # create a mock SlitModel
    model = create_model()

    del model.wavelength
    model.meta.exposure.type = "NIS_SOSS"

    wcs = model.meta.wcs
    new_wcs = NirissSOSSModel(
        [
            1,
        ],
        [
            wcs,
        ],
    )
    model.meta.wcs = new_wcs

    # calculate what the wavelength array should be
    wl_og = create_mock_wl()

    wl = get_wavelengths(model, order=1)
    assert_allclose(wl, wl_og)


def test_get_wavelength_wavecorr():
    # create a mock SlitModel
    model = create_model()

    wl_og = create_mock_wl()

    # Test use_wavecorr with no wavelength correction modificiation
    # get_wavelengths should return the same wavelengths for use_wavecorr
    # True and False

    wl_corr = get_wavelengths(model, use_wavecorr=True)
    assert_allclose(wl_corr, wl_og)

    wl_uncorr = get_wavelengths(model, use_wavecorr=False)
    assert_allclose(wl_corr, wl_uncorr)

    # Update the model wcs to add a wavelength corrected slit frame
    slit_spatial = cf.Frame2D(
        name="slit_spatial",
        axes_order=(0, 1),
        unit=("", ""),
        axes_names=("x_slit", "y_slit"),
    )
    spec = cf.SpectralFrame(
        name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
    )
    wcorr_frame = cf.CompositeFrame([slit_spatial, spec], name="wavecorr_frame")

    # Insert the new transform into the slit wcs object
    wave2wavecorr = Identity(2) & Shift(0.1)
    model.meta.wcs.insert_frame("slit_frame", wave2wavecorr, wcorr_frame)

    bounding_box = util.wcs_bbox_from_shape(model.data.shape)
    x, y = wcstools.grid_from_bounding_box(bounding_box, step=(1, 1))
    _, _, lam = model.meta.wcs(x, y)
    model.wavelength = lam

    # calculate what the corrected wavelength array should be
    wl_corr_og = wl_og + 0.1

    wl_corr = get_wavelengths(model, use_wavecorr=True)
    assert_allclose(wl_corr, wl_corr_og)

    wl_uncorr = get_wavelengths(model, use_wavecorr=False)
    assert_allclose(wl_uncorr, wl_og)
