"""
Test using a reference image in extract_1d for IFU data
"""
import numpy as np

from jwst import datamodels
from jwst.extract_1d import extract


data_shape = (941, 48, 46)
x_center = 23.
y_center = 26.                          # offset +2 from image center
radius = 11.5
inner_bkg = 11.5
outer_bkg = 16.5
method = "exact"

def test_ifu_3d():
    """Test 1"""

    input = make_ifu_cube(data_shape, source=5., background=3.7,
                          x_center=x_center, y_center=y_center,
                          radius=radius,
                          inner_bkg=inner_bkg, outer_bkg=outer_bkg)

    ref_image_2d = make_ref_image(data_shape[-2:],      # 2-D ref image
                                  x_center=x_center, y_center=y_center,
                                  radius=radius,
                                  inner_bkg=inner_bkg, outer_bkg=outer_bkg)

    ref_image_3d = make_ref_image(data_shape,           # 3-D ref image
                                  x_center=x_center, y_center=y_center,
                                  radius=radius,
                                  inner_bkg=inner_bkg, outer_bkg=outer_bkg)

    ref_dict_2d = {"ref_file_type": extract.FILE_TYPE_IMAGE,
                   "ref_model": ref_image_2d}
    truth = extract.do_extract1d(input, ref_dict_2d, smoothing_length=0,
                                 bkg_order=0, log_increment=50,
                                 subtract_background=True)

    ref_dict_3d = {"ref_file_type": extract.FILE_TYPE_IMAGE,
                   "ref_model": ref_image_3d}
    output = extract.do_extract1d(input, ref_dict_3d, smoothing_length=0,
                                  bkg_order=0, log_increment=50,
                                  subtract_background=True)

    true_wl = truth.spec[0].spec_table['wavelength']
    true_flux = truth.spec[0].spec_table['flux']
    true_bkg = truth.spec[0].spec_table['background']

    wavelength = output.spec[0].spec_table['wavelength']
    flux = output.spec[0].spec_table['flux']
    background = output.spec[0].spec_table['background']

    # These should all be the same because the reference image is the
    # same in every plane.

    assert np.allclose(wavelength, true_wl, rtol=1.e-14, atol=1.e-14)

    assert np.allclose(flux, true_flux, atol=1.e-14)

    assert np.allclose(background, true_bkg, atol=1.e-14)

    input.close()
    truth.close()
    output.close()
    del ref_dict_2d, ref_dict_3d
    ref_image_2d.close()
    ref_image_3d.close()


def make_ifu_cube(data_shape, source=None, background=None,
                  x_center=None, y_center=None,
                  radius=None, inner_bkg=None, outer_bkg=None):
    """Create "science" data for testing.

    Returns
    -------
    input_model : `~jwst.datamodels.ifucube.IFUCubeModel`
    """

    data = np.zeros(data_shape, dtype=np.float32)
    weightmap = np.zeros(data_shape, dtype=np.float32)
    r_2 = radius**2
    if inner_bkg is not None and outer_bkg is not None:
        create_background = True
        bkg = background
        inner_2 = inner_bkg**2
        outer_2 = outer_bkg**2
    elif inner_bkg is not None or outer_bkg is not None:
        raise RuntimeError("Specify both inner_bkg and outer_bkg or neither.")
    else:
        create_background = False
        bkg = 0.

    for j in range(data_shape[-2]):
        for i in range(data_shape[-1]):
            dist_2 = (float(i) - x_center)**2 + (float(j) - y_center)**2
            if dist_2 <= r_2:
                data[:, j, i] = source + bkg
                weightmap[:, j, i] = 1
            if create_background:
                if dist_2 > inner_2 and dist_2 <= outer_2:
                    data[:, j, i] = bkg
                    weightmap[:, j, i] = 1

    dq = np.zeros(data_shape, dtype=np.uint32)
    input_model = datamodels.IFUCubeModel(data=data, dq=dq, weightmap=weightmap)
    # Populate the BUNIT keyword so that in ifu.py the net will be moved
    # to the flux column.
    input_model.meta.bunit_data = 'MJy/sr'

    def mock_wcs(x, y, z):
        """Fake wcs method."""

        wavelength = np.linspace(0.5975, 5.2975, 941, endpoint=True,
                                 retstep=False, dtype=np.float64)

        if hasattr(z, 'dtype'):
            iz = np.around(z).astype(np.int64)
        else:
            iz = round(z)

        wl = wavelength[iz]
        ra = wl.copy()                          # dummy values
        dec = wl.copy()                         # dummy values
        return (ra, dec, wl)

    input_model.meta.wcs = mock_wcs
    input_model.meta.target.source_type = 'POINT'

    return input_model


# Functions user_bkg_spec_a, user_bkg_spec_b, and user_bkg_spec_c create
# "user-supplied" background data for the tests above.

def make_ref_image(shape,
                   x_center=None, y_center=None,
                   radius=None, inner_bkg=None, outer_bkg=None):

    """Create an image reference file for testing.

    Returns
    -------
    ref_image : `~jwst.datamodels.MultiExtract1dImageModel`
    """

    ref_image = datamodels.MultiExtract1dImageModel()

    if len(shape) < 2 or len(shape) > 3:
        raise RuntimeError("shape must be either 2-D or 3-D")

    mask = np.zeros(shape, dtype=np.float32)

    r_2 = radius**2
    if inner_bkg is not None and outer_bkg is not None:
        create_background = True
        inner_2 = inner_bkg**2
        outer_2 = outer_bkg**2
    elif inner_bkg is not None or outer_bkg is not None:
        raise RuntimeError("Specify both inner_bkg and outer_bkg or neither.")
    else:
        create_background = False

    for j in range(shape[-2]):
        for i in range(shape[-1]):
            dist_2 = (float(i) - x_center)**2 + (float(j) - y_center)**2
            if dist_2 <= r_2:
                mask[..., j, i] = 1.                    # source
            if create_background:
                if dist_2 > inner_2 and dist_2 <= outer_2:
                    mask[..., j, i] = -1.               # background

    ref_image = datamodels.MultiExtract1dImageModel()

    im = datamodels.Extract1dImageModel(data=mask)
    im.name = "ANY"
    ref_image.images.append(im)

    return ref_image
