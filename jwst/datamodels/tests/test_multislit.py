from astropy.io import fits
from astropy.time import Time
import jsonschema
import numpy as np
from numpy.testing import assert_array_equal
import pytest

from jwst.datamodels import MultiSlitModel, ImageModel


def test_multislit_model():
    array1 = np.asarray(np.random.rand(2, 2), dtype='float32')
    array2 = np.asarray(np.random.rand(2, 2), dtype='float32')

    with MultiSlitModel() as ms:
        assert len(ms.slits) == 0
        ms.slits.append(ms.slits.item())
        ms.slits[-1].data = array1
        assert len(ms.slits) == 1
        ms.slits.append(ms.slits.item())
        ms.slits[-1].data = array2
        assert len(ms.slits) == 2
        for i, slit in enumerate(ms.slits):
            assert slit == ms.slits[i]
        ms2 = ms.copy()
        assert len(ms2.slits) == 2
        assert_array_equal(ms.slits[-1].data, array2)
        del ms.slits[0]
        assert len(ms.slits) == 1
        assert_array_equal(ms.slits[0].data, array2)


def test_multislit_move_from_fits(tmp_path):
    path = tmp_path / "multislit_move.fits"

    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    for i in range(5):
        hdu = fits.ImageHDU(data=np.zeros((64, 64)), name='SCI')
        hdu.ver = i + 1
        hdulist.append(hdu)

    hdulist.writeto(path)

    n = MultiSlitModel()
    with MultiSlitModel(path) as m:
        n.slits.append(m.slits[2])

        assert len(n.slits) == 1


def test_multislit_append_string():
    with pytest.raises(jsonschema.ValidationError):
        m = MultiSlitModel(strict_validation=True)
        m.slits.append('junk')


def test_multislit():
    with MultiSlitModel() as dm:
        dm.slits.append(dm.slits.item())
        slit = dm.slits[-1]
        slit.data = np.random.rand(5, 5)
        slit.dm = np.random.rand(5, 5)
        slit.err = np.random.rand(5, 5)
        assert slit.wavelength.shape == (0, 0)
        assert slit.pathloss_point.shape == (0, 0)
        assert slit.pathloss_uniform.shape == (0, 0)
        assert slit.barshadow.shape == (0, 0)


def test_multislit_from_image():
    with ImageModel((64, 64)) as im:
        with MultiSlitModel(im) as ms:
            assert len(ms.slits) == 1
            assert ms.slits[0].data.shape == (64, 64)


def test_multislit_from_saved_imagemodel(tmp_path):
    path = tmp_path / "multislit_from_image.fits"
    with ImageModel((64, 64)) as im:
        im.save(path)

    with MultiSlitModel(path) as ms:
        assert len(ms.slits) == 1
        assert ms.slits[0].data.shape == (64, 64)

        for i, slit in enumerate(ms.slits):
            assert slit.data is ms.slits[i].data

        ms2 = ms.copy()
        ms2.slits = ms.slits
        assert len(ms2.slits) == 1


def test_multislit_metadata():
    with MultiSlitModel() as ms:
        with ImageModel((64, 64)) as im:
            ms.slits.append(ms.slits.item())
            ms.slits[-1].data = im.data
        slit = ms.slits[0]
        slit.name = "FOO"
        assert ms.slits[0].name == "FOO"


def test_multislit_metadata2():
    with MultiSlitModel() as ms:
        ms.slits.append(ms.slits.item())
        for key, val in ms.items():
            assert isinstance(val, (bytes, str, int, float, bool, Time))


def test_multislit_copy(tmp_path):
    path = tmp_path / "multislit.fits"
    with MultiSlitModel() as input:
        for i in range(4):
            input.slits.append(input.slits.item(data=np.empty((50, 50))))

        assert len(input.slits) == 4
        input.save(path)

        output = input.copy()
        assert len(output.slits) == 4

    with fits.open(path, memmap=False) as hdulist:
        assert len(hdulist) == 6

    with MultiSlitModel(path) as model:
        for i, slit in enumerate(model.slits):
            pass
        assert i+1 == 4

        output = model.copy()
        assert len(output.slits) == 4


def test_copy_multslit():
    model1 = MultiSlitModel()
    model2 = MultiSlitModel()

    model1.slits.append(ImageModel(np.ones((1024, 1024))))
    model2.slits.append(ImageModel(np.ones((1024, 1024)) * 2))

    # Create the ouput model as a copy of the first input
    output = model1.copy()

    assert len(model1.slits) == 1
    assert len(model2.slits) == 1
    assert len(output.slits) == 1

    assert model1.slits[0].data[330, 330] == 1
    assert output.slits[0].data[330, 330] == 1
    assert id(model1.slits[0].data) != id(output.slits[0].data)

    output.slits[0].data = model1.slits[0].data - model2.slits[0].data

    assert model1.slits[0].data[330, 330] == 1
    assert output.slits[0].data[330, 330] == -1
