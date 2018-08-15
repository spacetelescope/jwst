import os
import shutil
import tempfile

import pytest
from astropy.time import Time
import numpy as np
from numpy.testing import assert_allclose

from .. import (DataModel, ImageModel, MaskModel, QuadModel,
                MultiSlitModel, ModelContainer, SlitModel,
                SlitDataModel, IFUImageModel)
from ..util import open as open_model


ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')
FITS_FILE = os.path.join(ROOT_DIR, 'test.fits')
ASN_FILE = os.path.join(ROOT_DIR, 'association.json')


def reset_group_id(container):
    """Remove group_id from all models in container"""
    for m in container:
        try:
            del m.meta.group_id
        except AttributeError:
            pass


def setup():
    global FITS_FILE, TMP_DIR, TMP_FITS, TMP_YAML, TMP_JSON, TMP_FITS2

    TMP_DIR = tempfile.mkdtemp()
    TMP_FITS = os.path.join(TMP_DIR, 'tmp.fits')
    TMP_YAML = os.path.join(TMP_DIR, 'tmp.yaml')
    TMP_JSON = os.path.join(TMP_DIR, 'tmp.json')
    TMP_FITS2 = os.path.join(TMP_DIR, 'tmp2.fits')


def teardown():
    shutil.rmtree(TMP_DIR)


def test_set_shape():
    with pytest.raises(AttributeError):
        with ImageModel((50, 50)) as dm:
            assert dm.shape == (50, 50)
            dm.shape = (42, 23)


def test_broadcast():
    with ImageModel((50, 50)) as dm:
        data = np.empty((50,))
        dm.dq = data


def test_broadcast2():
    with ImageModel() as dm:
        data = np.empty((52, 50))
        dm.data = data

        dq = np.empty((50,))
        dm.dq = dq


def test_from_hdulist():
    from astropy.io import fits
    with fits.open(FITS_FILE) as hdulist:
        with open_model(hdulist) as dm:
            dm.data
        assert hdulist.fileinfo(0)['file'].closed == False


def delete_array():
    with open_model() as dm:
        del dm.data


def test_subarray():
    with DataModel(FITS_FILE) as dm:
        dm.meta.subarray.xstart


def roundtrip(func):
    def _create_source():
        dm = DataModel(FITS_FILE)

        assert dm.meta.instrument.name == 'MIRI'

        dm.meta.instrument.name = 'NIRCAM'
        dm.meta.subarray.xstart = 42
        return dm

    def _check_output(dm):
        assert dm.meta.instrument.name == 'NIRCAM'
        assert dm.meta.subarray.xstart == 42

    def test():
        with _create_source() as dm:
            with func(dm) as dm2:
                _check_output(dm2)

    test.__name__ = func.__name__

    return test


@roundtrip
def test_from_fits_write(dm):
    dm.to_fits(TMP_FITS, overwrite=True)
    return DataModel.from_fits(TMP_FITS)


def test_delete():
    with DataModel() as dm:
        dm.meta.instrument.name = 'NIRCAM'
        assert dm.meta.instrument.name == 'NIRCAM'
        del dm.meta.instrument.name
        assert dm.meta.instrument.name is None


def test_open():
    with open_model() as dm:
        pass

    with open_model((50, 50)) as dm:
        pass

    with open_model(FITS_FILE) as dm:
        assert isinstance(dm, QuadModel)


def test_copy():
    with ImageModel((50, 50)) as dm:
        dm.meta.instrument.name = "NIRCAM"
        dm.meta.foo = "BAR"

        with dm.copy() as dm2:
            dm2.data[0, 0] = 42
            assert np.sum(dm.data.flatten()) == 0

            assert dm2.meta.instrument.name == "NIRCAM"
            assert dm2.meta.foo == "BAR"
            dm2.meta.foo = "BAZ"
            assert dm.meta.foo == "BAR"
            dm2.meta.observation.obs_id = "FOO"
            assert dm.meta.observation.obs_id is None

def test_stringify():
    im = DataModel()
    assert str(im) == '<DataModel>'
    im.close()

    im = ImageModel((10,100))
    assert str(im) == '<ImageModel(10, 100)>'
    im.close()

    image = ROOT_DIR + "/nircam_mask.fits"
    im = MaskModel(image)
    assert str(im) ==  '<MaskModel(2048, 2048) from nircam_mask.fits>'
    im.close()

def test_section():
    with QuadModel((5, 35, 40, 32)) as dm:
        section = dm.get_section('data')[3:4, 1:3]
        assert section.shape == (1, 2, 40, 32)


def test_init_with_array():
    array = np.empty((50, 50))
    with open_model(array) as dm:
        assert dm.data.shape == (50, 50)
        assert isinstance(dm, ImageModel)


def test_init_with_array2():
    array = np.empty((50, 50))
    with ImageModel(array) as dm:
        assert dm.data.shape == (50, 50)


def test_init_with_array3():
    with pytest.raises(ValueError):
        array = np.empty((50,))
        with ImageModel(array) as dm:
            dm.data


def test_set_array():
    with pytest.raises(ValueError):
        with ImageModel() as dm:
            data = np.empty((50,))
            dm.data = data


def test_set_array2():
    with ImageModel() as dm:
        data = np.empty((50, 50))
        dm.data = data


def test_base_model_has_no_arrays():
    with pytest.raises(AttributeError):
        with DataModel() as dm:
            dm.data


def test_array_type():
    with ImageModel() as dm:
        assert dm.dq.dtype == np.uint32


def test_copy_model():
    with DataModel() as dm:
        with DataModel(dm) as dm2:
            assert hasattr(dm2, 'meta')


def test_dtype_match():
    with ImageModel() as dm:
        dm.data = np.array([[1, 2, 3]], np.float32)


def test_multislit():
    with MultiSlitModel() as dm:
        dm.slits.append(dm.slits.item())
        slit = dm.slits[-1]
        slit.data = np.random.rand(5, 5)
        slit.dm = np.random.rand(5, 5)
        slit.err = np.random.rand(5, 5)


def test_secondary_shapes():
    with ImageModel((256, 256)) as dm:
        assert dm.dq.shape == (256, 256)
        dm.dq = [1]


def test_multislit_from_image():
    with ImageModel((64, 64)) as im:
        with MultiSlitModel(im) as ms:
            assert len(ms.slits) == 1
            assert ms.slits[0].data.shape == (64, 64)


def test_multislit_from_fits_image():
    with ImageModel((64, 64)) as im:
        im.save(TMP_FITS, overwrite=True)

    with MultiSlitModel(TMP_FITS) as ms:
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
        im = ms.slits[0]
        im.name = "FOO"
        assert ms.slits[0].name == "FOO"


def test_multislit_metadata2():
    with MultiSlitModel() as ms:
        ms.slits.append(ms.slits.item())
        for key, val in ms.iteritems():
            assert isinstance(val, (bytes, str, int, float, bool, Time))


def test_multislit_copy():
    with MultiSlitModel() as input:
        for i in range(4):
            input.slits.append(input.slits.item(
                data=np.empty((50, 50), dtype=np.float32)))

        i = 0
        for slit in input.slits:
            i += 1
        assert i == 4

        input.save(TMP_FITS)

        output = input.copy()
        assert len(output.slits) == 4

        i = 0
        for slit in output.slits:
            i += 1
        assert i == 4

    from astropy.io import fits
    with fits.open(TMP_FITS) as hdulist:
        assert len(hdulist) == 6

    with MultiSlitModel(TMP_FITS) as input:
        i = 0
        for slit in input.slits:
            i += 1
        assert i == 4

        output = input.copy()
        assert len(output.slits) == 4


def test_model_with_nonstandard_primary_array():
    ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')

    class NonstandardPrimaryArrayModel(DataModel):
        schema_url = os.path.join(
            ROOT_DIR, "nonstandard_primary_array.schema.yaml")

        def __init__(self, init=None, wavelength=None, alpha=None, **kwargs):
            super(NonstandardPrimaryArrayModel, self).__init__(init=init, **kwargs)

            if wavelength is not None:
                self.wavelength = wavelength

            if alpha is not None:
                self.alpha = alpha

        # The wavelength array is the primary array.
        # Try commenting this function out and the problem goes away.
        def get_primary_array_name(self):
            return 'wavelength'

    m = NonstandardPrimaryArrayModel()
    list(m.keys())


def test_relsens():
    with ImageModel() as im:
        assert len(im.relsens.dtype) == 2


def test_image_with_extra_keyword_to_multislit():
    with ImageModel(data=np.empty((32, 32))) as im:
        im.save(TMP_FITS, overwrite=True)

    from astropy.io import fits
    with fits.open(TMP_FITS, mode="update") as hdulist:
        hdulist[1].header['BUNIT'] = 'x'

    with ImageModel(TMP_FITS) as im:
        with MultiSlitModel() as ms:
            ms.update(im)
            for i in range(3):
                ms.slits.append(ImageModel(data=np.empty((4, 4))))
            assert len(ms.slits) == 3

            ms.save(TMP_FITS2, overwrite=True)

    with MultiSlitModel(TMP_FITS2) as ms:
        assert len(ms.slits) == 3
        for slit in ms.slits:
            assert slit.data.shape == (4, 4)


@pytest.fixture(scope="module")
def container():
    with ModelContainer(ASN_FILE, persist=True) as c:
        for m in c:
            m.meta.observation.program_number = '0001'
            m.meta.observation.observation_number = '1'
            m.meta.observation.visit_number = '1'
            m.meta.observation.visit_group = '1'
            m.meta.observation.sequence_id = '01'
            m.meta.observation.activity_id = '1'
            m.meta.observation.exposure_number = '1'
            m.meta.instrument.name = 'NIRCAM'
            m.meta.instrument.channel = 'SHORT'
        yield c


def test_modelcontainer_iteration(container):
    for model in container:
        assert model.meta.telescope == 'JWST'


def test_modelcontainer_indexing(container):
    assert isinstance(container[0], DataModel)


def test_modelcontainer_group1(container):
    reset_group_id(container)
    for group in container.models_grouped:
        assert len(group) == 2
        for model in group:
            pass


def test_modelcontainer_group2(container):
    reset_group_id(container)
    container[0].meta.observation.exposure_number = '2'
    for group in container.models_grouped:
        assert len(group) == 1
        for model in group:
            pass
    container[0].meta.observation.exposure_number = '1'


def test_modelcontainer_group_names(container):
    reset_group_id(container)
    assert len(container.group_names) == 1
    reset_group_id(container)
    container[0].meta.observation.exposure_number = '2'
    assert len(container.group_names) == 2
    container[0].meta.observation.exposure_number = '1'


def test_object_node_iterator():
    im = ImageModel()
    items = []
    for i in im.meta.items():
        items.append(i[0])

    assert 'date' in items
    assert 'model_type' in items

def test_hasattr():
    model = DataModel()

    has_date = model.meta.hasattr('date')
    assert has_date, "Check that date exists"

    has_filename = model.meta.hasattr('filename')
    assert not has_filename, "Check that filename does not exist"

def test_multislit_model():
    data = np.arange(24, dtype=np.float32).reshape((6, 4))
    err = np.arange(24, dtype=np.float32).reshape((6, 4)) + 2
    wav = np.arange(24, dtype=np.float32).reshape((6, 4)) + 3
    dq = np.arange(24,dtype=np.uint32).reshape((6, 4)) + 1
    s0 = SlitDataModel(data=data, err=err, dq=dq, wavelength=wav)
    s1 = SlitDataModel(data=data+1, err=err+1, dq=dq+1, wavelength=wav+1)

    ms = MultiSlitModel()
    ms.slits.append(s0)
    ms.slits.append(s1)
    ms.meta.instrument.name = 'NIRSPEC'
    ms.meta.exposure.type = 'NRS_IMAGE'
    slit1 = ms[1]
    assert isinstance(slit1, SlitModel)
    assert slit1.meta.instrument.name == 'NIRSPEC'
    assert slit1.meta.exposure.type == 'NRS_IMAGE'
    assert_allclose(slit1.data, data + 1)


def test_slit_from_image():
    data = np.arange(24, dtype=np.float32).reshape((6, 4))
    im = ImageModel(data=data, err=data/2, dq=data)
    im.meta.instrument.name = "MIRI"
    slit_dm = SlitDataModel(im)
    assert_allclose(im.data, slit_dm.data)
    assert hasattr(slit_dm, 'pathloss_pointsource')
    # this should be enabled after gwcs starts using non-coordinate inputs
    #assert not hasattr(slit_dm, 'meta')

    slit = SlitModel(im)
    assert_allclose(im.data, slit.data)
    assert_allclose(im.err, slit.err)
    assert hasattr(slit, 'pathloss_pointsource')
    assert slit.meta.instrument.name == "MIRI"

    im = ImageModel(slit)
    assert type(im) == ImageModel
    im.close()

    im = ImageModel(slit_dm)
    assert type(im) == ImageModel
    im.close()


def test_ifuimage():
    data = np.arange(24, dtype=np.float32).reshape((6, 4))
    im = ImageModel(data=data, err=data/2, dq=data)
    ifuimage = IFUImageModel(im)
    assert_allclose(im.data, ifuimage.data)
    assert_allclose(im.err, ifuimage.err)
    assert_allclose(im.dq, ifuimage.dq)

    im = ImageModel(ifuimage)
    assert type(im) == ImageModel
    im.close()
