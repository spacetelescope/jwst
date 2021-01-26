import pytest
import numpy as np

from jwst.assign_wcs import pointing
from jwst.datamodels.image import ImageModel

from .test_nircam import create_hdul


def create_imaging_datamodel(v2_ref, v3_ref, va_scale):
    hdul = create_hdul()
    datamodel = ImageModel(hdul)
    datamodel.meta.velocity_aberration.scale_factor = va_scale
    datamodel.meta.wcsinfo.v2_ref = v2_ref
    datamodel.meta.wcsinfo.v3_ref = v3_ref
    return datamodel


def test_va_corr_valid_args():
    v2_ref = -380
    v3_ref = -770
    va_scale = 1.001
    dm = create_imaging_datamodel(v2_ref=v2_ref, v3_ref=v3_ref, va_scale=va_scale)

    with pytest.raises(KeyError) as e:
        pointing.va_corr_model(None, v2_ref=v2_ref)
    assert str(e.value.args[0]) ==  ("'v2_ref', 'v3_ref', and 'va_scale' are all "
                                     "required when 'datamodel' is set to None.")

    with pytest.raises(TypeError) as e:
        pointing.va_corr_model(None, v2_ref='I', v3_ref='II', va_scale='III')
    assert str(e.value.args[0]) ==  "'v2_ref', 'v3_ref', and 'va_scale' must be numbers."

    with pytest.raises(ValueError) as e:
        pointing.va_corr_model(dm, v2_ref=v2_ref, v3_ref=v3_ref, va_scale=va_scale)
    assert str(e.value.args[0]) ==  ("'v2_ref', 'v3_ref', and 'va_scale' cannot be "
                                     "provided when 'datamodel' is not None.")


def test_va_corr_noop_missing_meta_values():
    dm = create_imaging_datamodel(v2_ref=None, v3_ref=None, va_scale=None)
    assert pointing.va_corr_model(dm) is None


@pytest.mark.parametrize(('fcorr'), [True, False])
def test_va_corr_valid_match(fcorr):
    v2_ref = -380
    v3_ref = -770
    va_scale = 1.001
    dm1 = create_imaging_datamodel(v2_ref=v2_ref, v3_ref=v3_ref, va_scale=va_scale)
    m1 = pointing.va_corr_model(dm1, fast_corr=fcorr)
    m2 = pointing.va_corr_model(None, fast_corr=fcorr, v2_ref=v2_ref, v3_ref=v3_ref, va_scale=va_scale)
    assert np.allclose(m1(1, 10), m2(1, 10))


@pytest.mark.parametrize(('fcorr'), [True, False])
def test_va_corr_inverse(fcorr):
    v2_ref = -380
    v3_ref = -770
    va_scale = 1.001
    test_v2 = 2300
    test_v3 = 5600
    fdm = create_imaging_datamodel(v2_ref=v2_ref, v3_ref=v3_ref, va_scale=va_scale)
    idm = create_imaging_datamodel(v2_ref=v2_ref, v3_ref=v3_ref, va_scale=1 / va_scale)
    fm = pointing.va_corr_model(fdm, fast_corr=fcorr)
    im = pointing.va_corr_model(idm, fast_corr=fcorr)
    assert np.allclose(im(*fm(test_v2, test_v3)), (test_v2, test_v3))


def test_va_corr_fast():
    v2_ref = 120.67
    v3_ref = -527.39
    va_scale = 1.00103
    v2 = 88
    v3 = -490

    m1 = pointing.va_corr_model(None, fast_corr=True, v2_ref=v2_ref, v3_ref=v3_ref, va_scale=va_scale)
    m2 = pointing.va_corr_model(None, fast_corr=False, v2_ref=v2_ref, v3_ref=v3_ref, va_scale=va_scale)

    assert np.allclose((v2, v3), m1.inverse(*m2(v2, v3)), rtol=0, atol=1e-7)
