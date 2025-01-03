import numpy as np
import pytest
from stdatamodels.jwst.datamodels import MiriLrsPsfModel


from jwst.extract_1d import psf_profile as pp


@pytest.mark.parametrize('exp_type', ['MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT', 'UNKNOWN'])
def test_open_psf(psf_reference_file, exp_type):
    # for any exptype, a model that can be read
    # as MiriLrsPsfModel will be, since it's the only
    # one implemented so far
    with pp.open_psf(psf_reference_file, exp_type=exp_type) as model:
        assert isinstance(model, MiriLrsPsfModel)


def test_open_psf_fail():
    with pytest.raises(NotImplementedError, match='could not be read'):
        pp.open_psf('bad_file', 'UNKNOWN')


@pytest.mark.parametrize('dispaxis', [1, 2])
def test_normalize_profile(nod_profile, dispaxis):
    profile = 2 * nod_profile
    if dispaxis == 2:
        profile = profile.T
    pp._normalize_profile(profile, dispaxis)
    assert np.allclose(np.sum(profile, axis=dispaxis - 1), 1.0)


@pytest.mark.parametrize('dispaxis', [1, 2])
def test_normalize_profile_with_nans(nod_profile, dispaxis):
    profile = -1 * nod_profile
    profile[10, :] = np.nan
    if dispaxis == 2:
        profile = profile.T

    pp._normalize_profile(profile, dispaxis)
    assert np.allclose(np.sum(profile, axis=dispaxis - 1), 1.0)
    assert np.all(np.isfinite(profile))


def test_make_cutout_profile(mock_miri_lrs_fs, psf_reference):
    data_shape = mock_miri_lrs_fs.data.shape
    yidx, xidx = np.mgrid[:data_shape[0], :data_shape[1]]
    dispaxis = 2

    psf_subpix = psf_reference.meta.psf.subpix
    profiles = pp._make_cutout_profile(xidx, yidx, psf_subpix, psf_reference.data, dispaxis)
    assert len(profiles) == 1
    assert profiles[0].shape == data_shape

    # No shift, profile is uniform and normalized to cross-dispersion size
    assert np.all(profiles[0] == 1 / data_shape[1])
