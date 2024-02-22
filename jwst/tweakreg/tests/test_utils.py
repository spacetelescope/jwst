from os import path
from copy import deepcopy

import pytest
import numpy as np
import asdf
from astropy.coordinates import SkyCoord
from astropy import units as u
from tweakwcs.correctors import JWSTWCSCorrector
from tweakwcs.linearfit import build_fit_matrix

from jwst.tweakreg.tests import data
from stdatamodels.jwst.datamodels import ImageModel
from jwst.tweakreg.utils import (
    adjust_wcs,
    transfer_wcs_correction,
    _wcsinfo_from_wcs_transform
)


data_path = path.split(path.abspath(data.__file__))[0]


def _get_test_pts(model, npts=5):
    assert npts >= 5
    ysize, xsize = model.data.shape
    x = np.empty(npts, dtype=float)
    y = np.empty(npts, dtype=float)
    x[:5] = [-0.499, -0.499, xsize - 0.501, xsize - 0.501, (xsize - 1.0) / 2.0]
    y[:5] = [-0.499, ysize - 0.501, -0.499, ysize - 0.501, (ysize - 1.0) / 2.0]
    if npts > 5:
        npts -= 5
        np.random.seed(0)
        x[5:] = xsize * np.random.random(npts) - 0.5
        y[5:] = ysize * np.random.random(npts) - 0.5
    return x, y


@pytest.fixture
def nircam_rate():
    wcs_file = path.join(data_path, 'nrcb1-wcs.asdf')
    with asdf.open(wcs_file) as af:
        wcs = deepcopy(af['wcs'])

    xsize = 204
    ysize = 204
    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.var_rnoise += 0

    wcsinfo = _wcsinfo_from_wcs_transform(wcs)
    im.meta.wcsinfo = {
        'ctype1': 'RA---TAN',
        'ctype2': 'DEC--TAN',
        'v2_ref': wcsinfo['v2_ref'],
        'v3_ref': wcsinfo['v3_ref'],
        'roll_ref': wcsinfo['roll_ref'],
        'ra_ref': wcsinfo['ra_ref'],
        'dec_ref': wcsinfo['dec_ref'],
        'v3yangle': -0.07385127,
        'vparity': -1,
        'wcsaxes': 2
    }

    im.meta.wcs = wcs

    im.meta.instrument = {
        'channel': 'LONG',
        'detector': 'NRCALONG',
        'filter': 'F444W',
        'lamp_mode': 'NONE',
        'module': 'A',
        'name': 'NIRCAM',
        'pupil': 'CLEAR'
    }
    im.meta.subarray = {
        'fastaxis': -1,
        'name': 'FULL',
        'slowaxis': 2,
        'xsize': xsize,
        'xstart': 1,
        'ysize': ysize,
        'ystart': 1
    }

    im.meta.observation = {
        'activity_id': '01',
        'date': '2021-10-25',
        'exposure_number': '00001',
        'obs_id': 'V42424001001P0000000001101',
        'observation_label': 'nircam_ptsrc_only',
        'observation_number': '001',
        'program_number': '42424',
        'sequence_id': '1',
        'time': '16:58:27.258',
        'visit_group': '01',
        'visit_id': '42424001001',
        'visit_number': '001'
    }

    im.meta.exposure = {
        'duration': 161.05155,
        'end_time': 59512.70899968495,
        'exposure_time': 150.31478,
        'frame_time': 10.73677,
        'group_time': 21.47354,
        'groupgap': 1,
        'integration_time': 150.31478,
        'mid_time': 59512.70812980775,
        'nframes': 1,
        'ngroups': 7,
        'nints': 1,
        'nresets_at_start': 1,
        'nresets_between_ints': 1,
        'readpatt': 'BRIGHT1',
        'sample_time': 10,
        'start_time': 59512.70725993055,
        'type': 'NRC_IMAGE'
    }

    im.meta.photometry = {
        'pixelarea_steradians': 1e-13,
        'pixelarea_arcsecsq': 4e-3,
    }

    return im


def test_adjust_wcs():
    wcs_file = path.join(data_path, 'nrcb1-wcs.asdf')
    with asdf.open(wcs_file) as af:
        w0 = deepcopy(af['wcs'])

    crval20, crval10 = w0.pipeline[-2].transform.angles_3.value.tolist()[-2:]
    crval10 = -crval10

    crpix10, crpix20 = w0.numerical_inverse(crval10, crval20)

    wa = adjust_wcs(
        w0,
        delta_ra=0.0135,
        delta_dec=0.0208,
        delta_roll=25.7,
        scale_factor=1.003
    )

    crval1a, crval2a = wa(crpix10, crpix20)

    assert np.allclose(
        [crval1a - crval10, crval2a - crval20],
        [0.0135, 0.0208],
        rtol=0,
        atol=1e-13
    )

    offset_ra10, offset_dec10 = w0(crpix10, crpix20 + 10)
    offset_ra1a, offset_dec1a = wa(crpix10, crpix20 + 10)

    ca0 = SkyCoord(crval1a * u.deg, crval2a * u.deg, frame='icrs')
    ca1 = SkyCoord(offset_ra1a * u.deg, offset_dec1a * u.deg, frame='icrs')
    c0 = SkyCoord(crval10 * u.deg, crval20 * u.deg, frame='icrs')
    c01 = SkyCoord(offset_ra10 * u.deg, offset_dec10 * u.deg, frame='icrs')

    sep0 = c0.separation(c01).degree
    sepa = ca0.separation(ca1).degree

    # test scale:
    assert np.allclose(sepa / sep0, 1.003, rtol=0, atol=1e-8)

    # test roll:
    assert np.allclose(
        ca1.position_angle(ca0).degree - c01.position_angle(c0).degree,
        25.7,
        rtol=1.0e-5,
        atol=0.0
    )


@pytest.mark.parametrize("input_type", ["wcs", "model", "file"])
def test_transfer_wcs_correction(nircam_rate, tmp_path, input_type):
    m1 = nircam_rate.copy()
    m2 = nircam_rate.copy()

    # apply some correction to m2 (from_image)
    wcsinfo = _wcsinfo_from_wcs_transform(m2.meta.wcs)
    corr = JWSTWCSCorrector(m2.meta.wcs, wcsinfo)
    mat = build_fit_matrix(23, 1.0001456)
    shift = corr.tanp_center_pixel_scale * np.array([3.5, 2.0])
    corr.set_correction(matrix=mat, shift=shift)
    m2.meta.wcs = corr.wcs

    if input_type == "file":
        # write to files:
        model_file1 = str(tmp_path / 'model1.fits')
        model_file2 = str(tmp_path / 'model2.fits')

        m1.write(model_file1)
        m2.write(model_file2)
        m1 = model_file1
        m2 = model_file2
        from_model = model_file2

    elif input_type == "wcs":
        from_model = m2.meta.wcs

    else:
        from_model = m2

    # transfer correction to m1:
    transfer_wcs_correction(m1, from_model)

    # test:
    if input_type == "file":
        m1 = ImageModel(m1)
        m2 = ImageModel(m2)

    x, y = _get_test_pts(m1, 5)

    assert np.allclose(
        np.linalg.norm(
            np.subtract(m2.meta.wcs(x, y), m1.meta.wcs(x, y)),
            axis=0
        ),
        0,
        rtol=0,
        atol=5e-11
    )
