import os.path
import numpy as np
from numpy.testing import assert_allclose
import pytest

from gwcs import wcstools

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.transforms import models

import jwst
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d import Extract2dStep
from jwst.srctype import SourceTypeStep
from jwst.wavecorr import WavecorrStep
from jwst.wavecorr import wavecorr

from jwst.assign_wcs.tests.test_nirspec import (create_nirspec_mos_file,
                                                create_nirspec_fs_file)


def test_wavecorr():
    hdul = create_nirspec_mos_file()
    msa_meta = os.path.join(jwst.__path__[0], *['assign_wcs', 'tests', 'data', 'msa_configuration.fits'])
    hdul[0].header['MSAMETFL'] = msa_meta
    hdul[0].header['MSAMETID'] = 12
    im = datamodels.ImageModel(hdul)
    im_wcs = AssignWcsStep.call(im)
    im_ex2d = Extract2dStep.call(im_wcs)
    im_ex2d.slits[0].meta.wcs.bounding_box = ((-.5, 1432.5), (-.5, 37.5))
    im_src = SourceTypeStep.call(im_ex2d)
    im_wave = WavecorrStep.call(im_src)

    # test dispersion is of the correct order
    # there's one slit only
    slit = im_src.slits[0]
    x, y = wcstools.grid_from_bounding_box(slit.meta.wcs.bounding_box)
    dispersion = wavecorr.compute_dispersion(slit.meta.wcs, x, y)
    assert_allclose(dispersion[~np.isnan(dispersion)], 1e-9, atol=1e-10)

    # the difference in wavelength should be of the order of e-10 in um
    assert_allclose(im_src.slits[0].wavelength - im_wave.slits[0].wavelength, 1e-10)

    # test on both sides of the shutter
    source_xpos1 = -.2
    source_xpos2 = .2

    ra, dec, lam = slit.meta.wcs(x, y)
    ref_name = im_wave.meta.ref_file.wavecorr.name
    freference = datamodels.WaveCorrModel(WavecorrStep.reference_uri_to_cache_path(ref_name, im.crds_observatory))
    zero_point1 = wavecorr.compute_zero_point_correction(lam, freference, source_xpos1, 'MOS', dispersion)
    zero_point2 = wavecorr.compute_zero_point_correction(lam, freference, source_xpos2, 'MOS', dispersion)
    diff_correction = np.abs(zero_point1[1] - zero_point2[1])
    non_zero = diff_correction[diff_correction != 0]
    assert_allclose(np.nanmean(non_zero), 0.75, atol=0.01)


def test_ideal_to_v23_fs():
    """Test roundtripping between Ideal and V2V3 frames."""
    v3yangle = 138.78
    v2_ref = 321.87
    v3_ref = -477.94
    vparity = -1
    id2v = models.IdealToV2V3(v3yangle, v2_ref, v3_ref, vparity)
    assert_allclose(id2v(0, 0), (v2_ref, v3_ref))
    assert_allclose(id2v.inverse(v2_ref, v3_ref), (0, 0))


def test_skipped():
    """ Test all conditions that lead to skipping wavecorr."""

    hdul = create_nirspec_fs_file(grating="G140H", filter="F100LP")
    im = datamodels.ImageModel(hdul)

    # test a non-valid exp_type2transform
    im.meta.exposure.type = 'NRS_IMAGE'
    out = WavecorrStep.call(im)
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    # Test an error is raised if assign_wcs or extract_2d were not run.
    im.meta.exposure.type = 'NRS_FIXEDSLIT'
    with pytest.raises(AttributeError):
        WavecorrStep.call(im)

    outa = AssignWcsStep.call(im)
    oute = Extract2dStep.call(outa)
    outs = SourceTypeStep.call(oute)

    # Test step is skipped if no coverage in CRDS
    outs.slits[0].meta.observation.date = '2001-08-03'
    outw = WavecorrStep.call(outs)
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    outs.meta.observation.date = '2017-08-03'
    outw = WavecorrStep.call(outs)
    # Primary name not set
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    outs.meta.instrument.fixed_slit = "S400A1"

    # Test step is skipped if meta.dither is not populated
    outw = WavecorrStep.call(outs)
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    dither = {'x_offset': 0.0, 'y_offset': 0.0}
    ind = np.nonzero([s.name == 'S400A1' for s in outs.slits])[0].item()
    outs.slits[ind].meta.dither = dither
    outs.slits[ind].meta.wcsinfo.v3yangle = 138.78
    outs.slits[ind].meta.wcsinfo.vparity = -1
    outs.slits[ind].meta.wcsinfo.v2_ref = 321.87
    outs.slits[ind].meta.wcsinfo.v3_ref = -477.94
    outs.slits[ind].meta.wcsinfo.roll_ref = 15.1234

    # Test step is skipped if source is "EXTENDED"
    outw = WavecorrStep.call(outs)
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    outs.slits[ind].source_type = 'POINT'
    outw = WavecorrStep.call(outs)

    source_pos = (0.004938526981283373, -0.02795306204991911)
    assert_allclose((outw.slits[ind].source_xpos, outw.slits[ind].source_ypos), source_pos)


def test_wavecorr_fs():
    hdul = create_nirspec_fs_file(grating="PRISM", filter="CLEAR")
    im = datamodels.ImageModel(hdul)
    dither = {'x_offset': -0.0264, 'y_offset': 1.089798712}

    im.meta.dither = dither
    im.meta.instrument.fixed_slit = 'S200A1'
    im.meta.instrument.lamp_state = 'NONE'
    im.meta.instrument.lamp_mode = 'FIXEDSLIT'
    im.meta.instrument.gwa_tilt = 698.1256999999999
    im.meta.instrument.gwa_xtilt = 0.334691644
    im.meta.instrument.gwa_ytilt = 0.0349255353
    im.meta.exposure.type = 'NRS_FIXEDSLIT'
    im.meta.wcsinfo = {'dec_ref': -70.77497509320972,
                       'dispersion_direction': 1,
                       'ra_ref': 90.75414044948525,
                       'roll_ref': 38.2179067606001,
                       'specsys': 'BARYCENT',
                       'spectral_order': 0,
                       'v2_ref': 332.136078,
                       'v3_ref': -479.224213,
                       'v3yangle': 138.7605896,
                       'velosys': -1398.2,
                       'vparity': -1,
                       'waverange_end': 5.3e-06,
                       'waverange_start': 6e-07}

    result = AssignWcsStep.call(im)
    result = Extract2dStep.call(result)
    bbox = ((-.5, 428.5), (-.5, 38.5))
    result.slits[0].meta.wcs.bounding_box = bbox
    x, y = wcstools.grid_from_bounding_box(bbox)
    ra, dec, lam_before = result.slits[0].meta.wcs(x, y)
    result.slits[0].wavelength = lam_before
    src_result = SourceTypeStep.call(result)
    result = WavecorrStep.call(src_result)

    assert_allclose(result.slits[0].source_xpos, 0.127111, atol=1e-6)

    slit = result.slits[0]
    source_xpos = wavecorr.get_source_xpos(slit, slit.meta.wcs, lam=2)
    assert_allclose(result.slits[0].source_xpos, source_xpos, atol=1e-6)

    mean_correction = np.abs(src_result.slits[0].wavelength - result.slits[0].wavelength)
    assert_allclose(np.nanmean(mean_correction), 0.003, atol=.001)

    dispersion = wavecorr.compute_dispersion(slit.meta.wcs, x, y)
    assert_allclose(dispersion[~np.isnan(dispersion)], 1e-8, atol=1.04e-8)

    # test on both sides of the slit center
    source_xpos1 = -.2
    source_xpos2 = .2

    ref_name = result.meta.ref_file.wavecorr.name
    freference = datamodels.WaveCorrModel(WavecorrStep.reference_uri_to_cache_path(ref_name, im.crds_observatory))
    zero_point1 = wavecorr.compute_zero_point_correction(lam_before, freference, source_xpos1, 'S200A1', dispersion)
    zero_point2 = wavecorr.compute_zero_point_correction(lam_before, freference, source_xpos2, 'S200A1', dispersion)
    diff_correction = np.abs(zero_point1[1] - zero_point2[1])
    assert_allclose(np.nanmean(diff_correction), 0.45, atol=0.01)
