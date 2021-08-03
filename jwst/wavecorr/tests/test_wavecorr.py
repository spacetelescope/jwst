import os.path
import numpy as np
from numpy.testing import assert_allclose
import pytest

from gwcs import wcstools
import jwst
from jwst import datamodels
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d import Extract2dStep
from jwst.srctype import SourceTypeStep
from jwst.transforms import models
from jwst.wavecorr import WavecorrStep
from jwst.wavecorr import wavecorr
from crds import CrdsLookupError

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
    assert_allclose(diff_correction[diff_correction.nonzero()].mean(), 0.02, atol=0.01)


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
    im.meta.exposure.type='NRS_IMAGE'
    out = WavecorrStep.call(im)
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    # Test an error is raised if assign_wcs or extract_2d were not run.
    im.meta.exposure.type='NRS_FIXEDSLIT'
    with pytest.raises(AttributeError):
        WavecorrStep.call(im)

    outa = AssignWcsStep.call(im)
    oute = Extract2dStep.call(outa)
    outs = SourceTypeStep.call(oute)

    # Test step is skipped if no coverage in CRDS
    outs.slits[0].meta.observation.date='2001-08-03'
    outw = WavecorrStep.call(outs)
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    outs.meta.observation.date='2017-08-03'
    outw = WavecorrStep.call(outs)
    # Primary name not set
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    outs.meta.instrument.fixed_slit = "S400A1"

    # Test step is skipped if meta.dither is not populated
    outw = WavecorrStep.call(outs)
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    dither = {'x_offset': 0.0, 'y_offset': 0.0}
    ind = np.nonzero([s.name=='S400A1' for s in outs.slits])[0].item()
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
