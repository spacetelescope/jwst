import os.path
import numpy as np
from numpy.testing import assert_allclose

from gwcs import wcstools
import jwst
from jwst import datamodels
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d import Extract2dStep
from jwst.srctype import SourceTypeStep
from jwst.wavecorr import WavecorrStep
from jwst.wavecorr import wavecorr

from jwst.assign_wcs.tests.test_nirspec import create_nirspec_mos_file


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
    freference = datamodels.WaveCorrModel(WavecorrStep.reference_uri_to_cache_path(ref_name))
    zero_point1 = wavecorr.compute_zero_point_correction(lam, freference, source_xpos1, 'MOS', dispersion)
    zero_point2 = wavecorr.compute_zero_point_correction(lam, freference, source_xpos2, 'MOS', dispersion)
    diff_correction = np.abs(zero_point1[1] - zero_point2[1])
    assert_allclose(diff_correction[diff_correction.nonzero()].mean(), 0.02, atol=0.01)
