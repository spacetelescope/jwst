from pathlib import Path

import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_filename
from gwcs import wcstools
from numpy.testing import assert_allclose
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.transforms import models
from stpipe.crds_client import reference_uri_to_cache_path

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_mos_file, create_nirspec_fs_file
from jwst.extract_2d import Extract2dStep
from jwst.srctype import SourceTypeStep
from jwst.wavecorr import WavecorrStep
from jwst.wavecorr import wavecorr


@pytest.fixture(scope="module")
def nrs_fs_model():
    hdul = create_nirspec_fs_file(grating="G140H", filter="F100LP")
    im = datamodels.ImageModel(hdul)
    im_wcs = AssignWcsStep.call(im)
    im_ex2d = Extract2dStep.call(im_wcs)
    yield im_ex2d
    im_ex2d.close()
    im_wcs.close()
    im.close()
    hdul.close()


@pytest.fixture(scope="module")
def nrs_slit_model(nrs_fs_model):
    im_ex2d = nrs_fs_model.copy()

    # make a slit model to run through correction
    slit = datamodels.SlitModel(im_ex2d.slits[0].data)
    slit.update(im_ex2d)
    slit.meta.wcs = im_ex2d.slits[0].meta.wcs
    slit.source_type = "POINT"
    slit.name = "S1600A1"
    yield slit
    slit.close()


def test_wavecorr():
    hdul = create_nirspec_mos_file()
    msa_meta = get_pkg_data_filename("data/msa_configuration.fits", package="jwst.assign_wcs.tests")
    hdul[0].header["MSAMETFL"] = msa_meta
    hdul[0].header["MSAMETID"] = 12
    im = datamodels.ImageModel(hdul)
    im_wcs = AssignWcsStep.call(im)
    im_ex2d = Extract2dStep.call(im_wcs)
    bbox = ((-0.5, 1432.5), (-0.5, 37.5))
    im_ex2d.slits[0].meta.wcs.bounding_box = bbox
    x, y = wcstools.grid_from_bounding_box(bbox)
    ra, dec, lam_before = im_ex2d.slits[0].meta.wcs(x, y)
    im_ex2d.slits[0].wavelength = lam_before
    im_src = SourceTypeStep.call(im_ex2d)

    # the mock msa source is an extended source, change to point for testing
    im_src.slits[0].source_type = "POINT"
    im_wave = WavecorrStep.call(im_src)

    # test dispersion is of the correct order
    # there's one slit only
    slit = im_src.slits[0]
    dispersion = wavecorr.compute_dispersion(slit.meta.wcs, x, y)
    assert_allclose(dispersion[~np.isnan(dispersion)], 1e-9, atol=1e-10)

    # test that the wavelength is on the order of microns
    wavelength = wavecorr.compute_wavelength(slit.meta.wcs, x, y)
    assert_allclose(np.nanmean(wavelength), 2.5, atol=0.1)

    # Check that the stored wavelengths were corrected
    abs_wave_correction = np.abs(im_src.slits[0].wavelength - im_wave.slits[0].wavelength)
    assert_allclose(np.nanmean(abs_wave_correction), 0.00046, atol=0.0001)

    # Check that the slit wcs has been updated to provide corrected wavelengths
    corrected_wavelength = wavecorr.compute_wavelength(im_wave.slits[0].meta.wcs, x, y)
    assert_allclose(im_wave.slits[0].wavelength, corrected_wavelength)

    # test the round-tripping on the wavelength correction transform
    ref_name = im_wave.meta.ref_file.wavecorr.name
    freference = datamodels.WaveCorrModel(
        reference_uri_to_cache_path(ref_name, im.crds_observatory)
    )

    lam_uncorr = lam_before * 1e-6
    wave2wavecorr = wavecorr.calculate_wavelength_correction_transform(
        lam_uncorr, dispersion, freference, slit.source_xpos, "MOS"
    )
    lam_corr = wave2wavecorr(lam_uncorr)
    assert_allclose(lam_uncorr, wave2wavecorr.inverse(lam_corr))

    # test on both sides of the shutter
    source_xpos1 = -0.2
    source_xpos2 = 0.2

    wave_transform1 = wavecorr.calculate_wavelength_correction_transform(
        lam_uncorr, dispersion, freference, source_xpos1, "MOS"
    )
    wave_transform2 = wavecorr.calculate_wavelength_correction_transform(
        lam_uncorr, dispersion, freference, source_xpos2, "MOS"
    )

    zero_point1 = wave_transform1(lam_uncorr)
    zero_point2 = wave_transform2(lam_uncorr)

    diff_correction = np.abs(zero_point1 - zero_point2)
    assert_allclose(np.nanmean(diff_correction), 8.0e-10, atol=0.1e-10)


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
    """Test all conditions that lead to skipping wavecorr."""
    hdul = create_nirspec_fs_file(grating="G140H", filter="F100LP")
    im = datamodels.ImageModel(hdul)

    # test a non-valid exp_type2transform
    im.meta.exposure.type = "NRS_IMAGE"
    out = WavecorrStep.call(im)
    assert out.meta.cal_step.wavecorr == "SKIPPED"

    # Test an error is raised if assign_wcs or extract_2d were not run.
    im.meta.exposure.type = "NRS_FIXEDSLIT"
    with pytest.raises(TypeError):
        WavecorrStep.call(im)

    outa = AssignWcsStep.call(im)

    # wcs and other metadata necessary for a good (not skipped) result
    dither = {"x_offset": 0.0, "y_offset": 0.0}
    outa.meta.dither = dither
    outa.meta.wcsinfo.v3yangle = 138.78
    outa.meta.wcsinfo.vparity = -1
    outa.meta.wcsinfo.v2_ref = 321.87
    outa.meta.wcsinfo.v3_ref = -477.94
    outa.meta.wcsinfo.roll_ref = 15.1234
    outa.meta.instrument.fixed_slit = "S400A1"

    oute = Extract2dStep.call(outa)

    # ensure the source position is set as expected by extract2d
    ind = np.nonzero([s.name == "S400A1" for s in oute.slits])[0].item()
    source_pos = (0.004938526981283373, -0.02795306204991911)
    assert_allclose((oute.slits[ind].source_xpos, oute.slits[ind].source_ypos), source_pos)

    outs = SourceTypeStep.call(oute)

    # primary slit is found by source type step to be EXTENDED
    # which means its source position is reset to (0, 0)
    # hack it here to be a point source with the same offset as before here
    outs.slits[ind].source_type = "POINT"
    outs.slits[ind].source_xpos = source_pos[0]
    outs.slits[ind].source_ypos = source_pos[1]

    # Test step is skipped if no coverage in CRDS
    outs.meta.observation.date = "2001-08-03"
    outw = WavecorrStep.call(outs)
    assert outw.meta.cal_step.wavecorr == "SKIPPED"
    outs.meta.observation.date = "2017-08-03"

    # Run the step for real
    outw = WavecorrStep.call(outs)

    # Test if the corrected wavelengths are not monotonically increasing

    # This case is not expected with real data, test that no correction
    # transform is returned (a skip criterion) with simple case
    # of flipped wavelength solutions, which produces a monotonically
    # decreasing wavelengths
    lam = np.tile(np.flip(np.arange(0.6, 5.5, 0.01) * 1e-6), (22, 1))
    disp = np.tile(np.full(lam.shape[-1], -0.01) * 1e-6, (22, 1))

    ref_name = outw.meta.ref_file.wavecorr.name
    reffile = datamodels.WaveCorrModel(reference_uri_to_cache_path(ref_name, im.crds_observatory))
    source_xpos = 0.1
    aperture_name = "S200A1"

    transform = wavecorr.calculate_wavelength_correction_transform(
        lam, disp, reffile, source_xpos, aperture_name
    )
    assert transform is None


def test_mos_slit_status():
    """Test conditions that are skipped for mos slitlets."""
    hdul = create_nirspec_mos_file()
    msa_meta = get_pkg_data_filename("data/msa_configuration.fits", package="jwst.assign_wcs.tests")
    hdul[0].header["MSAMETFL"] = msa_meta
    hdul[0].header["MSAMETID"] = 12
    im = datamodels.ImageModel(hdul)
    im_wcs = AssignWcsStep.call(im)
    im_ex2d = Extract2dStep.call(im_wcs)
    bbox = ((-0.5, 1432.5), (-0.5, 37.5))
    im_ex2d.slits[0].meta.wcs.bounding_box = bbox
    x, y = wcstools.grid_from_bounding_box(bbox)
    ra, dec, lam_before = im_ex2d.slits[0].meta.wcs(x, y)
    im_ex2d.slits[0].wavelength = lam_before
    im_src = SourceTypeStep.call(im_ex2d)

    # test the mock msa source as an extended source
    im_src.slits[0].source_type = "EXTENDED"
    im_wave = WavecorrStep.call(im_src)

    # check that the step is recorded as skipped,
    # since no slits were corrected
    assert im_wave.meta.cal_step.wavecorr == "SKIPPED"

    # check that wavelength_corrected is False for extended mos sources
    assert im_wave.slits[0].wavelength_corrected is False

    # test the mock msa source as a point source
    im_src.slits[0].source_type = "POINT"
    im_wave = WavecorrStep.call(im_src)

    # check that the step is recorded as completed
    assert im_wave.meta.cal_step.wavecorr == "COMPLETE"

    # check that wavelength_corrected is True for mos point sources
    assert im_wave.slits[0].wavelength_corrected is True


def test_wavecorr_fs():
    hdul = create_nirspec_fs_file(grating="PRISM", filter="CLEAR")
    im = datamodels.ImageModel(hdul)
    dither = {"x_offset": -0.0264, "y_offset": 1.089798712}

    im.meta.dither = dither
    im.meta.instrument.fixed_slit = "S200A1"
    im.meta.instrument.lamp_state = "NONE"
    im.meta.instrument.lamp_mode = "FIXEDSLIT"
    im.meta.instrument.gwa_tilt = 698.1256999999999
    im.meta.instrument.gwa_xtilt = 0.334691644
    im.meta.instrument.gwa_ytilt = 0.0349255353
    im.meta.exposure.type = "NRS_FIXEDSLIT"
    im.meta.wcsinfo = {
        "dec_ref": -70.77497509320972,
        "dispersion_direction": 1,
        "ra_ref": 90.75414044948525,
        "roll_ref": 38.2179067606001,
        "specsys": "BARYCENT",
        "spectral_order": 0,
        "v2_ref": 332.136078,
        "v3_ref": -479.224213,
        "v3yangle": 138.7605896,
        "velosys": -1398.2,
        "vparity": -1,
        "waverange_end": 5.3e-06,
        "waverange_start": 6e-07,
    }

    result = AssignWcsStep.call(im)
    result = Extract2dStep.call(result)
    bbox = ((-0.5, 428.5), (-0.5, 38.5))
    result.slits[0].meta.wcs.bounding_box = bbox
    x, y = wcstools.grid_from_bounding_box(bbox)
    ra, dec, lam_before = result.slits[0].meta.wcs(x, y)
    result.slits[0].wavelength = lam_before
    src_result = SourceTypeStep.call(result)
    result = WavecorrStep.call(src_result)

    assert_allclose(result.slits[0].source_xpos, 0.127111, atol=1e-6)

    slit = result.slits[0]

    mean_correction = np.abs(src_result.slits[0].wavelength - result.slits[0].wavelength)
    assert_allclose(np.nanmean(mean_correction), 0.003, atol=0.001)

    dispersion = wavecorr.compute_dispersion(slit.meta.wcs, x, y)
    assert_allclose(dispersion[~np.isnan(dispersion)], 1e-8, atol=1.04e-8)

    # Check that the slit wavelengths are consistent with the slit wcs
    corrected_wavelength = wavecorr.compute_wavelength(slit.meta.wcs, x, y)
    assert_allclose(slit.wavelength, corrected_wavelength)

    # test the round-tripping on the wavelength correction transform
    ref_name = result.meta.ref_file.wavecorr.name
    freference = datamodels.WaveCorrModel(
        reference_uri_to_cache_path(ref_name, im.crds_observatory)
    )

    lam_uncorr = lam_before * 1e-6
    wave2wavecorr = wavecorr.calculate_wavelength_correction_transform(
        lam_uncorr, dispersion, freference, slit.source_xpos, "S200A1"
    )
    lam_corr = wave2wavecorr(lam_uncorr)
    assert_allclose(lam_uncorr, wave2wavecorr.inverse(lam_corr))

    # test on both sides of the slit center
    source_xpos1 = -0.2
    source_xpos2 = 0.2

    wave_transform1 = wavecorr.calculate_wavelength_correction_transform(
        lam_uncorr, dispersion, freference, source_xpos1, "S200A1"
    )
    wave_transform2 = wavecorr.calculate_wavelength_correction_transform(
        lam_uncorr, dispersion, freference, source_xpos2, "S200A1"
    )

    zero_point1 = wave_transform1(lam_uncorr)
    zero_point2 = wave_transform2(lam_uncorr)

    diff_correction = np.abs(zero_point1 - zero_point2)
    assert_allclose(np.nanmean(diff_correction), 6.3e-9, atol=0.1e-9)


def test_assign_wcs_skipped():
    hdul = create_nirspec_fs_file(grating="G140H", filter="F100LP")
    im = datamodels.ImageModel(hdul)

    # if assign_wcs was skipped, wavecorr is also skipped
    im.meta.cal_step.assign_wcs = "SKIPPED"
    result = WavecorrStep.call(im)
    assert result.meta.cal_step.wavecorr == "SKIPPED"

    hdul.close()
    im.close()
    result.close()


def test_missing_wcs():
    hdul = create_nirspec_fs_file(grating="G140H", filter="F100LP")
    im = datamodels.SlitModel(hdul)

    with pytest.raises(AttributeError, match="does not have a WCS"):
        WavecorrStep.call(im)

    hdul.close()
    im.close()


def test_invalid_exptype(nrs_fs_model):
    im_ex2d = nrs_fs_model.copy()

    # Test the skip at the do_correction level
    im_ex2d.meta.exposure.type = "ANY"
    result = wavecorr.do_correction(im_ex2d, None)
    assert result.meta.cal_step.wavecorr == "SKIPPED"


def test_invalid_slit(nrs_slit_model):
    slit = nrs_slit_model.copy()
    slit.name = None
    result = WavecorrStep.call(slit)
    assert result.meta.cal_step.wavecorr == "SKIPPED"


@pytest.mark.parametrize("source_type", ["POINT", "EXTENDED"])
def test_slitmodel(source_type, nrs_slit_model):
    slit = nrs_slit_model.copy()
    slit.source_type = source_type

    result = WavecorrStep.call(slit)
    if source_type == "POINT":
        assert result.meta.cal_step.wavecorr == "COMPLETE"
    else:
        assert result.meta.cal_step.wavecorr == "SKIPPED"

    slit.close()
    result.close()


@pytest.mark.parametrize("level", ["top", "meta", None])
@pytest.mark.parametrize("value", ["POINT", "EXTENDED", None])
def test_is_point_source(level, value):
    model = datamodels.SlitModel()
    if level == "top":
        model.source_type = value
        model.meta.target.source_type = None
    elif level == "meta":
        model.source_type = None
        model.meta.target.source_type = value
    else:
        model.source_type = None
        model.meta.target.source_type = None
    if level is not None and value == "POINT":
        assert wavecorr._is_point_source(model) is True
    else:
        assert wavecorr._is_point_source(model) is False
