import numpy as np
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.barshadow import BarShadowStep
import jwst.datamodels as dm
from jwst.flatfield import FlatFieldStep
from jwst.pathloss import PathLossStep
from jwst.photom import PhotomStep


@pytest.mark.bigdata
def test_flat_field_step_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test providing a user-supplied flat field to the FlatFieldStep"""
    data = rtdata.get_data('nirspec/mos/usf_wavecorr.fits')
    user_supplied_flat = rtdata.get_data('nirspec/mos/usf_flat.fits')

    data_flat_fielded = FlatFieldStep.call(data, user_supplied_flat=user_supplied_flat)
    rtdata.output = 'flat_fielded_step_user_supplied.fits'
    data_flat_fielded.write(rtdata.output)

    rtdata.get_truth('truth/test_nirspec_mos_spec2/flat_fielded_step_user_supplied.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_ff_inv(rtdata, fitsdiff_default_kwargs):
    """Test flat field inversion"""
    with dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits')) as data:
        flatted = FlatFieldStep.call(data)
        unflatted = FlatFieldStep.call(flatted, inverse=True)

        bad_slits = []
        for idx, slits in enumerate(zip(data.slits, unflatted.slits)):
            data_slit, unflatted_slit = slits
            if not np.allclose(data_slit.data, unflatted_slit.data):
                bad_slits.append(idx)

    assert not bad_slits, f'Inversion failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_pathloss_corrpars(rtdata):
    """Test PathLossStep using correction_pars"""
    with dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits')) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.use_correction_pars = True
        corrected_corrpars = pls.run(data)

    bad_slits = []
    for idx, slits in enumerate(zip(corrected.slits, corrected_corrpars.slits)):
        corrected_slit, corrected_corrpars_slit = slits
        if not np.allclose(corrected_slit.data, corrected_corrpars_slit.data, equal_nan=True):
            bad_slits.append(idx)
    assert not bad_slits, f'correction_pars failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_pathloss_inverse(rtdata):
    """Test PathLossStep using inversion"""
    with dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits')) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.inverse = True
        corrected_inverse = pls.run(corrected)

        bad_slits = []
        for idx, slits in enumerate(zip(data.slits, corrected_inverse.slits)):
            data_slit, corrected_inverse_slit = slits
            non_nan = ~np.isnan(corrected_inverse_slit.data)
            if not np.allclose(data_slit.data[non_nan], corrected_inverse_slit.data[non_nan]):
                bad_slits.append(idx)

    assert not bad_slits, f'Inversion failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_pathloss_source_type(rtdata):
    """Test PathLossStep forcing source type"""
    with dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits')) as data:
        pls = PathLossStep()
        pls.source_type = 'extended'
        pls.run(data)

    bad_slits = []
    for idx, slit in enumerate(pls.correction_pars.slits):
        if slit:
            if not np.allclose(slit.data, slit.pathloss_uniform, equal_nan=True):
                bad_slits.append(idx)
    assert not bad_slits, f'Force to uniform failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_barshadow_corrpars(rtdata):
    """BarShadowStep using correction_pars"""
    with dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits')) as data:
        pls = BarShadowStep()
        corrected = pls.run(data)

        pls.use_correction_pars = True
        corrected_corrpars = pls.run(data)

    bad_slits = []
    for idx, slits in enumerate(zip(corrected.slits, corrected_corrpars.slits)):
        corrected_slit, corrected_corrpars_slit = slits
        if not np.allclose(corrected_slit.data, corrected_corrpars_slit.data, equal_nan=True):
            bad_slits.append(idx)
    assert not bad_slits, f'correction_pars failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_barshadow_inverse(rtdata):
    """BarShadowStep using inversion"""
    with dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits')) as data:
        pls = BarShadowStep()
        corrected = pls.run(data)

        pls.inverse = True
        corrected_inverse = pls.run(corrected)

        bad_slits = []
        for idx, slits in enumerate(zip(data.slits, corrected_inverse.slits)):
            data_slit, corrected_inverse_slit = slits
            non_nan = ~np.isnan(corrected_inverse_slit.data)
            if not np.allclose(data_slit.data[non_nan], corrected_inverse_slit.data[non_nan]):
                bad_slits.append(idx)

    assert not bad_slits, f'Inversion failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_barshadow_source_type(rtdata):
    """Test BarShadowStep forcing source type"""
    with dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits')) as data:
        pls = BarShadowStep()
        pls.source_type = 'extended'
        corrected = pls.run(data)

    bad_slits = []
    for idx, slit in enumerate(corrected.slits):
        if np.allclose(slit.barshadow, np.ones(slit.data.shape), equal_nan=True):
            bad_slits.append(idx)
    assert not bad_slits, f'Force to uniform failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_photom_corrpars(rtdata):
    """Test for photom correction parameters"""
    with dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits')) as data:
        pls = PhotomStep()
        corrected = pls.run(data)

        pls.use_correction_pars = True
        corrected_corrpars = pls.run(data)

    bad_slits = []
    for idx, slits in enumerate(zip(corrected.slits, corrected_corrpars.slits)):
        corrected_slit, corrected_corrpars_slit = slits
        if not np.allclose(corrected_slit.data, corrected_corrpars_slit.data, equal_nan=True):
            bad_slits.append(idx)
    assert not bad_slits, f'correction_pars failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_photom_inverse(rtdata):
    """PhotomStep using inversion"""
    with dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits')) as data:
        pls = PhotomStep()
        corrected = pls.run(data)

        pls.inverse = True
        corrected_inverse = pls.run(corrected)

        bad_slits = []
        for idx, slits in enumerate(zip(data.slits, corrected_inverse.slits)):
            data_slit, corrected_inverse_slit = slits
            non_nan = ~np.isnan(corrected_inverse_slit.data)
            if not np.allclose(data_slit.data[non_nan], corrected_inverse_slit.data[non_nan]):
                bad_slits.append(idx)

    assert not bad_slits, f'Inversion failed for slits {bad_slits}'
