import os
import pytest

from astropy.io.fits.diff import FITSDiff
import numpy as np

from jwst.barshadow import BarShadowStep
import jwst.datamodels as dm
from jwst.flatfield import FlatFieldStep
from jwst.pathloss import PathLossStep
from jwst.photom import PhotomStep
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run the calwebb_spec2 pipeline on a single NIRSpec MOS exposure."""

    rtdata = rtdata_module

    # Get the cfg files
    collect_pipeline_cfgs("config")

    # Get the MSA metadata file referenced in the input exposure
    rtdata.get_data("nirspec/mos/jw95065006001_0_short_msa.fits")

    # Get the input exposure
    rtdata.get_data("nirspec/mos/f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod.fits")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = ["config/calwebb_spec2.cfg", rtdata.input,
            "--steps.assign_wcs.save_results=true",
            "--steps.msa_flagging.save_results=true",
            "--steps.extract_2d.save_results=true",
            "--steps.srctype.save_results=true",
            "--steps.wavecorr.save_results=true",
            "--steps.flat_field.save_results=true",
            "--steps.pathloss.save_results=true",
            "--steps.barshadow.save_results=true"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("output",[
    "assign_wcs", "msa_flagging", "extract_2d", "wavecorr", "flat_field", "srctype",
    "pathloss", "barshadow", "cal", "s2d", "x1d"])
def test_nirspec_mos_spec2(run_pipeline, fitsdiff_default_kwargs, output):
    """Regression test of the calwebb_spec2 pipeline on a
       NIRSpec MOS exposure."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    rtdata.output = "f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod_" + output + ".fits"

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_mos_spec2",
                                  "f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod_" + output + ".fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_flat_field_step_user_supplied_flat(jail, rtdata_module, fitsdiff_default_kwargs):
    """Test providing a user-supplied flat field to the FlatFieldStep"""
    rtdata = rtdata_module
    data = rtdata.get_data('nirspec/mos/usf_wavecorr.fits')
    user_supplied_flat = rtdata.get_data('nirspec/mos/usf_flat.fits')

    data_flat_fielded = FlatFieldStep.call(data, user_supplied_flat=user_supplied_flat)
    rtdata.output = 'flat_fielded_step_user_supplied.fits'
    data_flat_fielded.write(rtdata.output)

    rtdata.get_truth('truth/test_nirspec_mos_spec2/flat_fielded_step_user_supplied.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_ff_inv(jail, rtdata_module, fitsdiff_default_kwargs):
    """Test flat field inversion"""
    rtdata = rtdata_module
    data = dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits'))

    flatted = FlatFieldStep.call(data)
    unflatted = FlatFieldStep.call(flatted, inverse=True)

    bad_slits = []
    for idx, slits in enumerate(zip(data.slits, unflatted.slits)):
        data_slit, unflatted_slit = slits
        if not np.allclose(data_slit.data, unflatted_slit.data):
            bad_slits.append(idx)
    assert not bad_slits, f'Inversion failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_pathloss_corrpars(jail, rtdata_module):
    """Test PathLossStep using correction_pars"""
    rtdata = rtdata_module
    data = dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits'))

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
def test_pathloss_inverse(jail, rtdata_module):
    """Test PathLossStep using inversion"""
    rtdata = rtdata_module
    data = dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits'))

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
def test_pathloss_source_type(jail, rtdata_module):
    """Test PathLossStep forcing source type"""
    rtdata = rtdata_module
    data = dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits'))

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
def test_barshadow_corrpars(jail, rtdata_module):
    """BarShadowStep using correction_pars"""
    rtdata = rtdata_module
    data = dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits'))

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
def test_barshadow_inverse(jail, rtdata_module):
    """BarShadowStep using inversion"""
    rtdata = rtdata_module
    data = dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits'))

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
def test_barshadow_source_type(jail, rtdata_module):
    """Test BarShadowStep forcing source type"""
    rtdata = rtdata_module
    data = dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits'))

    pls = BarShadowStep()
    pls.source_type = 'extended'
    corrected = pls.run(data)

    bad_slits = []
    for idx, slit in enumerate(corrected.slits):
        if np.allclose(slit.barshadow, np.ones(slit.data.shape), equal_nan=True):
            bad_slits.append(idx)
    assert not bad_slits, f'Force to uniform failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_photom_corrpars(jail, rtdata_module):
    """Test for photom correction parameters"""
    rtdata = rtdata_module
    data = dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits'))

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
def test_photom_inverse(jail, rtdata_module):
    """PhotomStep using inversion"""
    rtdata = rtdata_module
    data = dm.open(rtdata.get_data('nirspec/mos/usf_wavecorr.fits'))

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
