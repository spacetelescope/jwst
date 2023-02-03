import pytest

from astropy.io.fits.diff import FITSDiff
import numpy as np

import stdatamodels.jwst.datamodels as dm

from jwst.master_background import MasterBackgroundMosStep, MasterBackgroundStep

from jwst.regtest import regtestdata as rt

pytestmark = pytest.mark.bigdata


@pytest.fixture(scope='module')
def run_spec2_mbkg(jail, rtdata_module):
    """Run Spec2 on MSA data with background slits"""
    rtdata = rtdata_module

    # Get data
    rtdata.get_data('nirspec/mos/jw01180025001_01_msa.fits')
    rtdata.get_data('nirspec/mos/jw01180025001_05101_00001_nrs2_rate.fits')

    # Run the pipeline
    step_params = {
        'step': 'calwebb_spec2',
        'args': [
            '--steps.master_background_mos.skip=false',
            '--steps.master_background_mos.save_background=true'
        ]
    }
    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


@pytest.fixture(scope='module')
def run_spec2_mbkg_user(jail, rtdata_module):
    """Run Spec2 on MSA data with a user-supplied master bkg spectrum"""
    rtdata = rtdata_module

    # Get data
    rtdata.get_data('nirspec/mos/nrs_mos_3pointnod_1_msa.fits')
    rtdata.get_data('nirspec/mos/jw00626030001_02103_00001_nrs1_masterbg1d.fits')
    rtdata.get_data('nirspec/mos/jw00626030001_02103_00001_nrs1_rate.fits')

    # Run the pipeline
    step_params = {
        'step': 'calwebb_spec2',
        'args': [
            '--steps.master_background_mos.skip=false',
            '--steps.master_background_mos.user_background=jw00626030001_02103_00001_nrs1_masterbg1d.fits'
        ]
    }
    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


def test_masterbkg_rerun(rtdata):
    """Test to ensure sequential runs of the step are consistent"""
    with dm.open(rtdata.get_data('nirspec/mos/nrs_mos_with_bkgslits_srctype.fits')) as data:
        mbs = MasterBackgroundMosStep()
        corrected = mbs.run(data)
        corrected_again = mbs.run(data)

    bad_slits = []
    for idx, slits in enumerate(zip(corrected.slits, corrected_again.slits)):
        corrected_slit, corrected_again_slit = slits
        if not np.allclose(corrected_slit.data, corrected_again_slit.data, equal_nan=True):
            bad_slits.append(idx)
    assert not bad_slits, f'rerun failed for slits {bad_slits}'


def test_masterbkg_corrpars(rtdata):
    """Test for correction parameters"""
    with dm.open(rtdata.get_data('nirspec/mos/nrs_mos_with_bkgslits_srctype.fits')) as data:
        mbs = MasterBackgroundMosStep()
        corrected = mbs.run(data)

        mbs.use_correction_pars = True
        corrected_corrpars = mbs.run(data)

    bad_slits = []
    for idx, slits in enumerate(zip(corrected.slits, corrected_corrpars.slits)):
        corrected_slit, corrected_corrpars_slit = slits
        if not np.allclose(corrected_slit.data, corrected_corrpars_slit.data, equal_nan=True):
            bad_slits.append(idx)
    assert not bad_slits, f'correction_pars failed for slits {bad_slits}'


@pytest.mark.parametrize(
    'suffix',
    ['masterbg1d', 'masterbg2d', 'cal', 's2d', 'x1d']
)
def test_nirspec_mos_mbkg(suffix, run_spec2_mbkg, fitsdiff_default_kwargs):
    """Run spec2 with master background"""
    rtdata = run_spec2_mbkg
    rt.is_like_truth(rtdata, fitsdiff_default_kwargs, suffix, truth_path='truth/test_nirspec_mos_mbkg')


@pytest.mark.parametrize(
    'suffix',
    ['cal', 's2d', 'x1d']
)
def test_nirspec_mos_mbkg_user(suffix, run_spec2_mbkg_user, fitsdiff_default_kwargs):
    """Run spec2 with master background and user-supplied mbkg"""
    rtdata = run_spec2_mbkg_user
    rt.is_like_truth(rtdata, fitsdiff_default_kwargs, suffix, truth_path='truth/test_nirspec_mos_mbkg_user')


def test_nirspec_fs_mbkg_user(rtdata, fitsdiff_default_kwargs):
    """Run a test for NIRSpec FS data with a user-supplied background file."""

    # Get user-supplied background
    user_background = "v2_nrs_bkg_user_clean_x1d.fits"
    rtdata.get_data(f"nirspec/fs/{user_background}")

    # Get input data
    rtdata.get_data("nirspec/fs/nrs_sci+bkg_cal.fits")

    MasterBackgroundStep.call(rtdata.input, save_results=True, suffix='master_background',
                              user_background=user_background)

    output = "nrs_sci+bkg_master_background.fits"
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_fs_mbkg_user/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nirspec_ifu_mbkg_user(rtdata, fitsdiff_default_kwargs):
    """Test NIRSpec IFU data with a user-supplied background file."""
    # Get user-supplied background
    user_background = "prism_bkg_x1d.fits"
    rtdata.get_data(f"nirspec/ifu/{user_background}")

    # Get input data
    rtdata.get_data("nirspec/ifu/prism_sci_bkg_cal.fits")

    MasterBackgroundStep.call(rtdata.input, user_background=user_background,
                              save_results=True, suffix='master_background')

    output = "prism_sci_bkg_master_background.fits"
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_ifu_mbkg_user/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    'output_file',
    ['ifu_prism_source_on_NRS1_master_background.fits',
     'ifu_prism_source_off_NRS1_o001_masterbg1d.fits',
     'ifu_prism_source_on_NRS1_o001_masterbg2d.fits'],
    ids=["on-source", "off-source", "on-source2d"]
)
def test_nirspec_ifu_mbkg_nod(rtdata, fitsdiff_default_kwargs, output_file):
    """Test NIRSpec IFU prism nodded data."""
    # Get input data
    rtdata.get_asn("nirspec/ifu/nirspec_spec3_asn.json")

    MasterBackgroundStep.call(rtdata.input, save_background=True, save_results=True,
                              suffix='master_background')

    rtdata.output = output_file

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_ifu_mbkg_nod/{output_file}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
