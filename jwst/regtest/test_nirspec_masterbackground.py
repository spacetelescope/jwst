import pytest

from astropy.io.fits.diff import FITSDiff
import numpy as np

import jwst.datamodels as dm
from jwst.master_background import MasterBackgroundNrsSlitsStep
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

from . import regtestdata as rt

pytestmark = pytest.mark.bigdata


@pytest.fixture(scope='module')
def run_spec2_mbkg(jail, rtdata_module):
    """Run Spec2 on MSA data with background slits"""
    rtdata = rtdata_module

    # Get data
    rtdata.get_data('nirspec/mos/nrs_mos_with_bkgslits_msa.fits')
    rtdata.get_data('nirspec/mos/nrs_mos_with_bkgslits_rate.fits')
    collect_pipeline_cfgs('config')

    # Run the pipeline
    step_params = {
        'step': 'calwebb_spec2.cfg',
        'args': [
            '--steps.master_background.skip=false',
            '--steps.master_background.save_background=true'
        ]
    }
    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


def test_masterbkg_rerun(rtdata):
    """Test to ensure sequential runs of the step are consistent"""
    data = dm.open(rtdata.get_data('nirspec/mos/nrs_mos_with_bkgslits_srctype.fits'))

    mbs = MasterBackgroundNrsSlitsStep()
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
    data = dm.open(rtdata.get_data('nirspec/mos/nrs_mos_with_bkgslits_srctype.fits'))

    mbs = MasterBackgroundNrsSlitsStep()
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
    ['cal', 'masterbg1d', 'masterbg2d']
)
def test_nirspec_spec2_mbkg(suffix, run_spec2_mbkg, fitsdiff_default_kwargs):
    """Run spec2 with master background"""
    rtdata = run_spec2_mbkg
    rt.is_like_truth(rtdata, fitsdiff_default_kwargs, suffix, truth_path='truth/test_nirspec_mos_mbkg_user')


def test_nirspec_fs_mbkg_user(rtdata, fitsdiff_default_kwargs):
    """Run a test for NIRSpec FS data with a user-supplied background file."""

    # Get user-supplied background
    user_background = "v2_nrs_bkg_user_clean_x1d.fits"
    rtdata.get_data(f"nirspec/fs/{user_background}")

    # Get input data
    rtdata.get_data("nirspec/fs/nrs_sci+bkg_cal.fits")

    collect_pipeline_cfgs("config")
    args = ["config/master_background.cfg", rtdata.input,
            "--user_background", user_background]
    Step.from_cmdline(args)

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

    collect_pipeline_cfgs("config")
    args = ["config/master_background.cfg", rtdata.input,
            "--user_background", user_background]
    Step.from_cmdline(args)

    output = "prism_sci_bkg_master_background.fits"
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_ifu_mbkg_user/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nirspec_mos_mbkg_user(rtdata, fitsdiff_default_kwargs):
    """Test NIRSpec MOS data with a user-supplied background file."""
    # Get user-supplied background
    user_background = "v2_nrs_mos_bkg_x1d.fits"
    rtdata.get_data(f"nirspec/mos/{user_background}")

    # Get input data
    rtdata.get_data("nirspec/mos/nrs_mos_sci+bkg_cal.fits")

    collect_pipeline_cfgs("config")
    args = ["config/master_background.cfg", rtdata.input,
            "--user_background", user_background]
    Step.from_cmdline(args)

    output = "nrs_mos_sci+bkg_master_background.fits"
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_mos_mbkg_user/{output}")

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

    collect_pipeline_cfgs("config")
    args = ["config/master_background.cfg", rtdata.input,
            "--save_background=True"]
    Step.from_cmdline(args)

    rtdata.output = output_file

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_ifu_mbkg_nod/{output_file}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
