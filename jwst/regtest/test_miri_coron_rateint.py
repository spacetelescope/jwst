"""Regression test of MIRI conronagraphic rateint data"""

import pytest
from jwst.regtest import regtestdata as rt


@pytest.fixture(scope='module')
def run_miri_coron(jail, rtdata_module):
    """Run the miri background step"""

    rtdata = rtdata_module
    rtdata.input = 'miri/coron/jw01386004001_04101_00001_mirimage_rateints.fits'

    # Get input data files ; 4 for background, 1 science
    rtdata.get_data('miri/coron/jw01386029001_03101_00001_mirimage_rateints.fits')
    rtdata.get_data('miri/coron/jw01386029001_02101_00001_mirimage_rateints.fits')
    rtdata.get_data('miri/coron/jw01386028001_03101_00001_mirimage_rateints.fits')
    rtdata.get_data('miri/coron/jw01386028001_02101_00001_mirimage_rateints.fits')
    rtdata.get_data('miri/coron/jw01386004001_04101_00001_mirimage_rateints.fits')

    bkg_list = ["jw01386029001_03101_00001_mirimage_rateints.fits",
                "jw01386029001_02101_00001_mirimage_rateints.fits",
                "jw01386028001_03101_00001_mirimage_rateints.fits",
                "jw01386028001_02101_00001_mirimage_rateints.fits"]

    # Run the step
    step_params = {
        'step': 'jwst.background.BackgroundStep',
        'args': [
            bkg_list,
            '--save_results=True',
            '--save_combined_background=True',
        ]
    }

    rtdata = rt.run_step_from_dict(rtdata, **step_params)

    return rtdata


@pytest.mark.bigdata
def test_miri_coron(run_miri_coron, fitsdiff_default_kwargs):
    """Regression test for miri_coron applied to MIRI coron rateint data"""

    rt.is_like_truth(
        run_miri_coron, fitsdiff_default_kwargs, 'backgroundstep', 'truth/test_miri_coron_rateint'
    )
