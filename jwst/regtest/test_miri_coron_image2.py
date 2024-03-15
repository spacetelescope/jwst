"""Regression test of MIRI conronagraphic data through the image2 pipeline,
   including multiple background exposures that have a mixture of NINTS values"""

import pytest
from jwst.regtest import regtestdata as rt

from jwst.stpipe import Step

@pytest.fixture(scope='module')
def run_image2(rtdata_module):
    """Run the calwebb_image2 pipeline"""

    rtdata = rtdata_module
    rtdata.get_asn('miri/coron/jw03254-c1009_20240118t185403_image2_00001_asn.json')

    args = ["calwebb_image2", rtdata.input,
            "--steps.bkg_subtract.save_results=true",
            "--steps.bkg_subtract.save_combined_background=true",
            "--steps.assign_wcs.save_results=true",
            "--steps.flat_field.save_results=true"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["bsubints", "combinedbackground", "assign_wcs", "flat_field", "calints"])
def test_miri_coron_image2(run_image2, fitsdiff_default_kwargs, suffix):
    """Regression test for image2 processing of MIRI coronagraphic data with background exposures"""

    rt.is_like_truth(
        run_image2, fitsdiff_default_kwargs, suffix, 'truth/test_miri_coron_image2'
    )
