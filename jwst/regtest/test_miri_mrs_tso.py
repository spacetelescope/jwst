"""Regression test for MIRI MRS TSO mode"""
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

# Define artifactory source and truth
INPUT_PATH = 'miri/mrs'
TRUTH_PATH = 'truth/test_miri_mrs_tso'


@pytest.fixture(scope='module')
def run_spec2(jail, rtdata_module):
    """Run the Spec2Pipeline on a single exposure"""
    rtdata = rtdata_module
    collect_pipeline_cfgs("config")

    # Setup the inputs
    file_name = 'jw80600018001_02101_00003_mirifushort_rateints.fits'
    rtdata.get_data(INPUT_PATH + '/' + file_name)

    # Run the pipeline
    args = ["config/calwebb_tso-spec2.cfg", rtdata.input,
        '--steps.assign_wcs.save_results=true',
        '--steps.flat_field.save_results=true',
        '--steps.srctype.save_results=true',
        '--steps.fringe.save_results=true',
        '--steps.photom.save_results=true',
        ]

    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix', ['assign_wcs', 'calints', 'flat_field', 'fringe', 'photom', 'srctype'])
def test_spec2(rtdata_module, run_spec2, fitsdiff_default_kwargs, suffix):
    """Test ensuring the calwebb_tso-spec2 is operating appropriately for MIRI MRS TSO data"""
    rtdata = rtdata_module
    output = f"jw80600018001_02101_00003_mirifushort_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"{TRUTH_PATH}/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
