"""Regression test for NIRSpec imprint subtraction"""
import pytest

from jwst.regtest import regtestdata as rt

# Define artifactory source and truth
INPUT_PATH = 'nirspec/imprint'
TRUTH_PATH = 'truth/test_nirspec_imprint'


@pytest.fixture(scope='module')
def run_spec2(jail, rtdata_module):
    """Run the Spec2Pipeline on a spec2 ASN containing a single exposure
    with multiple imprint exposures"""
    rtdata = rtdata_module

    # Setup the inputs
    asn_name = 'jw01802-o011_spec2_00046_asn.json'
    asn_path = INPUT_PATH + '/' + asn_name

    # Run the pipeline
    step_params = {
        'input_path': asn_path,
        'step': 'calwebb_spec2',
        'args': [
            '--steps.assign_wcs.save_results=true',
            '--steps.imprint_subtract.save_results=true',
            '--steps.msa_flagging.skip=true',
            '--steps.srctype.skip=true',
            '--steps.flat_field.skip=true',
            '--steps.pathloss.skip=true',
            '--steps.photom.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
        ]
    }

    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ['assign_wcs', 'imprint_subtract', 'cal']
)
def test_nirspec_imprint(run_spec2, fitsdiff_default_kwargs, suffix):
    """Regression test matching output files"""
    rt.is_like_truth(run_spec2, fitsdiff_default_kwargs, suffix,
                     truth_path=TRUTH_PATH)

