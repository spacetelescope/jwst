"""Regression tests for NIRSpec IFU"""
import pytest

from jwst.regtest import regtestdata as rt

# Define artifactory source and truth
INPUT_PATH = 'nirspec/ifu'
TRUTH_PATH = 'truth/test_nirspec_ifu'


@pytest.fixture(scope='module')
def run_spec2(jail, rtdata_module):
    """Run the Spec2Pipeline on a spec2 ASN containing a single exposure"""
    rtdata = rtdata_module

    # Setup the inputs
    asn_name = 'jw01251-o004_20221105t023552_spec2_031_asn.json'
    asn_path = INPUT_PATH + '/' + asn_name

    # Run the pipeline
    step_params = {
        'input_path': asn_path,
        'step': 'calwebb_spec2',
        'args': [
            '--steps.assign_wcs.save_results=true',
            '--steps.msa_flagging.save_results=true',
            '--steps.srctype.save_results=true',
            '--steps.flat_field.save_results=true',
            '--steps.pathloss.save_results=true',
        ]
    }

    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ['assign_wcs', 'cal', 'flat_field', 'msa_flagging',
     'pathloss', 's3d', 'srctype', 'x1d']
)
def test_spec2(run_spec2, fitsdiff_default_kwargs, suffix):
    """Regression test matching output files"""
    rt.is_like_truth(run_spec2, fitsdiff_default_kwargs, suffix,
                     truth_path=TRUTH_PATH)
