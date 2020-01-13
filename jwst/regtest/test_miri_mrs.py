"""Regression tests for MIRI MRS modes"""
from pathlib import Path
import pytest

from jwst.lib.suffix import replace_suffix

from . import regtestdata as rt


@pytest.fixture(scope='module')
def run_spec2(jail, rtdata_module):
    """Run the Spec2Pipeline on a single exposure"""
    rtdata = rtdata_module

    # Setup the inputs
    asn_name = 'ifushort_ch12_rate_asn3.json'
    rtdata.get_asn('miri/mrs/' + asn_name, get_members=False)
    asn_path = rtdata.input
    member_path = Path(rtdata.asn['products'][0]['members'][0]['expname'])
    rate_path = member_path.stem
    rate_path = replace_suffix(rate_path, 'rate')
    rate_path = 'miri/mrs/' + rate_path + member_path.suffix

    # Run the pipeline
    step_params = {
        'input_path': rate_path,
        'step': 'calwebb_spec2.cfg',
        'args': [
            '--steps.bkg_subtract.save_results=true',
            '--steps.assign_wcs.save_results=true',
            '--steps.imprint_subtract.save_results=true',
            '--steps.msa_flagging.save_results=true',
            '--steps.extract_2d.save_results=true',
            '--steps.flat_field.save_results=true',
            '--steps.srctype.save_results=true',
            '--steps.straylight.save_results=true',
            '--steps.fringe.save_results=true',
            '--steps.pathloss.save_results=true',
            '--steps.barshadow.save_results=true',
            '--steps.photom.save_results=true',
            '--steps.resample_spec.save_results=true',
            '--steps.cube_build.save_results=true',
            '--steps.extract_1d.save_results=true',
        ]
    }

    # return rt.run_step_from_dict(rtdata, **step_params)
    rt.run_step_from_dict_mock(
        rtdata,
        '/Users/eisenham/Downloads/artifactory/jwst-pipeline/newdev/truth/test_miri_mrs',
        **step_params
    )

    # Set RT input back to the association path for the next pipeline run.
    rtdata.input = asn_path

    # That's all folks
    return rtdata


@pytest.fixture(scope='module')
def run_spec3(jail, run_spec2):
    """Run the Spec3Pipeline on the results from the Spec2Pipeline run"""
    rtdata = run_spec2

    # The presumption is that `run_spec2` has set the input to the
    # original association. To use this default, and not re-download
    # the association, simply do not specify `step_params["input_path"]`
    step_params = {
        'step': 'calwebb_spec3.cfg',
        'args': [
            '--steps.master_background.save_results=true',
            '--steps.mrs_imatch.save_results=true',
            '--steps.outlier_detection.save_results=true',
            '--steps.resample_spec.save_results=true',
            '--steps.cube_build.save_results=true',
            '--steps.extract_1d.save_results=true',
            '--steps.combine_1d.save_results=true',
        ]
    }

    return rt.run_step_from_dict(rtdata, **step_params)


@pytest.fixture(scope='module')
def run_spec3_multi(jail, rtdata_module):
    """Run the Spec3Pipeline on multi channel/multi filter data"""
    step_params = {
        'input_path': 'miri/mrs/ifushort_set2_asn3.json',
        'step': 'calwebb_spec3.cfg',
        'args': [
            '--steps.master_background.save_results=true',
            '--steps.mrs_imatch.save_results=true',
            '--steps.outlier_detection.save_results=true',
            '--steps.resample_spec.save_results=true',
            '--steps.cube_build.save_results=true',
            '--steps.extract_1d.save_results=true',
            '--steps.combine_1d.save_results=true',
        ]
    }

    return rt.run_step_from_dict(rtdata_module, **step_params)

@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ['assign_wcs', 'cal', 'flat_field', 'fringe', 'photom', 's3d', 'srctype', 'straylight', 'x1d']
)
def test_spec2(run_spec2, fitsdiff_default_kwargs, suffix):
    """Test ensuring the callwebb_spec2 is operating appropriately for MIRI MRS data"""
    rt.is_like_truth(run_spec2, fitsdiff_default_kwargs, suffix,
                     truth_path='truth/test_miri_mrs')


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'ifushort_ch12_spec3_ch1-medium_s3d.fits',
        'ifushort_ch12_spec3_ch1-medium_x1d.fits',
        'ifushort_ch12_spec3_ch2-medium_s3d.fits',
        'ifushort_ch12_spec3_ch2-medium_x1d.fits',
        'ifushort_ch12_spec3_mrs_imatch.fits'
    ]
)
def test_spec3(run_spec3, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rt.is_like_truth(
        run_spec3, fitsdiff_default_kwargs, output,
        truth_path='truth/test_miri_mrs',
        is_suffix=False
    )
