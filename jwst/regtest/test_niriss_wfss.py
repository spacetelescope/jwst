"""Regression tests for NIRISS WFSS mode"""
import pytest
from astropy.io.fits.diff import FITSDiff

from . import regtestdata as rt


@pytest.fixture(scope='module')
def run_nis_wfss_spec2(jail, rtdata_module):
    """Run the calwebb_spec2 pipeline"""
    step_params = {
        'input_path': 'niriss/wfss/nir_wfss_spec2_asn.json',
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

    return rt.run_step_from_dict(rtdata_module, **step_params)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ['assign_wcs', 'bsub', 'cal', 'extract_2d', 'flat_field', 'photom', 'srctype', 'x1d']
)
def test_nis_wfss_spec2(run_nis_wfss_spec2, fitsdiff_default_kwargs, suffix):
    """Regression test for calwebb_spec2 applied to NIRISS WFSS data"""
    rt.is_like_truth(
        run_nis_wfss_spec2, fitsdiff_default_kwargs, suffix, 'truth/test_niriss_wfss'
    )


@pytest.fixture(scope='module')
def run_nis_wfss_spec3(jail, rtdata_module):
    """Run the calwebb_spec3 pipeline"""
    step_params = {
        'input_path': 'niriss/wfss/jw00625-o030_20191121t041727_spec3_001_asn.json',
        'step': 'calwebb_spec3.cfg',
        'args': [
            '--steps.extract_1d.save_results=true',
            '--steps.combine_1d.save_results=true',
        ]
    }

    return rt.run_step_from_dict(rtdata_module, **step_params)


@pytest.mark.bigdata
@pytest.mark.parametrize('suffix', ['cal', 'x1d', 'c1d'])
@pytest.mark.parametrize('source_id', ['s00001', 's00002'])
def test_nis_wfss_spec3(run_nis_wfss_spec3, suffix, source_id, fitsdiff_default_kwargs):
    """Regression test of the calwebb_spec3 pipeline applied to NIRISS WFSS data"""

    # Run the pipeline
    rtdata = run_nis_wfss_spec3

    # Get outputs and truths
    output = "jw00625-o030_" + source_id + "_niriss_f090w-gr150c-gr150r_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_niriss_wfss/" + output)

    # Compare the results
    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
