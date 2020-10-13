import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.lib.suffix import replace_suffix
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

ROOT = 'nrs_verify_nrs1'


@pytest.fixture(scope='module')
def run_detector1(rtdata_module):
    """Run NRS_VERIFY through detector1"""
    rtdata = rtdata_module
    rtdata.get_data('nirspec/imaging/' + replace_suffix(ROOT, 'uncal') + '.fits')

    collect_pipeline_cfgs('config')

    args = [
        'config/calwebb_detector1.cfg', rtdata.input,
        '--steps.dq_init.save_results=True',
        '--steps.saturation.save_results=True',
        '--steps.superbias.save_results=True',
        '--steps.refpix.save_results=True',
        '--steps.linearity.save_results=True',
        '--steps.dark_current.save_results=True',
        '--steps.jump.save_results=True',
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope='module')
def run_image2(run_detector1, rtdata_module):
    """Run NRS_VERIFY through image2"""
    rtdata = rtdata_module

    # Get some predefined references due to insufficient coverage
    # from CRDS.
    refs = ['jwst_nirspec_fpa_0005.asdf', 'jwst_nirspec_flat_0061.fits', 'jwst_nirspec_area_0001.fits']
    for ref in refs:
        rtdata.get_data('nirspec/imaging/' + ref)

    rtdata.input = replace_suffix(ROOT, 'rate') + '.fits'

    args = [
        'jwst.pipeline.Image2Pipeline', rtdata.input,
        '--steps.assign_wcs.save_results=true',
        '--steps.assign_wcs.override_fpa=jwst_nirspec_fpa_0005.asdf',
        '--steps.flat_field.save_results=true',
        '--steps.flat_field.override_flat=jwst_nirspec_flat_0061.fits',
        '--steps.flat_field.override_sflat=N/A',
        '--steps.flat_field.override_fflat=N/A',
        '--steps.flat_field.override_dflat=N/A',
        '--steps.photom.save_results=true',
        '--steps.photom.override_area=jwst_nirspec_area_0001.fits',
    ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix', [
        'dq_init', 'saturation', 'superbias', 'refpix', 'linearity',
        'dark_current', 'jump', 'rate',
    ])
def test_verify_detector1(run_detector1, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test results of the detector1 and image2 processing"""
    rtdata = rtdata_module
    rtdata.input = replace_suffix(ROOT, 'uncal') + '.fits'
    output = replace_suffix(ROOT, suffix) + '.fits'
    rtdata.output = output
    rtdata.get_truth('truth/test_nirspec_verify/' + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix', [
        'assign_wcs', 'flat_field', 'cal',
    ])
def test_verify_image2(run_image2, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test results of the detector1 and image2 processing"""
    rtdata = rtdata_module
    rtdata.input = replace_suffix(ROOT, 'uncal') + '.fits'
    output = replace_suffix(ROOT, suffix) + '.fits'
    rtdata.output = output
    rtdata.get_truth('truth/test_nirspec_verify/' + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
