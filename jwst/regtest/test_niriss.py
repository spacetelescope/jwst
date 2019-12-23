"""Regression tests for NIRISS"""
import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.lib.suffix import replace_suffix
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


def is_like_truth(rtdata, fitsdiff_default_kwargs, suffix, truth_path='truth/niriss/test_niriss'):
    """Compare step outputs with truth"""

    # If the input was an association, the output should be the name of the product
    # Otherwise, output is based on input.
    if rtdata.asn:
        output = rtdata.asn['products'][0]['name']
    else:
        output = os.path.splitext(os.path.basename(rtdata.input))[0]
    output = replace_suffix(output, suffix) + '.fits'
    rtdata.output = output

    rtdata.get_truth(os.path.join(truth_path, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.fixture(scope='module')
def run_nrs_wfss_spectral(jail, rtdata_module):
    """Run the pipelines"""
    rtdata = rtdata_module
    rtdata.get_asn('niriss/level2a/nir_wfss_spec2_asn.json')

    collect_pipeline_cfgs('config')
    args = [
        'config/calwebb_spec2.cfg',
        rtdata.input,
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
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ['assign_wcs', 'bsub', 'cal', 'extract_2d', 'flat_field', 'photom', 'srctype', 'x1d']
)
def test_nrs_wfss_spectral(run_nrs_wfss_spectral, fitsdiff_default_kwargs, suffix):
    """Regression test matching output files"""
    fitsdiff_default_kwargs['ignore_keywords'].append('SCATFILE')
    is_like_truth(run_nrs_wfss_spectral, fitsdiff_default_kwargs, suffix)
