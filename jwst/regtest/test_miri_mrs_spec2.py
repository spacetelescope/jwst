"""Regression tests for MIRI MRS modes in calwebb_spec2"""
import pytest

from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose

from stdatamodels.jwst import datamodels

from jwst.stpipe import Step

from jwst.regtest import regtestdata as rt

# Define artifactory source and truth file paths
INPUT_PATH = 'miri/mrs'
TRUTH_PATH = 'truth/test_miri_mrs'


@pytest.fixture(scope='module')
def run_spec2(jail, rtdata_module):
    """Run the Spec2Pipeline on a single exposure"""
    rtdata = rtdata_module

    # Get the input rate file
    rtdata.get_data(INPUT_PATH + '/' + 'jw01024001001_04101_00001_mirifulong_rate.fits')

    # Run the pipeline
    args = ["calwebb_spec2", rtdata.input,
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


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ['assign_wcs', 'cal', 'flat_field', 'fringe', 'photom', 's3d', 'srctype', 'straylight', 'x1d']
)
def test_miri_mrs_spec2(run_spec2, fitsdiff_default_kwargs, suffix, rtdata_module):
    """Regression test of the calwebb_spec2 pipeline on a MIRI MRS long-wave exposure"""

    rtdata = rtdata_module
    rt.is_like_truth(rtdata, fitsdiff_default_kwargs, suffix, truth_path=TRUTH_PATH)


@pytest.mark.slow
@pytest.mark.bigdata
def test_miri_mrs_wcs(run_spec2, fitsdiff_default_kwargs, rtdata_module):

    rtdata = rtdata_module

    # get input assign_wcs and truth file
    output = "jw01024001001_04101_00001_mirifulong_assign_wcs.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_miri_mrs/{output}")

    # Open the output and truth file
    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        x, y = grid_from_bounding_box(im.meta.wcs.bounding_box)
        ra, dec, lam = im.meta.wcs(x, y)
        ratruth, dectruth, lamtruth = im_truth.meta.wcs(x, y)
        assert_allclose(ra, ratruth)
        assert_allclose(dec, dectruth)
        assert_allclose(lam, lamtruth)

        # Test the inverse transform
        xtest, ytest = im.meta.wcs.backward_transform(ra, dec, lam)
        xtruth, ytruth = im_truth.meta.wcs.backward_transform(ratruth, dectruth, lamtruth)
        assert_allclose(xtest, xtruth)
        assert_allclose(ytest, ytruth)
