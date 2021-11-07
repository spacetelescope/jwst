"""Regression tests for MIRI MRS modes"""
import os
import pytest

from jwst.regtest import regtestdata as rt

from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step
from jwst.associations.asn_from_list import asn_from_list

# Define artifactory source and truth
INPUT_PATH = 'miri/mrs'
TRUTH_PATH = 'truth/test_miri_mrs'


@pytest.fixture(scope='module')
def run_spec3(jail, rtdata_module):
    """Run the Spec3Pipeline on the single cal result from the Spec2Pipeline run"""
    rtdata = rtdata_module

    # Note that we use the truth file from spec2 processing as the input to spec3
    rtdata.get_data(TRUTH_PATH + '/' + 'ifushort_ch12_cal.fits')

    asn = asn_from_list([rtdata.input], product_name='ifushort_ch12_spec3')
    asn.data["program"] = "00024"
    asn.data["asn_type"] = "spec3"
    asn.sequence = 1
    asn_name, serialized = asn.dump(format="json")
    with open(asn_name, "w") as f:
        f.write(serialized)

    rtdata.input = asn_name

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.master_background.save_results=true',
        '--steps.mrs_imatch.save_results=true',
        '--steps.outlier_detection.save_results=true',
        '--steps.resample_spec.save_results=true',
        '--steps.cube_build.save_results=true',
        '--steps.extract_1d.save_results=true',
        '--steps.combine_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope='module')
def run_spec3_multi(jail, rtdata_module):
    """Run the Spec3Pipeline on multi channel/multi filter data"""

    step_params = {
        'input_path': INPUT_PATH + '/' + 'ifushort_set2_asn3.json',
        'step': 'calwebb_spec3',
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


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'ifushort_ch12_spec3_mrs_imatch.fits',
        'ifushort_ch12_spec3_ch1-medium_s3d.fits',
        'ifushort_ch12_spec3_ch2-medium_s3d.fits',
        'ifushort_ch12_spec3_ch1-medium_x1d.fits',
        'ifushort_ch12_spec3_ch2-medium_x1d.fits',
    ],
    ids=["mrs_imatch", "ch1-s3d", "ch2-s3d", "ch1-x1d", "ch2-x1d"]
)
def test_spec3(run_spec3, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""

    rtdata = run_spec3
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, rtdata.output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'ifushort_set2_0_mrs_imatch.fits',
        'ifushort_set2_1_mrs_imatch.fits',
        'ifushort_set2_0_a3001_crf.fits',
        'ifushort_set2_1_a3001_crf.fits',
        'ifushort_set2_ch1-short_s3d.fits',
        'ifushort_set2_ch2-short_s3d.fits',
        'ifushort_set2_ch1-short_x1d.fits',
        'ifushort_set2_ch2-short_x1d.fits',
    ],
    ids=["ch1-mrs_imatch", "ch2-mrs_imatch", "ch1-crf", "ch2-crf",
         "ch1-s3d", "ch2-s3d", "ch1-x1d", "ch2-x1d"]
)
def test_spec3_multi(run_spec3_multi, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rt.is_like_truth(
        run_spec3_multi, fitsdiff_default_kwargs, output,
        truth_path=TRUTH_PATH,
        is_suffix=False
    )
