"""Test aspects of Spec2Pipline"""
import os
import subprocess

import pytest

from jwst.assign_wcs.util import NoDataOnDetectorError
from jwst.pipeline import Spec2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


@pytest.mark.xfail(
    reason='Input data not available',
    run=False
)
def test_nrs2_nodata_api(_bigdata):
    """

    Regression test of handling NRS2 detector that has no data.

    """

    # Only need to ensure that assing_wcs is run.
    # This still will fail and should cause the pipeline to halt.
    step = Spec2Pipeline()
    step.assign_wcs.skip = False

    with pytest.raises(NoDataOnDetectorError):
        step.run(os.path.join(
            _bigdata, 'nirspec', 'test_assignwcs', 'jw84700006001_02101_00001_nrs2_rate.fits'
        ))


@pytest.mark.xfail(
    reason='Input data not available',
    run=False
)
def test_nrs2_nodata_strun(_bigdata):
    """Ensure that the appropriate exit status is returned from strun"""

    data_file = os.path.join(
        _bigdata, 'nirspec', 'test_assignwcs', 'jw84700006001_02101_00001_nrs2_rate.fits'
    )

    cmd = [
        'strun',
        'jwst.pipeline.Spec2Pipeline',
        data_file,
        '--steps.assign_wcs.skip=false'
    ]
    status = subprocess.run(cmd)

    assert status.returncode == 64
