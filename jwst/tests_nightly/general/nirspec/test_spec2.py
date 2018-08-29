"""Test aspects of Spec2Pipline"""
import os
import pytest

from jwst.pipeline import Spec2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_nrs2_nodata(_bigdata):
    """

    Regression test of handling NRS2 detector that has no data.

    """

    # Only need to ensure that assing_wcs is run.
    # This still will fail and should cause the pipeline to halt.
    step = Spec2Pipeline()
    step.assign_wcs.skip = False

    step.run(os.path.join(
        _bigdata, 'nirspec', 'test_assignwcs', 'jw84700006001_02101_00001_nrs2_rate.fits'
    ))

    assert False
