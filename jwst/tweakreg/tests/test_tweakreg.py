from copy import deepcopy
import os

import asdf
from astropy.modeling.models import Shift
import pytest

from jwst.tweakreg import tweakreg_step
from stdatamodels.jwst.datamodels import ImageModel


@pytest.mark.parametrize("offset, is_good", [(1 / 3600, True), (11 / 3600, False)])
def test_is_wcs_correction_small(offset, is_good):
    path = os.path.join(os.path.dirname(__file__), "mosaic_long_i2d_gwcs.asdf")
    with asdf.open(path) as af:
        wcs = af.tree["wcs"]

    # Make a copy and add an offset at the end of the transform
    twcs = deepcopy(wcs)
    step = twcs.pipeline[0]
    step.transform = step.transform | Shift(offset) & Shift(offset)
    twcs.bounding_box = wcs.bounding_box

    step = tweakreg_step.TweakRegStep()

    assert step._is_wcs_correction_small(wcs, twcs) == is_good


@pytest.mark.parametrize(
    "groups, all_group_names, common_name",
    [
        ([['abc1_cal.fits', 'abc2_cal.fits']], [], ['abc']),
        (
            [
                ['abc1_cal.fits', 'abc2_cal.fits'],
                ['abc1_cal.fits', 'abc2_cal.fits'],
                ['abc1_cal.fits', 'abc2_cal.fits'],
                ['def1_cal.fits', 'def2_cal.fits'],
            ],
            [],
            ["abc", "abc1", "abc2", "def"],
        ),
        ([['cba1_cal.fits', 'abc2_cal.fits']], [], ['Group #1']),
        ([['cba1_cal.fits', 'abc2_cal.fits']], ['Group #1'], ['Group #2']),
        ([['cba1_cal.fits', 'abc2_cal.fits']], None, ['Unnamed Group']),
    ]
)
def test_common_name(groups, all_group_names, common_name):
    for g, cn_truth in zip(groups, common_name):
        group = []
        for fname in g:
            model = ImageModel()
            model.meta.filename = fname
            group.append(model)

        cn = tweakreg_step._common_name(group, all_group_names)
        assert cn == cn_truth
