from copy import deepcopy
import os

import asdf
from astropy.modeling.models import Shift
from astropy.table import Table
import pytest

from jwst.tweakreg import tweakreg_step
from jwst.tweakreg import tweakreg_catalog
from stdatamodels.jwst.datamodels import ImageModel


@pytest.fixture
def dummy_source_catalog():

    columns = ['id', 'xcentroid', 'ycentroid', 'flux']
    catalog = Table(names=columns, dtype=(int, float, float, float))
    catalog.add_row([1, 100.0, 100.0, 100.0])

    return catalog


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
        ([['abc1_cal.fits', 'abc2_cal.fits']], {}, ['abc #1']),
        (
            [
                ['abc1_cal.fits', 'abc2_cal.fits'],
                ['abc1_cal.fits', 'abc2_cal.fits'],
                ['abc1_cal.fits', 'abc2_cal.fits'],
                ['def1_cal.fits', 'def2_cal.fits'],
                ['cba1_cal.fits', 'cba2_cal.fits'],
            ],
            {'cba': 1},
            ["abc #1", "abc #2", "abc #3", "def #1", "cba #2"],
        ),
        ([['cba1_cal.fits', 'abc2_cal.fits']], {}, ['Group #1']),
        ([['cba1_cal.fits', 'abc2_cal.fits']], {'Group': 1}, ['Group #2']),
        ([['cba1_cal.fits', 'abc2_cal.fits']], None, ['Group #1']),
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


def test_common_name_special():
    all_group_names = {}
    group = []
    cn = tweakreg_step._common_name(group, all_group_names)
    assert cn == 'Group #1'

    group = []
    all_group_names = {'abc': 1, 'abc #2': 1}
    for fname in ['abc1', 'abc2']:
        model = ImageModel()
        model.meta.filename = fname
        group.append(model)
    cn = tweakreg_step._common_name(group, all_group_names)
    assert len(cn.rsplit('_', 1)[1]) == 6


def test_expected_failure_bad_starfinder():

    model = ImageModel()
    with pytest.raises(ValueError):
        tweakreg_catalog.make_tweakreg_catalog(model, 5.0, bkg_boxsize=400, starfinder='bad_value')


def test_write_catalog(dummy_source_catalog, tmp_cwd):
    '''
    Covers an issue where catalog write did not respect self.output_dir
    '''

    OUTDIR = 'outdir'
    model = ImageModel()
    step = tweakreg_step.TweakRegStep()
    os.mkdir(OUTDIR)
    step.output_dir = OUTDIR
    expected_outfile = os.path.join(OUTDIR, 'catalog.ecsv')
    step._write_catalog(model, dummy_source_catalog, 'catalog.ecsv')

    assert os.path.exists(expected_outfile)