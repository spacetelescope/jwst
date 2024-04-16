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


@pytest.mark.parametrize("inplace", [True, False])
def test_rename_catalog_columns(dummy_source_catalog, inplace):
    renamed_catalog = tweakreg_step._rename_catalog_columns(dummy_source_catalog)

    # if testing inplace, check the input catalog
    if inplace:
        catalog = dummy_source_catalog
    else:
        catalog = renamed_catalog

    assert 'xcentroid' not in catalog.colnames
    assert 'ycentroid' not in catalog.colnames
    assert 'y' in catalog.colnames
    assert 'y' in catalog.colnames


@pytest.mark.parametrize("missing", ["x", "y", "xcentroid", "ycentroid"])
def test_rename_catalog_columns_invalid(dummy_source_catalog, missing):
    # if the column we want to remove is not in the table, first run
    # rename to rename columns this should add the column we want to remove
    if missing not in dummy_source_catalog.colnames:
        tweakreg_step._rename_catalog_columns(dummy_source_catalog)
    dummy_source_catalog.remove_column(missing)
    with pytest.raises(ValueError, match="catalogs must contain"):
        tweakreg_step._rename_catalog_columns(dummy_source_catalog)


def test_rename_catalog_columns_inplace(dummy_source_catalog):
    catalog = tweakreg_step._rename_catalog_columns(dummy_source_catalog)
    assert 'xcentroid' not in catalog.colnames
    assert 'ycentroid' not in catalog.colnames
    assert 'y' in catalog.colnames
    assert 'y' in catalog.colnames


@pytest.mark.skip(reason="test will need to be refactored as _is_wcs_correction_small is gone in favor of comparing skycoords")
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


def test_expected_failure_bad_starfinder():

    model = ImageModel()
    with pytest.raises(ValueError):
        tweakreg_catalog.make_tweakreg_catalog(model, 5.0, bkg_boxsize=400, starfinder='bad_value')


def test_write_catalog(dummy_source_catalog, tmp_cwd):
    '''
    Covers an issue where catalog write did not respect self.output_dir
    '''

    OUTDIR = 'outdir'
    step = tweakreg_step.TweakRegStep()
    os.mkdir(OUTDIR)
    step.output_dir = OUTDIR
    expected_outfile = os.path.join(OUTDIR, 'catalog.ecsv')
    step._write_catalog(dummy_source_catalog, 'catalog.ecsv')

    assert os.path.exists(expected_outfile)
