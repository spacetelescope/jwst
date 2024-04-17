from copy import deepcopy
import os

import asdf
from astropy.modeling.models import Shift
from astropy.table import Table
import numpy as np
import pytest

from jwst.tweakreg import tweakreg_step
from jwst.tweakreg import tweakreg_catalog
from jwst.tweakreg.utils import _wcsinfo_from_wcs_transform
from stdatamodels.jwst.datamodels import ImageModel
from jwst.datamodels import ModelContainer


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


@pytest.fixture()
def example_wcs():
    path = os.path.join(
        os.path.dirname(__file__),
        "data",
        "nrcb1-wcs.asdf")
    with asdf.open(path, lazy_load=False) as af:
        return af.tree["wcs"]


@pytest.fixture()
def example_input(example_wcs):
    m0 = ImageModel((512, 512))

    # add a wcs and wcsinfo
    m0.meta.wcs = example_wcs
    m0.meta.wcsinfo = _wcsinfo_from_wcs_transform(example_wcs)

    # and a few 'sources'
    m0.data[:] = 0.001
    n_sources = 21  # a few more than default minobj
    rng = np.random.default_rng(26)
    xs = rng.choice(50, n_sources, replace=False) * 8 + 10
    ys = rng.choice(50, n_sources, replace=False) * 8 + 10
    for y, x in zip(ys, xs):
        m0.data[y-1:y+2, x-1:x+2] = [
            [0.1, 0.6, 0.1],
            [0.6, 0.8, 0.6],
            [0.1, 0.6, 0.1],
        ]

    m1 = m0.copy()
    # give each a unique filename
    m0.meta.filename = 'some_file_0'
    m1.meta.filename = 'some_file_1'
    c = ModelContainer([m0, m1])
    return c


def test_run_tweakreg(example_input):
    # shift 9 pixels
    example_input[1].data = np.roll(example_input[1].data, 9, axis=0)

    # assign images to different groups (so they are aligned to each other)
    example_input[0].meta.group_id = 'a'
    example_input[1].meta.group_id = 'b'
    step = tweakreg_step.TweakRegStep()
    result = step(example_input)

    # check that step completed
    for model in result:
        assert model.meta.cal_step.tweakreg == 'COMPLETE'

    # and that the wcses are different
    abs_delta = abs(result[1].meta.wcs(0, 0)[0] - result[0].meta.wcs(0, 0)[0])
    assert abs_delta > 1E-5


def test_abs_refcat():
    pass


def test_custom_catalog():
    """
    Options:
        if use_custom_catalogs is False, don't use a catalog
        if use_custom_catalogs is True...
            if catfile is defined...
                if catfile loads -> use_custom_catalogs (ignore asn table)
                if catfile fails to load -> warn and disable custom catalogs
            if catfile is not defined...
                if input doesn't have an asn table -> disable custom catalogs
                if input has an asn table...
                    if member has a tweakreg_catalog use it
                    if member doesn't have a tweakreg_catalog don't use a custom catalog
    """
    pass
