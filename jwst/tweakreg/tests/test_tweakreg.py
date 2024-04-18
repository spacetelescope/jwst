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


N_EXAMPLE_SOURCES = 21
N_CUSTOM_SOURCES = 15


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

    class FakeCorrector:
        def __init__(self, wcs, original_skycoord):
            self.wcs = wcs
            self._original_skycoord = original_skycoord

        @property
        def meta(self):
            return {'original_skycoord': self._original_skycoord}

    correctors = [FakeCorrector(twcs, tweakreg_step._wcs_to_skycoord(wcs))]

    assert step._is_wcs_correction_small(correctors) == is_good


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
    n_sources = N_EXAMPLE_SOURCES  # a few more than default minobj
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
    m0.meta.filename = 'some_file_0.fits'
    m1.meta.filename = 'some_file_1.fits'
    c = ModelContainer([m0, m1])
    return c


def test_tweakreg_step(example_input):
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


@pytest.fixture()
def custom_catalog_path(tmp_path):
    fn = tmp_path / "custom_catalog.ecsv"

    # it's important that the sources here don't match
    # those added in example_input but conform to the input
    # shape, wcs, etc used in example_input
    rng = np.random.default_rng(42)
    n_sources = N_CUSTOM_SOURCES
    xs = rng.choice(50, n_sources, replace=False) * 8 + 10
    ys = rng.choice(50, n_sources, replace=False) * 8 + 10
    catalog = Table(np.vstack((xs, ys)).T, names=['x', 'y'], dtype=[float, float])
    catalog.write(fn)
    return fn


@pytest.mark.parametrize(
    "catfile",
    ["skip", "valid", "invalid", "empty"],
    ids=["catfile_skip", "catfile_valid", "catfile_invalid", "catfile_empty"],
)
@pytest.mark.parametrize(
    "asn",
    ["skip", "valid", "empty"],
    ids=["asn_skip", "asn_valid", "asn_empty"],
)
@pytest.mark.parametrize(
    "meta",
    ["skip", "valid", "empty"],
    ids=["meta_skip", "meta_valid", "meta_empty"],
)
@pytest.mark.parametrize("custom", [True, False])
@pytest.mark.slow
def test_custom_catalog(custom_catalog_path, example_input, catfile, asn, meta, custom, monkeypatch):
    """
    Options:
        if use_custom_catalogs is False, don't use a catalog
        if use_custom_catalogs is True...
            if catfile is defined...
                if catfile loads -> use_custom_catalogs (ignore asn table, but not model.tweakreg_catalog?)
                if catfile fails to load -> warn and disable custom catalogs
            if catfile is not defined...
                if input doesn't have an asn table...
                    if model has tweakreg_catalog, use it
                if input has an asn table...
                    if member has a tweakreg_catalog use it
                    if not, use meta.tweakreg_catalog
                    if member doesn't have a tweakreg_catalog don't use a custom catalog

    If the entry in either catfile or asn is "", don't use a custom catalog.

    Inputs:
        use_custom_catalogs: True/False
        catfile: missing, non-valid file (load attempted), valid file
        asn_file/table: with tweakreg_catalog, without tweakreg_catalog
        model.meta.tweakreg_catalog (per model)

    Could run step with save_catalog to generate a catalog... or just make a fake one
    """
    example_input[0].meta.group_id = 'a'
    example_input[1].meta.group_id = 'b'

    # this worked because if use_custom_catalogs was true but
    # catfile was blank tweakreg still uses custom catalogs
    # which in this case is defined in model.meta.tweakreg_catalog
    if meta == "valid":
        example_input[0].meta.tweakreg_catalog = str(custom_catalog_path)
    elif meta == "empty":
        example_input[0].meta.tweakreg_catalog = ""

    # write out the ModelContainer and association (so the association table will be loaded)
    example_input.save(dir_path=str(custom_catalog_path.parent))
    asn_data = {
        'asn_id': 'foo',
        'asn_pool': 'bar',
        'products': [
            {
                'members': [{'expname': m.meta.filename, 'exptype': 'science'} for m in example_input],
            },
        ],
    }

    if asn == "empty":
        asn_data['products'][0]['members'][0]['tweakreg_catalog'] = ''
    elif asn == "valid":
        asn_data['products'][0]['members'][0]['tweakreg_catalog'] = str(custom_catalog_path.name)

    import json
    asn_path = custom_catalog_path.parent / 'example_input.json'
    with open(asn_path, 'w') as f:
        json.dump(asn_data, f)

    # write out a catfile
    if catfile != "skip":
        catfile_path = custom_catalog_path.parent / 'catfile.txt'
        with open(catfile_path, 'w') as f:
            if catfile == "valid":
                f.write(f"{example_input[0].meta.filename} {custom_catalog_path.name}")
            elif catfile == "empty":
                f.write(f"{example_input[0].meta.filename}")
            elif catfile == "invalid":
                pass

    # figure out how many sources to expect for the model in group 'a' 
    n_sources = N_EXAMPLE_SOURCES
    custom_number = N_CUSTOM_SOURCES
    if not custom:
        # if use_custom_catalog is False, the custom catalog shouldn't be used
        n_custom_sources = n_sources
    else:
        if catfile == "valid":
            # for a 'valid' catfile, expect the custom number
            n_custom_sources = custom_number
        elif catfile == "invalid":
            # for an 'invalid' catfile, use_custom_catalog should become disabled
            n_custom_sources = n_sources
        elif catfile == "empty":
            # for a catfile with an 'empty' entry, no custom catalog should be used
            n_custom_sources = n_sources
        else:  # catfile == "skip"
            assert catfile == "skip"  # sanity check
            # since catfile is not defined, now look at asn_
            if asn == "valid":
                # for a 'valid' asn entry, expect the custom number
                n_custom_sources = custom_number
            elif asn == "empty":
                # for a 'empty' asn entry, no custom catalog should be used
                n_custom_sources = n_sources
            else:  # asn == "skip"
                assert asn == "skip"  # sanity check
                if meta == "valid":
                    n_custom_sources = custom_number
                elif meta == "empty":
                    n_custom_sources = n_sources
                else:  # meta == "skip"
                    assert meta == "skip"
                    n_custom_sources = n_sources

    kwargs = {'use_custom_catalogs': custom}
    if catfile != "skip":
        kwargs["catfile"] = str(catfile_path)
    step = tweakreg_step.TweakRegStep(**kwargs)

    # patch _construct_wcs_corrector to check the correct catalog was loaded
    def patched_construct_wcs_corrector(model, catalog, _seen=[]):
        # we don't need to continue
        if model.meta.group_id == 'a':
            assert len(catalog) == n_custom_sources
        elif model.meta.group_id == 'b':
            assert len(catalog) == n_sources
        _seen.append(model)
        if len(_seen) == 2:
            raise ValueError("done testing")
        return None

    monkeypatch.setattr(tweakreg_step, "_construct_wcs_corrector", patched_construct_wcs_corrector)

    with pytest.raises(ValueError, match="done testing"):
        step(str(asn_path))
