import json
import os
from copy import deepcopy

import asdf
import numpy as np
import pytest
from astropy.wcs import WCS
from astropy.modeling.models import Shift
from astropy.table import Table
from gwcs.wcstools import grid_from_bounding_box
from stdatamodels.jwst.datamodels import ImageModel

from jwst.datamodels import ModelContainer
from jwst.tweakreg import tweakreg_step
from jwst.tweakreg import tweakreg_catalog
from jwst.tweakreg.utils import _wcsinfo_from_wcs_transform


BKG_LEVEL = 0.001
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
    """
    Test that a catalog with 'xcentroid' and 'ycentroid' columns
    passed to _renamed_catalog_columns successfully renames those columns
    to 'x' and 'y' (and does so "inplace" modifying the input catalog)
    """
    renamed_catalog = tweakreg_step._rename_catalog_columns(dummy_source_catalog)

    # if testing inplace, check the input catalog
    if inplace:
        catalog = dummy_source_catalog
    else:
        catalog = renamed_catalog

    assert 'xcentroid' not in catalog.colnames
    assert 'ycentroid' not in catalog.colnames
    assert 'x' in catalog.colnames
    assert 'y' in catalog.colnames


@pytest.mark.parametrize("missing", ["x", "y", "xcentroid", "ycentroid"])
def test_rename_catalog_columns_invalid(dummy_source_catalog, missing):
    """
    Test that passing a catalog that is missing either "x" or "y"
    (or "xcentroid" and "ycentroid" which is renamed to "x" or "y")
    results in an exception indicating that a required column is missing
    """
    # if the column we want to remove is not in the table, first run
    # rename to rename columns this should add the column we want to remove
    if missing not in dummy_source_catalog.colnames:
        tweakreg_step._rename_catalog_columns(dummy_source_catalog)
    dummy_source_catalog.remove_column(missing)
    with pytest.raises(ValueError, match="catalogs must contain"):
        tweakreg_step._rename_catalog_columns(dummy_source_catalog)


@pytest.mark.parametrize("offset, is_good", [(1 / 3600, True), (11 / 3600, False)])
def test_is_wcs_correction_small(offset, is_good):
    """
    Test that the _is_wcs_correction_small method returns True for a small
    wcs correction and False for a "large" wcs correction. The values in this
    test are selected based on the current step default parameters:
        - use2dhist
        - searchrad
        - tolerance
    Changes to the defaults for these parameters will likely require updating the
    values uses for parametrizing this test.
    """
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
    m0.data[:] = BKG_LEVEL
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


@pytest.mark.parametrize("with_shift", [True, False])
def test_tweakreg_step(example_input, with_shift):
    """
    A simplified unit test for basic operation of the TweakRegStep
    when run with or without a small shift in the input image sources
    """
    if with_shift:
        # shift 9 pixels so that the sources in one of the 2 images
        # appear at different locations (resulting in a correct wcs update)
        example_input[1].data[:-9] = example_input[1].data[9:]
        example_input[1].data[-9:] = BKG_LEVEL

    # assign images to different groups (so they are aligned to each other)
    example_input[0].meta.group_id = 'a'
    example_input[1].meta.group_id = 'b'

    # make the step with default arguments
    step = tweakreg_step.TweakRegStep()

    # run the step on the example input modified above
    result = step(example_input)

    # check that step completed
    for model in result:
        assert model.meta.cal_step.tweakreg == 'COMPLETE'

    # and that the wcses differ by a small amount due to the shift above
    # by projecting one point through each wcs and comparing the difference
    abs_delta = abs(result[1].meta.wcs(0, 0)[0] - result[0].meta.wcs(0, 0)[0])
    if with_shift:
        assert abs_delta > 1E-5
    else:
        assert abs_delta < 1E-12


@pytest.mark.parametrize("alignment_type", ['', 'abs_'])
def test_src_confusion_pars(example_input, alignment_type):
    # assign images to different groups (so they are aligned to each other)
    example_input[0].meta.group_id = 'a'
    example_input[1].meta.group_id = 'b'

    # make the step with arguments that may cause source confusion in match
    pars = {
        f"{alignment_type}separation": 1.0,
        f"{alignment_type}tolerance": 1.0,
    }
    step = tweakreg_step.TweakRegStep(**pars)
    result = step(example_input)

    # check that step was skipped
    for model in result:
        assert model.meta.cal_step.tweakreg == 'SKIPPED'


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
    ["no_catfile", "valid_catfile", "invalid_catfile", "empty_catfile_row"],
)
@pytest.mark.parametrize(
    "asn",
    ["no_cat_in_asn", "cat_in_asn", "empty_asn_entry"],
)
@pytest.mark.parametrize(
    "meta",
    ["no_meta", "cat_in_meta", "empty_meta"],
)
@pytest.mark.parametrize("custom", [True, False])
@pytest.mark.slow
def test_custom_catalog(custom_catalog_path, example_input, catfile, asn, meta, custom, monkeypatch):
    """
    Test that TweakRegStep uses a custom catalog provided by the user
    when the correct set of options are provided. The combinations here can be confusing
    and this test attempts to test all likely combinations of:
        - a catalog in a `catfile`
        - a catalog in the asn
        - a catalog in the metadata
    combined with step options:
        - `use_custom_catalogs` (True/False)
        - a "valid" file passed as `catfile`
    """
    example_input[0].meta.group_id = 'a'
    example_input[1].meta.group_id = 'b'

    # this worked because if use_custom_catalogs was true but
    # catfile was blank tweakreg still uses custom catalogs
    # which in this case is defined in model.meta.tweakreg_catalog
    if meta == "cat_in_meta":
        example_input[0].meta.tweakreg_catalog = str(custom_catalog_path)
    elif meta == "empty_meta":
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

    if asn == "empty_asn_entry":
        asn_data['products'][0]['members'][0]['tweakreg_catalog'] = ''
    elif asn == "cat_in_asn":
        asn_data['products'][0]['members'][0]['tweakreg_catalog'] = str(custom_catalog_path.name)

    asn_path = custom_catalog_path.parent / 'example_input.json'
    with open(asn_path, 'w') as f:
        json.dump(asn_data, f)

    # write out a catfile
    if catfile != "no_catfile":
        catfile_path = custom_catalog_path.parent / 'catfile.txt'
        with open(catfile_path, 'w') as f:
            if catfile == "valid_catfile":
                f.write(f"{example_input[0].meta.filename} {custom_catalog_path.name}")
            elif catfile == "empty_catfile_row":
                f.write(f"{example_input[0].meta.filename}")
            elif catfile == "invalid_catfile":
                pass

    # figure out how many sources to expect for the model in group 'a'
    n_custom_sources = N_EXAMPLE_SOURCES
    if custom:
        if catfile == "valid_catfile":
            # for a 'valid' catfile, expect the custom number
            n_custom_sources = N_CUSTOM_SOURCES
        elif catfile == "no_catfile":
            # since catfile is not defined, now look at asn_
            if asn == "cat_in_asn":
                # for a 'valid' asn entry, expect the custom number
                n_custom_sources = N_CUSTOM_SOURCES
            elif asn == "no_cat_in_asn" and meta == "cat_in_meta":
                n_custom_sources = N_CUSTOM_SOURCES

    kwargs = {'use_custom_catalogs': custom}
    if catfile != "no_catfile":
        kwargs["catfile"] = str(catfile_path)
    step = tweakreg_step.TweakRegStep(**kwargs)

    # patch _construct_wcs_corrector to check the correct catalog was loaded
    def patched_construct_wcs_corrector(model, catalog, _seen=[]):
        # we don't need to continue
        if model.meta.group_id == 'a':
            assert len(catalog) == n_custom_sources
        elif model.meta.group_id == 'b':
            assert len(catalog) == N_EXAMPLE_SOURCES
        _seen.append(model)
        if len(_seen) == 2:
            raise ValueError("done testing")
        return None

    monkeypatch.setattr(tweakreg_step, "_construct_wcs_corrector", patched_construct_wcs_corrector)

    with pytest.raises(ValueError, match="done testing"):
        step(str(asn_path))

@pytest.mark.parametrize("with_shift", [True, False])
def test_sip_approx(example_input, with_shift):
    """
    Test the output FITS WCS.
    """
    if with_shift:
        # shift 9 pixels so that the sources in one of the 2 images
        # appear at different locations (resulting in a correct wcs update)
        example_input[1].data[:-9] = example_input[1].data[9:]
        example_input[1].data[-9:] = BKG_LEVEL

    # assign images to different groups (so they are aligned to each other)
    example_input[0].meta.group_id = 'a'
    example_input[1].meta.group_id = 'b'

    # call th step with override SIP approximation parameters
    step = tweakreg_step.TweakRegStep()
    step.sip_approx = True
    step.sip_max_pix_error = 0.1
    step.sip_degree = 3
    step.sip_max_inv_pix_error = 0.1
    step.sip_inv_degree = 3
    step.sip_npoints = 12

    # run the step on the example input modified above
    result = step(example_input)

    # output wcs differs by a small amount due to the shift above:
    # project one point through each wcs and compare the difference
    abs_delta = abs(result[1].meta.wcs(0, 0)[0] - result[0].meta.wcs(0, 0)[0])
    if with_shift:
        assert abs_delta > 1E-5
    else:
        assert abs_delta < 1E-12

    # the first wcs is identical to the input and
    # does not have SIP approximation keywords --
    # they are normally set by assign_wcs
    assert np.allclose(result[0].meta.wcs(0, 0)[0], example_input[0].meta.wcs(0, 0)[0])
    for key in ['ap_order', 'bp_order']:
        assert key not in result[0].meta.wcsinfo.instance

    # for the second, SIP approximation should be present
    for key in ['ap_order', 'bp_order']:
        assert result[1].meta.wcsinfo.instance[key] == 3

    # evaluate fits wcs and gwcs for the approximation, make sure they agree
    wcs_info = result[1].meta.wcsinfo.instance
    grid = grid_from_bounding_box(result[1].meta.wcs.bounding_box)
    gwcs_ra, gwcs_dec = result[1].meta.wcs(*grid)
    fits_wcs = WCS(wcs_info)
    fitswcs_res = fits_wcs.pixel_to_world(*grid)

    assert np.allclose(fitswcs_res.ra.deg, gwcs_ra)
    assert np.allclose(fitswcs_res.dec.deg, gwcs_dec)
