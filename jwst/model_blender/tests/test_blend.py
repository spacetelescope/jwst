"""Test blend_models"""
import pytest

from astropy.table import Table
import numpy as np

from stdatamodels.schema import walk_schema
from stdatamodels.jwst.datamodels import ImageModel

from jwst.model_blender import blendmeta


# Setup various input meta data
N_MODELS = 3  # Number of input models made. All below lists should be this length.
START_TIMES = [57877.00359994354, 57877.0168373584, 57877.03126958496]
EXP_TIMES = [107.3676, 107.3676, 107.3676]
END_TIMES = [57877.0048426241, 57877.01808003896, 57877.03251226551]
FILENAMES = ['image1_cal.fits', 'image2_cal.fits', 'image3_cal.fits']
DATETIMES = ['2017-11-30T13:52:20.367', '2017-11-11T15:14:29.176',
             '2017-11-11T15:15:06.118']
DATES = ['2017-11-30', '2017-11-11', '2017-12-10']
INSTRUMENT_NAMES = ['NIRCAM'] * 3
CORONAGRAPHS = ['4QPM', '4QPM_1065', '4QPM_1140']


def _make_data():
    """Create a set of input models to blendmeta

    Returns
    -------
    (models, input_values, output_values): 3-tuple
        Return a 3-tuple of:
        - models: [model[,...]]
          List of `DataModel` with various meta data.
        - input_values: dict
          Input meta data.
        - output_values: dict
          Expected output values of the blended meta data.
    """
    models = [ImageModel() for i in range(N_MODELS)]

    input_values = {
        'meta.exposure.start_time': START_TIMES,
        'meta.exposure.exposure_time': EXP_TIMES,
        'meta.exposure.end_time': END_TIMES,
        'meta.filename': FILENAMES,
        'meta.instrument.coronagraph': CORONAGRAPHS,
        'meta.instrument.name': INSTRUMENT_NAMES,
        'meta.date': DATETIMES,
        'meta.observation.date': DATES,
        'meta.observation.date_beg': DATETIMES,
    }
    output_values = {
        'meta.exposure.start_time': START_TIMES[0],
        'meta.exposure.exposure_time': np.sum(EXP_TIMES),
        'meta.exposure.end_time': END_TIMES[-1],
        'meta.filename': FILENAMES[0],
        'meta.instrument.coronagraph': CORONAGRAPHS[0],
        'meta.instrument.name': INSTRUMENT_NAMES[0],
        'meta.date': DATETIMES[0],
        'meta.observation.date': DATES[1],
        'meta.observation.date_beg': DATETIMES[1]
    }

    for i, model in enumerate(models):
        for attr in input_values:
            model[attr] = input_values[attr][i]

    return models, input_values, output_values


# The fixture decorator is separated from the actual
# function to allow access to `_make_data` outside the
# `pytest` framework.
@pytest.fixture(scope='module')
def make_data():
    return _make_data()


@pytest.fixture(scope='module')
def blend(make_data):
    """Blend the meta data

    Parameters
    ----------
    make_data: (models, input_values, output_values)
        Results either from the `pytest.fixture.make_data`
        or from the results of `_make_data` directly.
    """
    models, input_values, output_values = make_data
    newmeta, newtab = blendmeta.get_blended_metadata(models)
    return newmeta, newtab, models, input_values, output_values


def test_blendmeta(blend):
    """Test blended metadata

    Parameters
    ----------
    blend: (newmeta, newtab, input_values, output_values)
        Results from `pytest.fixture.blend`
    """
    newmeta, newtab, models, input_values, output_values = blend

    for attr in input_values:
        assert newmeta[attr] == output_values[attr]


def build_fits_dict(schema):
    """
    Utility function to create a dict that maps FITS keywords to their
    metadata attribute in a input schema.

    Parameters
    ----------
    schema : JSON schema fragment
        The schema in which to search.

    Returns
    -------
    results : dict
        Dictionary with FITS keywords as keys and schema metadata
        attributes as values

    """
    def build_fits_dict(subschema, path, combiner, ctx, recurse):
        if len(path) and path[0] == 'extra_fits':
            return True
        kw = subschema.get('fits_keyword')
        if kw is not None:
            results[kw] = '.'.join(path)

    results = {}
    walk_schema(schema, build_fits_dict, results)

    return results


def test_blendtab(blend):
    """Test blended table

    Parameters
    ----------
    blend: (newmeta, newtab, input_values, output_values)
        Results from `pytest.fixture.blend`
    """
    newmeta, newtab, models, input_values, output_values = blend

    # Since the table is of FITS keywords, the meta-to-FITS mapping
    # needs to be determined.
    fits_to_meta = build_fits_dict(models[0].schema)
    meta_to_fits = dict(map(reversed, fits_to_meta.items()))
    fits_expected = set(
        meta_to_fits[meta]
        for meta in input_values
    )

    # Ensure all the expected FITS keywords are in the table.
    table = Table(newtab.data)
    assert not fits_expected.difference(table.colnames)
