"""Test blend_models"""
import pytest

import numpy as np

from jwst.datamodels import ImageModel

from .. import blendmeta


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


@pytest.fixture(scope='module')
def make_data(jail):
    """Create a set of input models to blendmeta

    Parameters
    ----------
    jail: None
        A `pytest` fixture that changes the current working directory
        to a temporary folder

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
        'meta.instrument.name': INSTRUMENT_NAMES,
        'meta.date': DATETIMES,
        'meta.observation.date': DATES,
        'meta.observation.date_beg': DATETIMES
    }
    output_values = {
        'meta.exposure.start_time': START_TIMES[0],
        'meta.exposure.exposure_time': np.sum(EXP_TIMES),
        'meta.exposure.end_time': END_TIMES[-1],
        'meta.filename': FILENAMES[0],
        'meta.instrument.name': INSTRUMENT_NAMES[0],
        'meta.date': DATETIMES[0],
        'meta.observation.date': DATES[1],
        'meta.observation.date_beg': DATETIMES[1]
    }

    for i, model in enumerate(models):
        for attr in input_values:
            model[attr] = input_values[attr][i]

    return models, input_values, output_values


def test_blendmeta(make_data):
    """Blend data and compare to expected output"""
    models, input_values, output_values = make_data

    newmeta, newtab = blendmeta.get_blended_metadata(models)
    for attr in input_values:
        assert newmeta[attr] == output_values[attr]
