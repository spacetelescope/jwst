import shutil
import tempfile

import numpy as np

from jwst.datamodels import ImageModel

from .. import blendmeta

#ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')
start_times = [57877.00359994354, 57877.0168373584, 57877.03126958496]
exp_times = [107.3676, 107.3676, 107.3676]
end_times = [57877.0048426241, 57877.01808003896, 57877.03251226551]
filenames = ['image1_cal.fits', 'image2_cal.fits', 'image3_cal.fits']
dates = ['2017-11-30T13:52:20.367', '2017-11-11T15:14:29.176',
         '2017-11-11T15:15:06.118']
instrument_names = ['NIRCAM'] * 3


def setup():
    global TMP_DIR, TMP_NAMES, TMP_FILES, INPUT_VALUES, OUTPUT_VALUES

    TMP_DIR = tempfile.mkdtemp()
    TMP_FILES = [ImageModel() for i in range(3)]

    INPUT_VALUES = {'meta.exposure.start_time': start_times,
                    'meta.exposure.exposure_time': exp_times,
                    'meta.exposure.end_time': end_times,
                    'meta.filename': filenames,
                    'meta.instrument.name': instrument_names,
                    'meta.date': dates}
    OUTPUT_VALUES = {'meta.exposure.start_time': start_times[0],
                    'meta.exposure.exposure_time': np.sum(exp_times),
                    'meta.exposure.end_time': end_times[-1],
                    'meta.filename': filenames[0],
                    'meta.instrument.name': instrument_names[0],
                    'meta.date': dates[0]}

    for i, tfile in enumerate(TMP_FILES):
        for attr in INPUT_VALUES:
            tfile[attr] = INPUT_VALUES[attr][i]

def teardown():
    shutil.rmtree(TMP_DIR)


def test_blendmeta():
    newmeta, newtab = blendmeta.get_blended_metadata(TMP_FILES)
    for attr in INPUT_VALUES:
        assert newmeta[attr] == OUTPUT_VALUES[attr]
