# Licensed under a 3-clause BSD style license - see LICENSE.rst

import datetime
import os
import shutil
import tempfile

import numpy as np
from astropy.io import fits
from astropy.time import Time

from asdf.tags.core import HistoryEntry

from jwst.datamodels import DataModel

TMP_FITS = None
TMP_DIR = None


def setup():
    global TMP_DIR, TMP_FITS

    TMP_DIR = tempfile.mkdtemp()
    TMP_FITS = os.path.join(TMP_DIR, 'tmp.fits')


def teardown():
    shutil.rmtree(TMP_DIR)


def test_history_from_model_to_fits():
    m = DataModel()
    m.history['entries']= [HistoryEntry({
        'description': 'First entry',
        'time': Time(datetime.datetime.now())})]
    m.history['entries'].append(HistoryEntry({
        'description': 'Second entry',
        'time': Time(datetime.datetime.now())
    }))
    m.save(TMP_FITS)

    hdulist = fits.open(TMP_FITS)
    assert list(hdulist[0].header['HISTORY']) == ["First entry", "Second entry"]
    hdulist.close()

    m = DataModel(TMP_FITS)
    m2 = DataModel()
    m2.update(m)
    m2.history = m.history

    assert m2.history['entries'] == [{'description': "First entry"},
                                     {'description': "Second entry"}]

    m2.save(TMP_FITS)

    hdulist = fits.open(TMP_FITS)
    assert list(hdulist[0].header['HISTORY']) == ["First entry", "Second entry"]
    hdulist.close()


def test_history_from_fits():
    header = fits.Header()
    header['HISTORY'] = "First entry"
    header['HISTORY'] = "Second entry"
    fits.writeto(TMP_FITS, np.array([]), header, overwrite=True)

    m = DataModel(TMP_FITS)
    assert m.history['entries'] == [{'description': 'First entry'},
                                    {'description': 'Second entry'}]

    del m.history['entries'][0]
    m.history['entries'].append(HistoryEntry({'description': "Third entry"}))
    assert m.history['entries'] == [{'description': "Second entry"},
                                    {'description': "Third entry"}]

    m.save(TMP_FITS)

    m = DataModel(TMP_FITS)
    assert m.history['entries'] == [{'description': "Second entry"},
                                    {'description': "Third entry"}]
