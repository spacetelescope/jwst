# Licensed under a 3-clause BSD style license - see LICENSE.rst

import datetime
import os
import shutil
import tempfile

import numpy as np
from astropy.io import fits
from astropy.time import Time

from asdf.tags.core import HistoryEntry

from .. import DataModel

TMP_FITS = None
TMP_DIR = None


def setup():
    global TMP_DIR, TMP_FITS

    TMP_DIR = tempfile.mkdtemp()
    TMP_FITS = os.path.join(TMP_DIR, 'tmp.fits')


def teardown():
    shutil.rmtree(TMP_DIR)

def test_historylist_methods():
    m = DataModel()
    h1 = m.history

    info = "First entry"
    h1.append(info)
    assert h1 == info, "Append new history entry"

    h2 = m.history
    assert h2 == info, "Two history lists point to the same object"

    assert len(h1) == 1, "Length of a history list"

    entry = h1[0]
    assert entry["description"] == info, "Get history list item"

    info += " for real"
    h1[0] = info
    assert h1 == info, "Set history list item"

    del h1[0]
    assert len(h1) == 0, "Delete history list item"

    info = ("First entry", "Second_entry", "Third entry")
    h1.extend(info)
    assert len(h1) == 3, "Length of extended history list"
    assert h1 == info, "Contents of extended history list"

    for entry, item in zip(h1, info):
        assert entry["description"]  == item, "Iterate over history list"

    h1.clear()
    assert len(h1) == 0, "Clear history list"

def test_history_from_model_to_fits():
    m = DataModel()
    m.history = [HistoryEntry({
        'description': 'First entry',
        'time': Time(datetime.datetime.now())})]
    m.history.append(HistoryEntry({
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

    assert m2.history == [{'description': "First entry"},
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
    assert m.history == [{'description': 'First entry'},
                         {'description': 'Second entry'}]

    del m.history[0]
    m.history.append(HistoryEntry({'description': "Third entry"}))
    assert m.history == [{'description': "Second entry"},
                         {'description': "Third entry"}]

    m.save(TMP_FITS)

    m = DataModel(TMP_FITS)
    assert m.history == [{'description': "Second entry"},
                         {'description': "Third entry"}]

