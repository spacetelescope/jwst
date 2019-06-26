import datetime

import numpy as np
from astropy.io import fits
from astropy.time import Time

from asdf.tags.core import HistoryEntry

from .. import DataModel


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

def test_history_from_model_to_fits(tmpdir):
    tmpfits = str(tmpdir.join('tmp.fits'))
    m = DataModel()
    m.history = [HistoryEntry({
        'description': 'First entry',
        'time': Time(datetime.datetime.now())})]
    m.history.append(HistoryEntry({
        'description': 'Second entry',
        'time': Time(datetime.datetime.now())
    }))
    m.save(tmpfits)

    with fits.open(tmpfits, memmap=False) as hdulist:
        assert list(hdulist[0].header['HISTORY']) == ["First entry",
                                                      "Second entry"]

    with DataModel(tmpfits) as m2:
        m2 = DataModel()
        m2.update(m)
        m2.history = m.history

        assert m2.history == [{'description': "First entry"},
                              {'description': "Second entry"}]

        m2.save(tmpfits)

    with fits.open(tmpfits, memmap=False) as hdulist:
        assert list(hdulist[0].header['HISTORY']) == ["First entry",
                                                      "Second entry"]


def test_history_from_fits(tmpdir):
    tmpfits = str(tmpdir.join('tmp.fits'))
    header = fits.Header()
    header['HISTORY'] = "First entry"
    header['HISTORY'] = "Second entry"
    fits.writeto(tmpfits, np.array([]), header, overwrite=True)

    with DataModel(tmpfits) as m:
        assert m.history == [{'description': 'First entry'},
                             {'description': 'Second entry'}]

        del m.history[0]
        m.history.append(HistoryEntry({'description': "Third entry"}))
        assert m.history == [{'description': "Second entry"},
                             {'description': "Third entry"}]
        m.save(tmpfits)

    with DataModel(tmpfits) as m:
        assert m.history == [{'description': "Second entry"},
                             {'description': "Third entry"}]
