import os
import tempfile
import os.path as op
from glob import glob
from shutil import copyfile, rmtree

from .. import asn_edit


ROOT_DIR = None
DATA_DIR = None
TMP_DIR = None
JSON_FILE = None
TMP_FITS = None
TMP_FITS2 = None
TMP_JSON = None

def setup():
    global ROOT_DIR, DATA_DIR, TMP_DIR, JSON_FILE
    global TMP_FITS, TMP_FITS2, TMP_JSON, TMP_JSON2
    ROOT_DIR = op.abspath(op.dirname(__file__))
    DATA_DIR = op.join(ROOT_DIR, 'data')
    TMP_DIR = tempfile.mkdtemp()

    JSON_FILE = os.path.join(DATA_DIR, 'asn_level2.json')
    TMP_FITS = os.path.join(TMP_DIR, 'test_lrs1_rate.fits')
    TMP_FITS2 = os.path.join(TMP_DIR, 'test_lrs5_rate.fits')
    TMP_JSON = os.path.join(TMP_DIR, 'asn_level2.json')

    copyfile(JSON_FILE, TMP_JSON)
    for fname in (TMP_FITS, TMP_FITS2):
        with open(fname, 'a'):
            os.utime(fname)
    os.chdir(TMP_DIR)


def teardown():
    os.chdir(ROOT_DIR)
    rmtree(TMP_DIR)


def test_add_asn():
    cmd = ["-a", TMP_JSON, TMP_FITS2]
    asn = asn_edit.reader(TMP_JSON)
    asn = asn_edit.add(asn, [TMP_FITS2], "science")

    nproducts = len(asn['products'])
    assert nproducts == 6, "Add a file"
    found = asn_edit._lookup(asn, "test_lrs5_rate.fits")
    assert len(found) == nproducts, "Add a file"


def test_remove_asn():
    for ignore in (False, True):
        asn = asn_edit.reader(TMP_JSON)
        asn = asn_edit.remove(asn, [TMP_FITS], ignore)

        assert len(asn['products']) == 5, "Remove a file"
        found = asn_edit._lookup(asn, "test_lrs1_rate.fits")
        assert len(found) == 0, "Remove a file"


for test_function in (test_add_asn, test_remove_asn):
    setup()
    test_function()
    teardown()
