import os.path as op

from jwst.associations import asn_edit

DATA_DIR = op.abspath(op.join(op.dirname(__file__), 'data'))
JSON_FILE = op.join(DATA_DIR, 'asn_level2.json')


def test_add_asn():
    """Test adding a product to an association"""
    asn = asn_edit.reader(JSON_FILE)
    asn = asn_edit.add(asn, ['test_lrs5_rate.fits'], "science")

    nproducts = len(asn['products'])
    assert nproducts == 6
    found = asn_edit._lookup(asn, 'test_lrs5_rate.fits')
    assert len(found) == nproducts


def test_remove_asn():
    """Test removing a product from an association"""
    for ignore in (False, True):
        asn = asn_edit.reader(JSON_FILE)
        assert len(asn['products']) == 6
        asn = asn_edit.remove(asn, ['test_lrs1_rate.fits'], ignore)

        assert len(asn['products']) == 5
        found = asn_edit._lookup(asn, 'test_lrs1_rate.fits')
        assert len(found) == 0
