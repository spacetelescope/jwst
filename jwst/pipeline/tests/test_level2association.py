"""Test Level2Association"""

from ...associations.tests import helpers as asn_helpers
from ...datamodels import ImageModel
from ..level2association import Level2Association


def test_datamodel():
    model = ImageModel()
    model.meta.filename = 'modelfile.fits'
    asn = Level2Association.open(model)
    assert asn['program'] == 'lvl2asncreated'
    assert asn['target'] == 'lvl2asncreated'
    assert asn['asn_pool'] == 'lvl2asncreated'
    assert len(asn['products']) == 1
    assert asn['products'][0]['name'] == 'modelfile'


def test_asn():
    asn_file = asn_helpers.t_path('data/asn_level2.json')
    asn = Level2Association.open(asn_file)
    assert len(asn['products']) == 6


def test_obj():
    obj = 'some funny object'
    asn = Level2Association.open(obj)
    assert len(asn['products']) == 1
    assert asn['program'] == 'lvl2asncreated'
    assert asn['target'] == 'lvl2asncreated'
    assert asn['asn_pool'] == 'lvl2asncreated'


def test_obj_list():
    model = ImageModel()
    model.meta.filename = 'modelfile.fits'
    objs = [
        'some funny object',
        model,
    ]
    asn = Level2Association.open(objs)
    assert len(asn['products']) == 2
    assert asn['program'] == 'lvl2asncreated'
    assert asn['target'] == 'lvl2asncreated'
    assert asn['asn_pool'] == 'lvl2asncreated'
