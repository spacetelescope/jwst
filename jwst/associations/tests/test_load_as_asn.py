"""Test LoadAsAssociation"""

import pytest

from . import helpers
from ...datamodels import ImageModel
from ..load_as_asn import (
    LoadAsLevel2Asn
)

DEFAULT_NAME = 'singleton'


@pytest.mark.parametrize(
    'test_input, expected',
    [
        ('file_rate.fits', 'file'),
        ('file_rate_ref.fits', 'file_rate_ref')
    ]
)
def test_suffix_removal(test_input, expected):
    """Ensure appropriate suffix removal is occurring"""

    asn = LoadAsLevel2Asn.load([test_input])
    assert asn['products'][0]['name'] == expected


def test_lv2_datamodel():
    model = ImageModel()
    model.meta.filename = 'modelfile.fits'
    asn = LoadAsLevel2Asn.load(model)
    assert asn.filename == DEFAULT_NAME
    assert asn['program'] == DEFAULT_NAME
    assert asn['target'] == DEFAULT_NAME
    assert asn['asn_pool'] == DEFAULT_NAME
    assert len(asn['products']) == 1
    assert asn['products'][0]['name'] == 'modelfile'


def test_lv2_asn():
    asn_file = helpers.t_path('data/asn_level2.json')
    asn = LoadAsLevel2Asn.load(asn_file)
    assert asn.filename == asn_file
    assert len(asn['products']) == 6


def test_lv2_obj():
    obj = 'some funny object'
    asn = LoadAsLevel2Asn.load(obj)
    assert asn.filename == DEFAULT_NAME
    assert len(asn['products']) == 1
    assert asn['program'] == DEFAULT_NAME
    assert asn['target'] == DEFAULT_NAME
    assert asn['asn_pool'] == DEFAULT_NAME


def test_lv2_obj_list():
    model = ImageModel()
    model.meta.filename = 'modelfile.fits'
    objs = [
        'some funny object',
        model,
    ]
    asn = LoadAsLevel2Asn.load(objs)
    assert asn.filename == DEFAULT_NAME
    assert len(asn['products']) == 2
    assert asn['program'] == DEFAULT_NAME
    assert asn['target'] == DEFAULT_NAME
    assert asn['asn_pool'] == DEFAULT_NAME
