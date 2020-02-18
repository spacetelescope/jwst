"""Tests for the schema editor

Notes
-----
The tests do not have coverage. The original code was written without tests.
The tests here are for later modifications.
"""
from os import path as ospath
import pytest

from ..schema_editor import Schema_editor


def t_path(partial_path):
    """Construct the full path for test files"""
    test_dir = ospath.dirname(__file__)
    return ospath.join(test_dir, partial_path)


@pytest.fixture
def keyword_db():
    """Define the keyword database"""
    keyword_db = t_path(ospath.join('data', 'jwstkd'))
    return keyword_db


def test_just_run(_jail, keyword_db):
    """Just run the editor"""

    editor = Schema_editor(
        input=keyword_db,
        add=True, delete=True, edit=True, rename=True, list=True
    )
    editor.change()


def test_no_option_warning(_jail, keyword_db):
    """If no operations are given, warn"""
    editor = Schema_editor(
        input=keyword_db
    )
    with pytest.raises(RuntimeError):
        editor.change()
