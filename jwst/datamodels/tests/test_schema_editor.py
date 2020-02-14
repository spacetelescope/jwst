"""Tests for the schema editor

Notes
-----
The tests do not have coverage. The original code was written without tests.
The tests here are for later modifiations.
"""
from os import path as ospath
import pytest

from ..schema_editor import Schema_editor


def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = ospath.dirname(__file__)
    return ospath.join(test_dir, partial_path)


def test_just_run(_jail):
    """Just run the editor"""

    keyword_db = t_path(ospath.join('data', 'jwstkd'))
    editor = Schema_editor(
        input=keyword_db,
        add=True, delete=True, edit=True, rename=True, list=True
    )
    editor.change()
