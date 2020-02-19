"""Tests for the schema editor

Notes
-----
The tests do not have coverage. The original code was written without tests.
The tests here are for later modifications.
"""
import os
from pathlib import Path
import pytest

from jwst.datamodels.schema_editor import Schema_editor
from jwst.lib.file_utils import pushdir

# Define data locations in artifactory
KEYWORD_DB = 'datamodels/keyword_db'

@pytest.fixture(scope='module')
def keyword_db(jail, rtdata_module):
    """Define the keyword database"""
    keyword_db = Path('keyword_db')
    keyword_db.mkdir()
    with pushdir(str(keyword_db)):
        schemas = rtdata_module.data_glob(KEYWORD_DB, glob='*.json')
        for schema in schemas:
            rtdata_module.get_data(schema)

    return keyword_db


def test_just_run(jail, keyword_db):
    """Just run the editor"""

    editor = Schema_editor(
        input=str(keyword_db),
        add=True, delete=True, edit=True, rename=True, list=True
    )
    editor.change()


def test_no_option_warning(jail, keyword_db):
    """If no operations are given, warn"""
    editor = Schema_editor(
        input=keyword_db
    )
    with pytest.raises(RuntimeError):
        editor.change()
