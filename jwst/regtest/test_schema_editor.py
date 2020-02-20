"""Tests for the schema editor

Notes
-----
The tests do not have coverage. The original code was written without tests.
The tests here are for later modifications.
"""
import os
from pathlib import Path
import sys

import pytest

from jwst.regtest.regtestdata import text_diff
import jwst.datamodels.schema_editor as se
from jwst.lib.file_utils import pushdir

# Define data locations in artifactory
KEYWORD_DB = 'datamodels/keyword_db'
FIXED_SCHEMA = 'fixed'
SCHEMA_TRUTH = 'truth/test_schema_editor'


@pytest.fixture(scope='module')
def keyword_db(jail, rtdata_module):
    """Define the keyword database"""
    rt = rtdata_module

    keyword_db = Path('keyword_db')
    keyword_db.mkdir()
    with pushdir(str(keyword_db)):
        schemas = rt.data_glob(KEYWORD_DB, glob='*.json')
        for schema in schemas:
            rt.get_data(schema)

    return keyword_db


@pytest.fixture(scope='module')
def run_editor_full(jail, keyword_db):
    """Just run the editor"""

    editor = se.Schema_editor(
        input=str(keyword_db), output=FIXED_SCHEMA,
        add=True, delete=True, edit=True, rename=True
    )
    editor.change()

    # The fixed schema live in a date-stamped folder under FIXED_SCHEMA.
    # It will be the only folder, so we can just return the first from the list.
    return os.path.join(FIXED_SCHEMA, os.listdir(FIXED_SCHEMA)[0])


def test_limit_datamodels():
    """Limit the datamodel list"""

    # Load the `DataModel` schemas as normal.
    models = se.Model_db()

    # Now reaload the models but exclude everything except the first model.
    exclude_list = list(models)[1:]
    models_limited = se.Model_db(exclude=exclude_list)
    assert len(models_limited) == 1


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'schema',
    ['container.schema.yaml', 'core.schema.yaml',
     'extract1dimage.schema.yaml',
     'guider_cal.schema.yaml', 'guider_raw.schema.yaml',
     'ifucube.schema.yaml',
     'keyword_exptype.schema.yaml', 'keyword_pband.schema.yaml', 'keyword_readpatt.schema.yaml',
     'lev3_prod.schema.yaml',
     'multiextract1d.schema.yaml', 'multispec.schema.yaml',
     'pathloss.schema.yaml',
     'referencefile.schema.yaml',
     'wcsinfo.schema.yaml',]
)
def test_full_run(jail, schema, run_editor_full, rtdata_module):
    """Check fixed schema files"""
    rt = rtdata_module
    rt.output = os.path.join(run_editor_full, schema)
    rt.get_truth(SCHEMA_TRUTH + '/' + schema)

    diff = list(text_diff(rt.output, rt.truth))
    assert not diff, sys.stderr.writelines(diff)


@pytest.mark.bigdata
def test_no_option_warning(jail, keyword_db):
    """If no operations is given, raise an error"""
    editor = se.Schema_editor(
        input=keyword_db
    )
    with pytest.raises(RuntimeError):
        editor.change()
