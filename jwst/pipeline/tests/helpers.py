"""test helpers"""

import json
import os
from os import path
import pytest
import tempfile

# Import from the common helpers module
# simply to make available from this module.
from ...tests.helpers import abspath # noqa: F401

from ...associations import load_asn

SCRIPT_PATH = path.dirname(__file__)
SCRIPT_DATA_PATH = path.join(SCRIPT_PATH, 'data')


@pytest.fixture
def update_asn_basedir(asn_file, root=None):
    """Create an association with filenames update
    for a different directory

    Parameters
    ----------
    asn_file: str
        The original association file

    root: str
        The root directory where the data actually resides

    Returns
    -------
    updated_asn_path: str
        The updated association file path
    """
    if root is None:
        root = os.getcwd()

    with open(asn_file) as fd:
        asn = load_asn(fd)

    for product in asn['products']:
        for member in product['members']:
            expname = member['expname']
            old_root, expname = path.split(expname)
            member['expname'] = path.join(root, expname)

    _, updated_asn_path = tempfile.mkstemp(suffix='.json')
    with open(updated_asn_path, 'w') as fd:
        json.dump(asn, fd)

    return updated_asn_path
