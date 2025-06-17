"""Build fixtures for testing pipelines."""

import json
from pathlib import Path
import pytest
import tempfile

# Import from the common helpers module
# simply to make available from this module.
from jwst.tests.helpers import abspath  # noqa: F401

from jwst.associations import load_asn

SCRIPT_PATH = Path(__file__).parent
SCRIPT_DATA_PATH = SCRIPT_PATH / "data"


@pytest.fixture
def update_asn_basedir(asn_file, root=None):
    """
    Create an association with filenames update for a different directory.

    Parameters
    ----------
    asn_file : str
        The original association file
    root : str
        The root directory where the data actually resides

    Returns
    -------
    updated_asn_path : str
        The updated association file path
    """
    if root is None:
        root = Path.cwd()

    with Path(asn_file).open() as fd:
        asn = load_asn(fd)

    for product in asn["products"]:
        for member in product["members"]:
            expname = member["expname"]
            expname = Path(expname).name
            member["expname"] = root / expname

    _, updated_asn_path = tempfile.mkstemp(suffix=".json")
    with Path(updated_asn_path).open("w") as fd:
        json.dump(asn, fd)

    return updated_asn_path
