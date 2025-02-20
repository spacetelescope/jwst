# Testing helpers
from pathlib import Path


INPUT_FILES_GLOB = "data/jwst_nod?_cal.fits"


def t_path(partial_path):
    """Construct the full path for test files."""
    test_dir = Path(__file__)
    return Path(test_dir) / partial_path


# Calculate some extra constants
INPUT_FILES = t_path(INPUT_FILES_GLOB).glob()
