# Testing helpers
from pathlib import Path


def t_path(partial_path):
    """Construct the full path for test files."""
    test_dir = Path(__file__).parent
    return Path(test_dir) / partial_path


# Calculate some extra constants
INPUT_FILES = list(t_path("data").glob("jwst_nod?_cal.fits"))
