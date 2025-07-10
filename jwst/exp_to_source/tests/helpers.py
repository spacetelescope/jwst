"""Testing helpers."""

from astropy.utils.data import get_pkg_data_filename

# Calculate some extra constants
INPUT_FILES = [
    get_pkg_data_filename(f"data/jwst_nod{i}_cal.fits", package="jwst.exp_to_source.tests")
    for i in (1, 2, 3)
]
