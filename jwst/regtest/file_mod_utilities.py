"""
Callback functions, which will be run by the Regression Test DaTa Updating tool (rtdtu).

These functions can have any inputs and outputs, and they must be wrapped in a function that has
no inputs or outputs.

Information about the tool and how it interacts with the jwst package can be found in the wiki at
https://github.com/spacetelescope/jwst/wiki/RTFile%3A-Setting-up-a-mapping-between-regression-test-data-and-MAST
"""

from pathlib import Path

from stdatamodels.jwst import datamodels

from jwst.lib.suffix import replace_suffix
from jwst.regtest.regtestdata import find_suffix


def trim_tso_data(file, ints_to_keep, intstart, ints_offset):
    """
    Trim TSO data and save into new file in the same directory as the input file.

    Parameters
    ----------
    file : str
        Full path and name of the file name to trim.
    ints_to_keep : int
        Number of integrations to keep.
    intstart : int
        Starting integration number.
    ints_offset : int
        Offset of integration number to start the trim in the array.
    """
    file_path = Path(file)
    root, suffix = find_suffix(file_path.name)
    modfname = replace_suffix(root, "mod_" + suffix) + ".fits"
    with datamodels.open(file) as dm:
        dm.meta.filename = modfname
        dm.meta.exposure.integration_start = intstart
        dm.meta.exposure.integration_end = intstart + ints_to_keep - 1
        dm.data = dm.data[ints_offset : ints_offset + ints_to_keep, ...]
        if dm.hasattr("refout"):
            dm.refout = dm.refout[ints_offset : ints_offset + ints_to_keep, ...]
        if dm.hasattr("error"):
            dm.error = dm.error[ints_offset : ints_offset + ints_to_keep, ...]
        if dm.hasattr("var_poisson"):
            dm.var_poisson = dm.var_poisson[ints_offset : ints_offset + ints_to_keep, ...]
        if dm.hasattr("var_noise"):
            dm.var_rnoise = dm.var_rnoise[ints_offset : ints_offset + ints_to_keep, ...]
        if dm.hasattr("int_times"):
            dm.int_times = dm.int_times[ints_offset : ints_offset + ints_to_keep]
        dm.save(modfname)
