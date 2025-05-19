import numpy as np
from numpy.testing import assert_allclose
from astropy.io import fits
import stdatamodels.jwst.datamodels as dm

from jwst.combine_1d._fileio import save_wfss_c1d


def test_save_wfss_c1d(tmp_cwd):
    """Test flat file format saving of WFSS c1d data."""

    # Set up a list of MultiSlitModel objects that look like outputs from combine_1d
    n_sources = 5
    n_rows = 3

    results_list = []
    for i in range(n_sources):
        multi = dm.MultiCombinedSpecModel()
        # Only one spec per model in this case
        spec = dm.CombinedSpecModel()
        spectable_dtype = spec.schema["properties"]["spec_table"]["datatype"]
        recarray_dtype = [(d["name"], d["datatype"]) for d in spectable_dtype]
        spec.source_id = n_sources - i # reverse the order to test sorting
        spec_table = np.recarray((n_rows,), dtype=recarray_dtype)
        spec_table["WAVELENGTH"] = [1.0, 2.0, 3.0]
        spec_table["FLUX"] = [1.0, 2.0, 3.0]
        spec.spec_table = spec_table
        multi.spec.append(spec)
        results_list.append(multi)

    # Save
    save_wfss_c1d(results_list, "test_c1d.fits")

    # Re-open as FITS and check the data
    hdul = fits.open("test_c1d.fits")
    assert len(hdul) == 3  # one primary HDU, one BinTableHDU, one ASDF extension
    bintable = hdul[1].data
    assert bintable.shape == (n_sources,)

    # ensure the source_id is the first column in the table
    assert_allclose(bintable["SOURCE_ID"], bintable.field(0))
    # Ensure the source_ids are sorted
    source_ids = bintable["SOURCE_ID"]
    assert np.all(np.diff(source_ids) >= 0)

    hdul.close()
